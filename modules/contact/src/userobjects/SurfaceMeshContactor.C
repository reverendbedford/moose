//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SurfaceMeshContactor.h"

#include "KDTree.h"
#include "TriangleManifold.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/replicated_mesh.h"

#include <algorithm>
#include <cmath>

registerMooseObject("ContactApp", SurfaceMeshContactor);

InputParameters
SurfaceMeshContactor::validParams()
{
  InputParameters params = LevelSetContactor::validParams();
  params.addClassDescription("Rigid contactor described by a closed, oriented, triangulated "
                             "surface mesh (e.g. STL).  Provides pointwise signed distance and "
                             "outward normal for use with RigidBodyNodalNCPKernel and "
                             "RigidBodyNormalMechanicalContact.");
  params.addRequiredParam<FileName>(
      "file",
      "Surface mesh file.  Format is dispatched from the extension via libMesh's mesh readers "
      "(e.g. .stl, .e).  Must be a closed, consistently outward-oriented 2-manifold of Tri3 "
      "elements after read+transform.");
  params.addParam<Point>("translation", Point(0, 0, 0), "Applied after scaling.");
  params.addRangeCheckedParam<Real>(
      "scale", 1.0, "scale > 0", "Uniform scale applied to the read mesh before translation.");
  params.addRangeCheckedParam<Real>(
      "surface_tolerance",
      1e-8,
      "surface_tolerance > 0",
      "Absolute tolerance used to validate the mesh (via TriangleManifold) at load time and "
      "to detect on-surface queries.  Choose relative to the mesh length scale.");
  return params;
}

SurfaceMeshContactor::SurfaceMeshContactor(const InputParameters & p)
  : LevelSetContactor(p),
    _file(getParam<FileName>("file")),
    _translation(getParam<Point>("translation")),
    _scale(getParam<Real>("scale")),
    _surface_tolerance(getParam<Real>("surface_tolerance"))
{
}

SurfaceMeshContactor::~SurfaceMeshContactor() = default;

void
SurfaceMeshContactor::initialSetup()
{
  // Serial (replicated) mesh so every rank has the full surface for pointwise queries.
  _mesh = std::make_unique<libMesh::ReplicatedMesh>(_communicator);
  _mesh->set_mesh_dimension(2);
  _mesh->read(_file);

  if (_scale != 1.0)
    libMesh::MeshTools::Modification::scale(*_mesh, _scale);
  if (_translation.norm() > 0.0)
    libMesh::MeshTools::Modification::translate(
        *_mesh, _translation(0), _translation(1), _translation(2));

  _mesh->prepare_for_use();

  // TriangleManifold's ctor validates Tri3-only + closed + consistently oriented.
  // We only need the validation (mooseError-on-invalid); we do not keep the
  // manifold instance because queryAt() derives the SDF sign from the closest
  // triangle's face normal, which is O(1) and does not need a runtime ray cast.
  TriangleManifold validate(*_mesh, _surface_tolerance);
  (void)validate;

  const auto n_active = _mesh->n_active_elem();
  _centroids.reserve(n_active);
  _triangles.reserve(n_active);
  for (const auto * elem : _mesh->active_element_ptr_range())
  {
    _centroids.push_back(elem->vertex_average());
    _triangles.push_back(elem);
  }

  _kd_tree = std::make_unique<KDTree>(_centroids, /*leaf_max_size=*/10);
}

RealVectorValue
SurfaceMeshContactor::faceNormal(const libMesh::Elem & tri)
{
  const Point & a = tri.point(0);
  const Point & b = tri.point(1);
  const Point & c = tri.point(2);
  RealVectorValue n = (b - a).cross(c - a);
  const Real nn = n.norm();
  if (nn > 0.0)
    n /= nn;
  return n;
}

Point
SurfaceMeshContactor::closestPointOnTriangle(const Point & p, const libMesh::Elem & tri)
{
  // Barycentric closest-point on a triangle: classical seven-region test with
  // clamping to edges/vertices.  Based on Ericson, "Real-Time Collision Detection".
  const Point & a = tri.point(0);
  const Point & b = tri.point(1);
  const Point & c = tri.point(2);

  const Point ab = b - a;
  const Point ac = c - a;
  const Point ap = p - a;

  const Real d1 = ab * ap;
  const Real d2 = ac * ap;
  if (d1 <= 0.0 && d2 <= 0.0)
    return a;

  const Point bp = p - b;
  const Real d3 = ab * bp;
  const Real d4 = ac * bp;
  if (d3 >= 0.0 && d4 <= d3)
    return b;

  const Real vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
  {
    const Real v = d1 / (d1 - d3);
    return a + v * ab;
  }

  const Point cp = p - c;
  const Real d5 = ab * cp;
  const Real d6 = ac * cp;
  if (d6 >= 0.0 && d5 <= d6)
    return c;

  const Real vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
  {
    const Real w = d2 / (d2 - d6);
    return a + w * ac;
  }

  const Real va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
  {
    const Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return b + w * (c - b);
  }

  // Interior of the triangle.
  const Real denom = 1.0 / (va + vb + vc);
  const Real v = vb * denom;
  const Real w = vc * denom;
  return a + ab * v + ac * w;
}

Point
SurfaceMeshContactor::closestSurfacePoint(const Point & x, const libMesh::Elem *& closest_tri) const
{
  std::vector<std::size_t> indices(1);
  _kd_tree->neighborSearch(x, 1, indices);

  const libMesh::Elem * hit = _triangles[indices.front()];
  Point best_cp = closestPointOnTriangle(x, *hit);
  Real best_d2 = (x - best_cp).norm_sq();
  closest_tri = hit;

  // Robustness sweep: also test the hit's edge neighbors.  Handles the case
  // where the true closest triangle is not the one whose centroid is nearest
  // (possible when triangles have very different sizes or when x sits closer
  // to a neighboring triangle's edge than to hit's).  Up to 3 extra
  // closest-point tests for Tri3 — cheap relative to the KDTree lookup.
  for (const auto s : make_range(hit->n_sides()))
  {
    const libMesh::Elem * neigh = hit->neighbor_ptr(s);
    if (!neigh)
      continue; // safety; a validated closed manifold has no null neighbors.
    const Point cp = closestPointOnTriangle(x, *neigh);
    const Real d2 = (x - cp).norm_sq();
    if (d2 < best_d2)
    {
      best_d2 = d2;
      best_cp = cp;
      closest_tri = neigh;
    }
  }
  return best_cp;
}

LevelSetContactor::Query
SurfaceMeshContactor::queryAtRaw(const Point & x) const
{
  // One KDTree search + neighbor sweep serves gap, normal, and hessian.  The
  // per-quantity accessors below defer to this method, so there is no slow
  // path.
  const libMesh::Elem * tri = nullptr;
  const Point cp = closestSurfacePoint(x, tri);
  const RealVectorValue v = x - cp;
  const Real d = v.norm();

  // Sign of g_LS: derived from the closest triangle's face normal rather than
  // a global point-in-solid ray cast.  For any consistently outward-oriented
  // closed manifold, v is essentially parallel to the outward normal of the
  // closest surface feature, so sign(v · n_face) tells us which side we're
  // on.  This is O(1) per query and avoids the discrete sign flip that a
  // TriangleManifold::contains fallback would introduce on the medial axis.
  const RealVectorValue n_face = faceNormal(*tri);
  const Real vdotn = v * n_face;
  const Real sign = vdotn >= 0.0 ? 1.0 : -1.0;

  Query q;
  q.gap = sign * d;

  if (d > _surface_tolerance)
    // grad(g_LS(x)) = sign(x) * (x - CP(x)) / |x - CP(x)|.  Always the outward
    // surface normal at CP: v/d itself on the outside, flipped on the inside.
    q.normal = (sign / d) * v;
  else
    // On-surface: v ≈ 0, use the closest triangle's face normal directly.
    q.normal = n_face;

  // Piecewise-flat facets ⇒ true Hessian is zero on facet interiors.
  q.hessian = RealTensorValue();
  return q;
}

Real
SurfaceMeshContactor::signedDistanceRaw(const Point & x) const
{
  return queryAtRaw(x).gap;
}

RealVectorValue
SurfaceMeshContactor::normalRaw(const Point & x) const
{
  return queryAtRaw(x).normal;
}
