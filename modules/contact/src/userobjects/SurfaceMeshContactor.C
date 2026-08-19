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
      "Absolute tolerance used by TriangleManifold for validation and near-surface "
      "classification; choose relative to the mesh length scale.");
  params.addRangeCheckedParam<unsigned int>(
      "nearest_neighbors",
      1,
      "nearest_neighbors >= 1",
      "Number of KDTree candidates to check per query.  1 is right for the vast majority of "
      "meshes; increase for meshes with highly non-uniform triangle sizes.");
  return params;
}

SurfaceMeshContactor::SurfaceMeshContactor(const InputParameters & p)
  : LevelSetContactor(p),
    _file(getParam<FileName>("file")),
    _translation(getParam<Point>("translation")),
    _scale(getParam<Real>("scale")),
    _surface_tolerance(getParam<Real>("surface_tolerance")),
    _K(getParam<unsigned int>("nearest_neighbors"))
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

  // TriangleManifold's ctor validates Tri3-only + closed + oriented and builds the inside/outside
  // acceleration structure.
  _manifold = std::make_unique<TriangleManifold>(*_mesh, _surface_tolerance);

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
  std::vector<std::size_t> indices(_K);
  _kd_tree->neighborSearch(x, _K, indices);

  Real best_d2 = std::numeric_limits<Real>::max();
  Point best_cp;
  closest_tri = nullptr;
  for (const auto idx : indices)
  {
    const auto * tri = _triangles[idx];
    const Point cp = closestPointOnTriangle(x, *tri);
    const Real d2 = (x - cp).norm_sq();
    if (d2 < best_d2)
    {
      best_d2 = d2;
      best_cp = cp;
      closest_tri = tri;
    }
  }
  return best_cp;
}

Real
SurfaceMeshContactor::signedDistance(const Point & x) const
{
  const libMesh::Elem * tri = nullptr;
  const Point cp = closestSurfacePoint(x, tri);
  const Real d = (x - cp).norm();
  return _manifold->contains(x) ? -d : d;
}

RealVectorValue
SurfaceMeshContactor::normal(const Point & x) const
{
  const libMesh::Elem * tri = nullptr;
  const Point cp = closestSurfacePoint(x, tri);
  const RealVectorValue v = x - cp;
  const Real nrm = v.norm();
  if (nrm > _surface_tolerance)
  {
    // grad(g_LS(x)) = sign(x) * (x - CP(x)) / |x - CP(x)| — always points
    // from interior toward exterior (the outward surface normal).  For an
    // outside query, x - CP already points outward; for an inside query, we
    // must flip it.
    const Real sign = _manifold->contains(x) ? -1.0 : 1.0;
    return (sign / nrm) * v;
  }

  // On-surface fallback: use the closest triangle's outward face normal.
  const Point & a = tri->point(0);
  const Point & b = tri->point(1);
  const Point & c = tri->point(2);
  RealVectorValue n = (b - a).cross(c - a);
  const Real nn = n.norm();
  if (nn > 0.0)
    n /= nn;
  return n;
}
