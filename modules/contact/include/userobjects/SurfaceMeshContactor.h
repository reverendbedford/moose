//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LevelSetContactor.h"
#include <memory>
#include <vector>

namespace libMesh
{
class ReplicatedMesh;
class Elem;
}

class KDTree;
class TriangleManifold;

/**
 * Rigid contactor described by a user-supplied closed, oriented, triangulated
 * surface mesh (typically STL).  Provides the same LevelSetContactor
 * `signedDistance` / `normal` interface as SphereContactor, computed as:
 *
 *   sign(x)   = -1 if x is inside the manifold, +1 outside
 *               (from TriangleManifold::contains)
 *   |g_LS(x)| = distance from x to the closest point on any triangle
 *               (KDTree nearest-centroid search + closest-point-on-triangle)
 *   n(x)      = grad(g_LS) = (x - closest_pt) / |x - closest_pt|
 *
 * The Hessian is inherited zero from LevelSetContactor: piecewise-flat facets
 * have zero true Hessian in facet interiors, so the finite-strain geometric
 * Jacobian term is dropped.
 */
class SurfaceMeshContactor : public LevelSetContactor
{
public:
  static InputParameters validParams();
  SurfaceMeshContactor(const InputParameters & parameters);
  virtual ~SurfaceMeshContactor();

  virtual void initialSetup() override;

  virtual Real signedDistance(const Point & x) const override;
  virtual RealVectorValue normal(const Point & x) const override;

private:
  /**
   * Find the closest point on `tri` (assumed Tri3) to `p`.  Standard
   * barycentric projection with clamping to edges/vertices.
   */
  static Point closestPointOnTriangle(const Point & p, const libMesh::Elem & tri);

  /**
   * KDTree candidate search + closest-point evaluation; returns the closest
   * point across all searched candidate triangles, and populates `closest_tri`
   * with the triangle it lies on.
   */
  Point closestSurfacePoint(const Point & x, const libMesh::Elem *& closest_tri) const;

  const FileName _file;
  const Point _translation;
  const Real _scale;
  const Real _surface_tolerance;
  const unsigned int _K;

  std::unique_ptr<libMesh::ReplicatedMesh> _mesh;
  std::unique_ptr<TriangleManifold> _manifold;
  std::unique_ptr<KDTree> _kd_tree;
  std::vector<Point> _centroids;
  std::vector<const libMesh::Elem *> _triangles;
};
