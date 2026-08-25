//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "UzawaSolveObject.h"

#include "Executioner.h"
#include "FEProblemBase.h"
#include "MooseVariableScalar.h"
#include "NonlinearSystemBase.h"
#include "RigidBodyLoadControl.h"
#include "SystemBase.h"

#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_nonlinear_solver.h"

#include "MooseApp.h"
#include "OutputWarehouse.h"
#include "PetscOutput.h"

#include <petscksp.h>
#include <petscsnes.h>

#include <cmath>

UzawaSolveObject::UzawaSolveObject(Executioner & ex,
                                   unsigned int outer_max_iter,
                                   Real outer_abs_tol,
                                   Real outer_rel_tol,
                                   Real max_step,
                                   unsigned int damp_max_retries,
                                   bool outer_verbose,
                                   bool inner_verbose)
  : SolveObject(ex),
    _load_control(nullptr),
    _outer_max_iter(outer_max_iter),
    _outer_abs_tol(outer_abs_tol),
    _outer_rel_tol(outer_rel_tol),
    _max_step(max_step),
    _damp_max_retries(damp_max_retries),
    _outer_verbose(outer_verbose),
    _inner_verbose(inner_verbose)
{
}

Real
UzawaSolveObject::readScalarValue() const
{
  // The scalar variable driven by _load_control is `variable` on that
  // kernel.  Reach it via the kernel's own `variable()` accessor.  The
  // scalar variable's `sln()` returns a VariableValue (size 1 for
  // order=FIRST).
  const auto & scalar_var = _load_control->variable();
  return scalar_var.sln()[0];
}

void
UzawaSolveObject::writeScalarValue(Real new_value)
{
  const auto & scalar_var = _load_control->variable();
  // Scalar DoFs live on the (usually last) rank that owns them.  Guard
  // the set() with a local-DoF check so non-owning ranks no-op cleanly
  // in parallel.
  auto & sys = _problem.getNonlinearSystemBase(0);
  auto & soln = sys.solution();
  const auto & dof_indices = scalar_var.dofIndices();
  const auto & dof_map = sys.dofMap();
  const auto first_local = dof_map.first_dof(_problem.mesh().comm().rank());
  const auto end_local = dof_map.end_dof(_problem.mesh().comm().rank());
  for (const auto d : dof_indices)
    if (d >= first_local && d < end_local)
      soln.set(d, new_value);
  soln.close();
  // Refresh MOOSE's cached variable values from the modified solution.
  sys.update();
}

bool
UzawaSolveObject::solve()
{
  // Pass-through when no load-control kernel is present.
  if (!_load_control)
    return _inner_solve->solve();

  Real r_s_initial = -1.0;
  Real s_current = readScalarValue();

  for (unsigned int outer = 0; outer < _outer_max_iter; ++outer)
  {
    // ------------------------------------------------------------------
    // (a) Toggle to PinScalar mode with s_pin = current s.  Primal solve
    //     will then hold s at s_current and only move (u, lambda).
    // ------------------------------------------------------------------
    _load_control->setMode(RigidBodyLoadControl::Mode::PinScalar, s_current);

    // ------------------------------------------------------------------
    // (b) Primal solve.  Optionally damp/retry on failure.  Silence
    //     the inner primal SNES/KSP monitor output unless
    //     `_inner_verbose` was set; the inner solve fires many times
    //     per outer iter and its per-iter `Nonlinear |R|` /
    //     `Linear |R|` lines drown out the outer Uzawa progress.
    // ------------------------------------------------------------------
    if (!_inner_verbose)
      silenceInnerSolveMonitors();
    bool primal_ok = _inner_solve->solve();
    if (!_inner_verbose)
      restoreInnerSolveMonitors();
    if (!primal_ok && _damp_max_retries > 0)
    {
      // If the primal solve failed, we can't rescue it from here (the
      // damping is on the outer scalar step, not the inner) -- report
      // failure and let the transient time-stepper cut back dt.
      if (_outer_verbose)
        _console << "Uzawa: primal solve failed on outer iter " << outer
                 << ", giving up so outer transient can cut back dt" << std::endl;
      _load_control->setMode(RigidBodyLoadControl::Mode::ForceBalance, 0.0);
      return false;
    }
    if (!primal_ok)
    {
      _load_control->setMode(RigidBodyLoadControl::Mode::ForceBalance, 0.0);
      return false;
    }

    // ------------------------------------------------------------------
    // (c) Read R_s at the converged primal state.  The most recent
    //     residual assembly (inside SNES's convergence check) populated
    //     `_cached_reaction_minus_F` in PinScalar mode -- the reaction
    //     sum is computed identically in both modes.
    // ------------------------------------------------------------------
    const Real r_s = _load_control->currentReactionMinusF();
    if (outer == 0)
      r_s_initial = std::abs(r_s);

    if (_outer_verbose)
      _console << "Uzawa outer iter " << outer << ": s = " << s_current
               << ", |R_s| = " << std::abs(r_s) << std::endl;

    // ------------------------------------------------------------------
    // (d) Convergence check on |R_s|.
    // ------------------------------------------------------------------
    if (std::abs(r_s) < _outer_abs_tol ||
        (r_s_initial > 0.0 && std::abs(r_s) < _outer_rel_tol * r_s_initial))
    {
      _load_control->setMode(RigidBodyLoadControl::Mode::ForceBalance, 0.0);
      if (_outer_verbose)
        _console << "Uzawa converged in " << (outer + 1) << " outer iters" << std::endl;
      return true;
    }

    // ------------------------------------------------------------------
    // (e) Scalar Newton step: ds = -R_s / kss_stiffness, clipped to
    //     +/- max_step (a trust-region bound so a wildly-off
    //     kss_stiffness cannot make s explode).  Use the physically-
    //     motivated Kss = kss_stiffness * total_nodal_area *
    //     sign(direction . axis_hat) already assembled in ForceBalance
    //     mode.  Since we cannot cheaply extract that from a
    //     PinScalar assembly, approximate as kss_stiffness directly --
    //     off by an O(1) area factor, but the trust region and re-
    //     iteration absorb the mismatch.
    // ------------------------------------------------------------------
    // signedKssApprox() = kss_stiffness * sign(direction . axis_hat), so
    // it carries the correct sign for `ds = -R_s / signedKss` to move s
    // in the reaction-reducing direction on either axis orientation.
    const Real signed_kss = _load_control->signedKssApprox();
    if (signed_kss == 0.0)
    {
      mooseError("UzawaSolveObject: kss_stiffness on load-control kernel '",
                 _load_control->name(),
                 "' must be > 0 for the outer scalar Newton step to make "
                 "sense.  Set kss_stiffness to (approximately) the "
                 "deformable body's Young's modulus.");
    }

    Real ds = -r_s / signed_kss;
    if (std::abs(ds) > _max_step)
      ds = (ds > 0.0 ? _max_step : -_max_step);

    s_current += ds;
    writeScalarValue(s_current);
  }

  // Fell off the end without meeting the tolerance -- outer non-
  // convergence.  Restore ForceBalance so the outer transient's own
  // post-processing sees the physical residual.
  _load_control->setMode(RigidBodyLoadControl::Mode::ForceBalance, 0.0);
  _console << "Uzawa: outer loop did not converge in " << _outer_max_iter
           << " iters, propagating failure so transient can cut back dt"
           << std::endl;
  return false;
}

// A no-op SNESSetUpdate callback: PETSc invokes this at the start of
// each SNES iteration, BEFORE the SNES monitors fire.  Inside it we
// cancel the SNES and KSP monitors that MOOSE's PetscOutput installed
// via its `solveSetup` earlier in the same inner-solve entry (which
// runs unconditionally on every FEProblemBase::solve()).  Doing the
// cancel here rather than pre-solve dodges the reinstall race:
// MOOSE calls PetscOutput::solveSetup BEFORE any SNES iterations
// start, so a cancel pre-solve is undone before we reach iteration
// 0's monitor.  Cancelling from SetUpdate at step 0 short-circuits
// PETSc's own monitor pass for step 0 too.
static PetscErrorCode
uzawaSilenceMonitorsUpdate(SNES snes, PetscInt /*step*/)
{
  PetscFunctionBegin;
  KSP ksp = nullptr;
  PetscCall(SNESGetKSP(snes, &ksp));
  PetscCall(SNESMonitorCancel(snes));
  PetscCall(KSPMonitorCancel(ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

void
UzawaSolveObject::silenceInnerSolveMonitors()
{
  // MOOSE's `PetscOutputInterface` re-installs the SNES/KSP monitors
  // on every `FEProblemBase::solve()` (via
  // `initPetscOutputAndSomeSolverSettings -> OutputWarehouse::solveSetup`),
  // which happens INSIDE `_inner_solve->solve()` -- so a cancel here
  // is undone before Newton even starts.  Instead, register a
  // per-iteration `SNESSetUpdate` callback that cancels monitors on
  // every SNES iter (including step 0).  `restoreInnerSolveMonitors`
  // clears the update callback after the inner solve.
  auto & nl = _problem.getNonlinearSystemBase(/*sys=*/0);
  SNES snes = nl.getSNES();
  auto ierr = SNESSetUpdate(snes, uzawaSilenceMonitorsUpdate);
  LibmeshPetscCallA(_problem.comm().get(), ierr);
}

void
UzawaSolveObject::restoreInnerSolveMonitors()
{
  // Remove the update callback so subsequent (non-Uzawa-inner) solves
  // do not accidentally lose their monitors.  Passing nullptr as the
  // callback clears any previously-registered SetUpdate.
  auto & nl = _problem.getNonlinearSystemBase(/*sys=*/0);
  SNES snes = nl.getSNES();
  auto ierr = SNESSetUpdate(snes, nullptr);
  LibmeshPetscCallA(_problem.comm().get(), ierr);
}
