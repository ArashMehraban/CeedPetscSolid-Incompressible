// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.

/// @file
/// Boundary condition functions for solid mechanics example using PETSc

#include "../elasticity.h"

// -----------------------------------------------------------------------------
// Boundary Functions
// -----------------------------------------------------------------------------
// Note: If additional boundary conditions are added, an update is needed in
//         elasticity.h for the boundaryOptions variable.

// BCMMS - boundary function
// Values on all points of the mesh is set based on given solution below
// for u[0], u[1], u[2]
PetscErrorCode BCMMS(PetscInt dim, PetscReal loadIncrement,
                     const PetscReal coords[], PetscInt ncompu,
                     PetscScalar *u, void *ctx) {
  PetscScalar x = coords[0];
  PetscScalar y = coords[1];
  PetscScalar z = coords[2];

  PetscFunctionBeginUser;

  u[0] = exp(2*x)*sin(3*y)*cos(4*z) / 1e8 * loadIncrement;
  u[1] = exp(3*y)*sin(4*z)*cos(2*x) / 1e8 * loadIncrement;
  u[2] = exp(4*z)*sin(2*x)*cos(3*y) / 1e8 * loadIncrement;

  PetscFunctionReturn(0);
};

// BCZero - fix boundary values at zero
PetscErrorCode BCZero(PetscInt dim, PetscReal loadIncrement,
                      const PetscReal coords[], PetscInt ncompu,
                      PetscScalar *u, void *ctx) {
  PetscFunctionBeginUser;

  u[0] = 0;
  u[1] = 0;
  u[2] = 0;

  PetscFunctionReturn(0);
};

// BCClampTranslate - translate boundary values at fraction of load increment
PetscErrorCode BCClampTranslate(PetscInt dim, PetscReal loadIncrement,
                                const PetscReal coords[], PetscInt ncompu,
                                PetscScalar *u, void *ctx) {
  PetscScalar (*clampMax) = (PetscScalar(*))ctx;

  PetscFunctionBeginUser;

  u[0] = clampMax[0]*loadIncrement;
  u[1] = clampMax[1]*loadIncrement;
  u[2] = clampMax[2]*loadIncrement;

  PetscFunctionReturn(0);
};

#ifndef M_PI
#  define M_PI    3.14159265358979323846
#endif

// BCClampRotate - rotate boundary values at fraction of load increment
PetscErrorCode BCClampRotate(PetscInt dim, PetscReal loadIncrement,
                             const PetscReal coords[], PetscInt ncompu,
                             PetscScalar *u, void *ctx) {
  PetscScalar x = coords[0];
  PetscScalar y = coords[1];
  PetscScalar z = coords[2];
  PetscScalar (*clampMax) = (PetscScalar(*))ctx;

  PetscFunctionBeginUser;

  PetscScalar theta = clampMax[3]*M_PI*loadIncrement,
              kx = clampMax[0], ky = clampMax[1], kz = clampMax[2];
  PetscScalar c = cos(theta), s = sin(theta);

  u[0] = s*(-kz*y + ky*z) + (1-c)*(-ky*ky+kz*kz*x + kx*ky*y + kx*kz*z);
  u[1] = s*(kz*x + -kx*z) + (1-c)*(kx*ky*x - (kx*kx+kz*kz)*y + ky*kz*z);
  u[2] = s*(-ky*x + kx*y) + (1-c)*(kx*kz*x + ky*kz*y - (kx*kx+ky*ky)*z);

  PetscFunctionReturn(0);
};
