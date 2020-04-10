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
/// Linear elasticity for solid mechanics example using PETSc

#ifndef INCOMPLIN_ELAS_H
#define INCOMPLIN_ELAS_H

#ifndef __CUDACC__
#  include <math.h>
#endif

#ifndef PHYSICS_STRUCT
#define PHYSICS_STRUCT
typedef struct Physics_private *Physics;
struct Physics_private {
  CeedScalar   nu;      // Poisson's ratio
  CeedScalar   E;       // Young's Modulus
};
#endif

// -----------------------------------------------------------------------------
// Residual evaluation for linear elasticity
// -----------------------------------------------------------------------------
CEED_QFUNCTION(IncompLinElasF)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                         CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*ug)[3][CEED_Q_VLA] = (const CeedScalar(*)[3][CEED_Q_VLA])in[0],
                   (*qdata)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[1];

  // Outputs
  CeedScalar (*dvdX)[3][CEED_Q_VLA] = (CeedScalar(*)[3][CEED_Q_VLA])out[0];
             // gradu not used for linear elasticity
             // (*gradu)[3][CEED_Q_VLA] = (CeedScalar(*)[3][CEED_Q_VLA])out[1];
  // *INDENT-ON*

  // Context
  const Physics context = (Physics)ctx;
  const CeedScalar E  = context->E;
  const CeedScalar nu = context->nu;
  const CeedScalar mu = E/2*(1.+2*nu);

  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    // Read spatial derivatives of u
    // *INDENT-OFF*
    const CeedScalar du[3][3]   = {{ug[0][0][i],
                                    ug[1][0][i],
                                    ug[2][0][i]},
                                   {ug[0][1][i],
                                    ug[1][1][i],
                                    ug[2][1][i]},
                                   {ug[0][2][i],
                                    ug[1][2][i],
                                    ug[2][2][i]}
                                  };
    // -- Qdata
    const CeedScalar wJ         =   qdata[0][i];
    const CeedScalar dXdx[3][3] = {{qdata[1][i],
                                    qdata[2][i],
                                    qdata[3][i]},
                                   {qdata[4][i],
                                    qdata[5][i],
                                    qdata[6][i]},
                                   {qdata[7][i],
                                    qdata[8][i],
                                    qdata[9][i]}
                                  };
    // *INDENT-ON*

    // Compute gradu
    //   dXdx = (dx/dX)^(-1)
    // Apply dXdx to du = gradu
    CeedScalar gradu[3][3];
    for (CeedInt j = 0; j < 3; j++)     // Component
      for (CeedInt k = 0; k < 3; k++) { // Derivative
        gradu[j][k] = 0;
        for (CeedInt m = 0; m < 3; m++)
          gradu[j][k] += dXdx[m][k] * du[j][m];
      }

    // Compute Strain : e (epsilon)
    // e = 1/2 (grad u + (grad u)^T)

    // *INDENT-OFF*
    const CeedScalar e[3][3] = {{(gradu[0][0] + gradu[0][0])/2.,
                                 (gradu[0][1] + gradu[1][0])/2.,
                                 (gradu[0][2] + gradu[2][0])/2.},
                                {(gradu[1][0] + gradu[0][1])/2.,
                                 (gradu[1][1] + gradu[1][1])/2.,
                                 (gradu[1][2] + gradu[2][1])/2.},
                                {(gradu[2][0] + gradu[0][2])/2.,
                                 (gradu[2][1] + gradu[1][2])/2.,
                                 (gradu[2][2] + gradu[2][2])/2.}
                               };
    // *INDENT-ON*

    //
    // Formulation uses Voigt notation:
    //   Deviatoric
    //  stress (sigma)      strain (epsilon)
    //
    //    [sigmaDev00]             [e00]
    //    [sigmaDev11]             [e11]
    //    [sigmaDev22]   =  S   *  [e22]
    //    [sigmaDev12]             [e12]
    //    [sigmaDev02]             [e02]
    //    [sigmaDev01]             [e01]
    //
    //        where
    //          [  4/3  -1/3    -1/3          ]
    //          [ -1/3   4/3    -1/3          ]
    //          [ -1/3  -1/3     4/3          ]
    // S = mu * [                     1       ]
    //          [                        1    ]
    //          [                           1 ]
    //  with mu = E/2*(1+2*nu)

    // Above Voigt Notation is placed in a 3x3 matrix:
    const CeedScalar sigmaDev00 = mu*( 4./3*e[0][0] - 1./3*e[1][1] - 1./3*e[2][2]),
                     sigmaDev11 = mu*(-1./3*e[0][0] + 4./3*e[1][1] - 1./3*e[2][2]),
                     sigmaDev22 = mu*(-1./3*e[0][0] - 1./3*e[1][1] + 4./3*e[2][2]),
                     sigmaDev12 = mu*e[1][2],
                     sigmaDev02 = mu*e[0][2],
                     sigmaDev01 = mu*e[0][1];
    // *INDENT-OFF*
    const CeedScalar sigmaDev[3][3] = {{sigmaDev00, sigmaDev01, sigmaDev02},
                                       {sigmaDev01, sigmaDev11, sigmaDev12},
                                       {sigmaDev02, sigmaDev12, sigmaDev22}
                                      };
    // *INDENT-ON*

    // Apply dXdx^T and weight to sigma
    for (CeedInt j = 0; j < 3; j++)     // Component
      for (CeedInt k = 0; k < 3; k++) { // Derivative
        dvdX[k][j][i] = 0;
        for (CeedInt m = 0; m < 3; m++)
          dvdX[k][j][i] += dXdx[k][m] * sigmaDev[j][m] * wJ;
      }

  } // End of Quadrature Point Loop

  return 0;
}

// -----------------------------------------------------------------------------
// Jacobian evaluation for linear elasticity
// -----------------------------------------------------------------------------
CEED_QFUNCTION(IncompLinElasdF)(void *ctx, CeedInt Q, const CeedScalar *const *in,
                          CeedScalar *const *out) {
  // *INDENT-OFF*
  // Inputs
  const CeedScalar (*deltaug)[3][CEED_Q_VLA] = (const CeedScalar(*)[3][CEED_Q_VLA])in[0],
                   (*qdata)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[1];
                   // gradu not used for linear elasticity
                   // (*gradu)[3][Q] = (CeedScalar(*)[3][Q])in[2];

  // Outputs
  CeedScalar (*deltadvdX)[3][CEED_Q_VLA] = (CeedScalar(*)[3][CEED_Q_VLA])out[0];
  // *INDENT-ON*

  // Context
  const Physics context = (Physics)ctx;
  const CeedScalar E  = context->E;
  const CeedScalar nu = context->nu;
  const CeedScalar mu = E/2*(1.+2*nu);

  // Quadrature Point Loop
  CeedPragmaSIMD
  for (CeedInt i=0; i<Q; i++) {
    // Read spatial derivatives of u
    // *INDENT-OFF*
    const CeedScalar deltadu[3][3] = {{deltaug[0][0][i],
                                       deltaug[1][0][i],
                                       deltaug[2][0][i]},
                                      {deltaug[0][1][i],
                                       deltaug[1][1][i],
                                       deltaug[2][1][i]},
                                      {deltaug[0][2][i],
                                       deltaug[1][2][i],
                                       deltaug[2][2][i]}
                                     };
    // -- Qdata
    const CeedScalar wJ         =      qdata[0][i];
    const CeedScalar dXdx[3][3] =    {{qdata[1][i],
                                       qdata[2][i],
                                       qdata[3][i]},
                                      {qdata[4][i],
                                       qdata[5][i],
                                       qdata[6][i]},
                                      {qdata[7][i],
                                       qdata[8][i],
                                       qdata[9][i]}
                                     };
    // *INDENT-ON*

    // Compute graddeltau
    //   dXdx = (dx/dX)^(-1)
    // Apply dXdx to deltadu = graddeltau
    CeedScalar graddeltau[3][3];
    for (CeedInt j = 0; j < 3; j++)     // Component
      for (CeedInt k = 0; k < 3; k++) { // Derivative
        graddeltau[j][k] = 0;
        for (CeedInt m = 0; m < 3; m++)
          graddeltau[j][k] += dXdx[m][k] * deltadu[j][m];
      }

    // Compute Strain : e (epsilon)
    // e = 1/2 (grad u + (grad u)^T)
    // *INDENT-OFF*
    const CeedScalar de[3][3] = {{(graddeltau[0][0] + graddeltau[0][0])/2.,
                                  (graddeltau[0][1] + graddeltau[1][0])/2.,
                                  (graddeltau[0][2] + graddeltau[2][0])/2.},
                                 {(graddeltau[1][0] + graddeltau[0][1])/2.,
                                  (graddeltau[1][1] + graddeltau[1][1])/2.,
                                  (graddeltau[1][2] + graddeltau[2][1])/2.},
                                 {(graddeltau[2][0] + graddeltau[0][2])/2.,
                                  (graddeltau[2][1] + graddeltau[1][2])/2.,
                                  (graddeltau[2][2] + graddeltau[2][2])/2.}
                                };

    // *INDENT-ON*
    // Formulation uses Voigt notation:
    //    Deviatoric
    //  stress (sigma)      strain (epsilon)
    //
    //    [dsigmaDev00]             [de00]
    //    [dsigmaDev11]             [de11]
    //    [dsigmaDev22]   =  S   *  [de22]
    //    [dsigmaDev12]             [de12]
    //    [dsigmaDev02]             [de02]
    //    [dsigmaDev01]             [de01]
    //
    //        where
    //          [  4/3  -1/3    -1/3          ]
    //          [ -1/3   4/3    -1/3          ]
    //          [ -1/3  -1/3     4/3          ]
    // S = mu * [                     1       ]
    //          [                        1    ]
    //          [                           1 ]
    //  with mu = E/2*(1+2*nu)

    // Above Voigt Notation is placed in a 3x3 matrix:
    const CeedScalar dsigmaDev00 = mu*( 4./3*de[0][0] - 1./3*de[1][1] - 1./3*de[2][2]),
                     dsigmaDev11 = mu*(-1./3*de[0][0] + 4./3*de[1][1] - 1./3*de[2][2]),
                     dsigmaDev22 = mu*(-1./3*de[0][0] - 1./3*de[1][1] + 4./3*de[2][2]),
                     dsigmaDev12 = mu*de[1][2],
                     dsigmaDev02 = mu*de[0][2],
                     dsigmaDev01 = mu*de[0][1];

    // *INDENT-OFF*
    const CeedScalar dsigmaDev[3][3] = {{dsigmaDev00, dsigmaDev01, dsigmaDev02},
                                        {dsigmaDev01, dsigmaDev11, dsigmaDev12},
                                        {dsigmaDev02, dsigmaDev12, dsigmaDev22}
                                       };
    // *INDENT-ON*

    // Apply dXdx^T and weight
    for (CeedInt j = 0; j < 3; j++)     // Component
      for (CeedInt k = 0; k < 3 ; k++) { // Derivative
        deltadvdX[k][j][i] = 0;
        for (CeedInt m = 0; m < 3; m++)
          deltadvdX[k][j][i] += dXdx[k][m] * dsigmaDev[j][m] * wJ;
      }

  } // End of Quadrature Point Loop

  return 0;
}
// -----------------------------------------------------------------------------
#endif //End of LIN_ELAS_H
