//
// Created by pc on 7.4.2020.
//

#ifndef PARTICLE_TRAJECTORIES_CALCULATIONS_BIRK_TOT_H
#define PARTICLE_TRAJECTORIES_CALCULATIONS_BIRK_TOT_H

void birk_tot(const int *IOPB, const float *PS,
              const double *X, const double *Y, const double *Z, const double *XKAPPA1, const double *XKAPPA2,
              double *BX11, double *BY11, double *BZ11, double *BX12, double *BY12, double *BZ12,
              double *BX21, double *BY21, double *BZ21, double *BX22, double *BY22, double *BZ22);

void birk_1n2(const int *NUMB, const int *MODE, const float *PS, const double *X, const double *Y, const double *Z,
              const double *XKAPPA, double *BX, double *BY, double *BZ);

void twocones(const float *A, const double *X, const double *Y, const double *Z, const int *MODE, const double *DTHETA,
              double *BX, double *BY, double *BZ);

void one_cone(const float *A, const double *X, const double *Y, const double *Z, const int *MODE, const double *DTHETA,
              double *BX, double *BY, double *BZ);

double r_s(const float *A, const double *R, const double *THETA);

double theta_s(const float *A, const double *R, const double *THETA);

void fialcos(const double *R, const double *THETA, const double *PHI, const int *MODE, const double *THETA0,
             const double *DT, double *BTHETA, double *BPHI);

void birk_shl(const float *A, const float *PS, const double *X_SC, const double *X, const double *Y, const double *Z,
              double *BX, double *BY, double *BZ);

#endif //PARTICLE_TRAJECTORIES_CALCULATIONS_BIRK_TOT_H
