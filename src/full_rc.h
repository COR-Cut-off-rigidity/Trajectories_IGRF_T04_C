//
// Created by pc on 7.4.2020.
//

#ifndef PARTICLE_TRAJECTORIES_CALCULATIONS_FULL_RC_H
#define PARTICLE_TRAJECTORIES_CALCULATIONS_FULL_RC_H

void full_rc(const int *IOPR, const float *PS, const double *X, const double *Y, const double *Z, const double *PHI,
             const double *SC_SY, const double *SC_PR, double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC,
             double *BYPRC, double *BZPRC);

void src_prc(const int *IOPR, const double *SC_SY, const double *SC_PR, const double *PHI, const float *PS,
             const double *X, const double *Y, const double *Z, double *BXSRC, double *BYSRC, double *BZSRC,
             double *BXPRC, double *BYPRC, double *BZPRC);

void rc_symm(const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ);

double ap(const double *R, const double *SINT, const double *COST);

void prc_symm(const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ);

double apprc(const double *R, const double *SINT, const double *COST);

void prc_quad(const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ);

double br_prc_q(const double *R, const double *SINT, const double *COST);

double bt_prc_q(const double *R, const double *SINT, const double *COST);

void calc_ffs(const double *A, const double *A0, const double *DA, double *F, double *FA, double *FS);

void rc_shield(const double *A, const float *PS, const double *X_SC, const double *X, const double *Y, const double *Z,
               double *BX, double *BY, double *BZ);

#endif //PARTICLE_TRAJECTORIES_CALCULATIONS_FULL_RC_H
