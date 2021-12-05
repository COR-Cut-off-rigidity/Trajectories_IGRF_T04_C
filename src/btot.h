//
// Created by pc on 7.4.2020.
//

#ifndef PARTICLE_TRAJECTORIES_CALCULATIONS_BTOT_H
#define PARTICLE_TRAJECTORIES_CALCULATIONS_BTOT_H

void btot(const float *X, const float *Y, const float *Z, const float *G, const float *H, const float *REC,
          const float *A1, const float *A2, const float *A3, const float *PARMOD, float *BX, float *BY, float *BZ);

void T04_s(const float *PARMOD, const float *PS, const float *X, const float *Y, const float *Z,
           float *BX, float *BY, float *BZ);

void externfunc(const float *PARMOD, const int *IOPGEN, const int *IOPT, const int *IOPB, const int *IOPR,
                const float *PS, const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ);

void shlcar3x3(const double *X, const double *Y, const double *Z, const float *PS, double *BX, double *BY, double *BZ);

void dipole(const float *PS, const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ);

void igrf_gsm(float XGSM, float YGSM, float ZGSM, const float *G, const float *H, const float *REC,
              const float *A1, const float *A2, const float *A3, float *HXGSM, float *HYGSM, float *HZGSM);

#endif //PARTICLE_TRAJECTORIES_CALCULATIONS_BTOT_H
