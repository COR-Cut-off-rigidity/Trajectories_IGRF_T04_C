//
// Created by pc on 7.4.2020.
//

#ifndef PARTICLE_TRAJECTORIES_CALCULATIONS_DEFORMED_H
#define PARTICLE_TRAJECTORIES_CALCULATIONS_DEFORMED_H

void deformed(const int *IOPT, float PS, const double *X, const double *Y, const double *Z,
              const double *DXSHIFT1, const double *DXSHIFT2, const double *D, const float *DELTADY,
              double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2);

void warped(const int *IOPT, const float *PS, const double *X, const double *Y, const double *Z,
            const double *DXSHIFT1, const double *DXSHIFT2, const double *D, const float *DELTADY,
            double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2);

void unwarped(const int *IOPT, const double *X, const double *Y, const double *Z,
              const double *DXSHIFT1, const double *DXSHIFT2, const double *D0, const float *DELTADY,
              double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2);

void taildisk(const double *D0, const double *DELTADX, const float *DELTADY,
              const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ);

void shlcar5x5(const float *A, const double *X, const double *Y, const double *Z, const double *DSHIFT,
               double *HX, double *HY, double *HZ);

#endif //PARTICLE_TRAJECTORIES_CALCULATIONS_DEFORMED_H
