//
// Created by pc on 7.4.2020.
//

#ifndef PARTICLE_TRAJECTORIES_CALCULATIONS_RECALC_H
#define PARTICLE_TRAJECTORIES_CALCULATIONS_RECALC_H

void recalc(int igrf_coefs_version, int iyear, const int *iday, const int *ihour, const int *min, const int *isec,
            float *G, float *H, float *REC, float *A1, float *A2, float *A3);

void interpolateOrExtrapolate(int igrf_coefs_version, int *year, const int *iday, float *G, float *H);

void interpolate(float *G, float *H, const float *Gx, const float *Hx, const float *Gy, const float *Hy, int edgeYear, int year, int iday);

void extrapolate(float *G, float *H, const float *Gx, const float *Hx, const float *Gy, const float *Hy, int maxYear, int year, int iday);

void sun(int *iyear, const int *iday, const int *ihour, const int *min, const int *isec,
         float *GST, float *SLONG, float *SRASN, float *SDEC);

#endif //PARTICLE_TRAJECTORIES_CALCULATIONS_RECALC_H
