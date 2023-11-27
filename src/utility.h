//
// Created by pc on 7.4.2020.
//

#ifndef PARTICLE_TRAJECTORIES_CALCULATIONS_UTILITY_H
#define PARTICLE_TRAJECTORIES_CALCULATIONS_UTILITY_H

#include <stdio.h>

#ifndef DEF_TU
#define DEF_TU 0.01f
#endif

void throwError(char *error);

double custom_pow(double x, unsigned int m);

float custom_powf(float x, unsigned int m);

void set_rmx_value(int *nza, float *rmx_value, const double *rig);

void reassign_previous_coordinates(int *nk0, float *x, float *y, float *z,
                                   const float *xp, const float *yp, const float *zp,
                                   float *vx, float *vy, float *vz,
                                   const float *vxp, const float *vyp, const float *vzp);

void set_double_coordinates(float *xx, float *yy, float *zz, const float *xs, const float *ys, const float *zs,
                            const float *Re);

void sphcar(float *R, float *THETA, float *PHI, float *X, float *Y, float *Z, int J);

void geogsm(float *XGEO, float *YGEO, float *ZGEO, float *XGSM, float *YGSM, float *ZGSM, int J,
            const float *A1, const float *A2, const float *A3);
            
void geogsm_b(float *XGEO, float *YGEO, float *ZGEO, float *XGSM, float *YGSM, float *ZGSM, int J,
            const float *A1, const float *A2, const float *A3);

void calculate_trajectory_point(FILE* outfil, const double *rig, const float *th,
                                const float *f, const float *vx, const float *vy, const float *vz, const float *rad,
                                const float *vv, const float *c, const float *r, const float *td, const float *fei,
                                const float *alength, const float *time);

void carbsp(const float *th, const float *f, const float *bx, const float *by, const float *bz,
            float *br, float *bt, float *bf);

void asymp(const float *th, const float *f, const float *vr, const float *vt, const float *vf, float *at, float *af);

void calculate_cutoff(FILE* outfil, const float *rmx1, const float *rmx2, const int *nza,
                      const float *rmi, const float *rni, const double *del);

#endif //PARTICLE_TRAJECTORIES_CALCULATIONS_UTILITY_H
