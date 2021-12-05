#ifndef PARTICLE_TRAJECTORIES_CALCULATIONS_POSITIONS_H
#define PARTICLE_TRAJECTORIES_CALCULATIONS_POSITIONS_H

void calculate_position_and_speed(const float *h, const float *d, const float *bx, const float *by, const float *bz,
                                  float *vx, float *vy, float *vz, float *x, float *y, float *z, float *vv);

float FNX(const float *d, const float *y, const float *z, const float *by, const float *bz);

float FNY(const float *d, const float *x, const float *z, const float *bx, const float *bz);

float FNZ(const float *d, const float *x, const float *y, const float *bx, const float *by);

#endif //PARTICLE_TRAJECTORIES_CALCULATIONS_POSITIONS_H