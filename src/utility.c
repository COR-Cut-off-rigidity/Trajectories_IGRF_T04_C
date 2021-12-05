//
// Created by pc on 7.4.2020.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utility.h"
#include "test.h"

void throwError(char *error){
    printf("%s", error);
    exit(-1);
}

double custom_pow(double x, int m) {
    unsigned int n = m < 0 ? -m : m;
    double y = n % 2 ? x : 1;
    while (n >>= 1u){
        x = x * x;
        if (n % 2){
            y = y * x;
        }
    }
    return m < 0 ? 1/y : y;
}

float custom_powf(float x, int m) {
    unsigned int n = m < 0 ? -m : m;
    float y = n % 2 ? x : 1;
    while (n >>= 1u){
        x = x * x;
        if (n % 2){
            y = y * x;
        }
    }
    return m < 0 ? 1/y : y;
}

void set_rmx_value(int *nza, float *rmx_value, const double *rig){
    *nza += 1;
    #pragma omp atomic write
    *rmx_value = fmax(*rmx_value, *rig);
}

void reassign_previous_coordinates(int *nk0, float *x, float *y, float *z,
                                   const float *xp, const float *yp, const float *zp,
                                   float *vx, float *vy, float *vz,
                                   const float *vxp, const float *vyp, const float *vzp){
    *nk0 = *nk0*2;
    *x = *xp;
    *y = *yp;
    *z = *zp;
    *vx = *vxp;
    *vy = *vyp;
    *vz = *vzp;
}

void set_double_coordinates(float *xx, float *yy, float *zz, const float *xs, const float *ys, const float *zs,
                            const float *Re){
    *xx = *xs / *Re;
    *yy = *ys / *Re;
    *zz = *zs / *Re;
}

void sphcar(float *R, float *THETA, float *PHI, float *X, float *Y, float *Z, int J){

    if(J > 0) {
        const float SQ = *R * sinf(*THETA);
        *X = SQ * cosf(*PHI);
        *Y = SQ * sinf(*PHI);
        *Z = *R * cosf(*THETA);
    } else {
        float SQ = custom_powf(*X, 2) + custom_powf(*Y, 2);
        *R = sqrtf(SQ + custom_powf(*Z, 2));
        if(SQ != 0.f){
            SQ = sqrtf(SQ);
            *PHI = atan2f(*Y, *X);
            *THETA = atan2f(SQ, *Z);
            if(*PHI < 0){
                *PHI = *PHI + 6.28318531f;
            }
        } else {
            *PHI = 0.f;
            if(*Z < 0.f){
                *THETA = 3.141592654f;
            } else {
                *THETA = 0.f;
            }
        }
    }

}

void geogsm(float *XGEO, float *YGEO, float *ZGEO, float *XGSM, float *YGSM, float *ZGSM, int J,
            const float *A1, const float *A2, const float *A3){
 
    if(J > 0){
        *XGSM = A1[0]**XGEO + A1[1]**YGEO + A1[2]**ZGEO;
        *YGSM = A2[0]**XGEO + A2[1]**YGEO + A2[2]**ZGEO;
        *ZGSM = A3[0]**XGEO + A3[1]**YGEO + A3[2]**ZGEO;
    } else {
        *XGEO = A1[0]**XGSM + A2[0]**YGSM + A3[0]**ZGSM;
        *YGEO = A1[1]**XGSM + A2[1]**YGSM + A3[1]**ZGSM;
        *ZGEO = A1[2]**XGSM + A2[2]**YGSM + A3[2]**ZGSM;
    }

}

void geogsm_b(float *XGEO, float *YGEO, float *ZGEO, float *XGSM, float *YGSM, float *ZGSM, int J,
            const float *A1, const float *A2, const float *A3){

    #ifndef TRAJ_TEST
    
    if(J > 0){
        *XGSM = A1[0]**XGEO + A1[1]**YGEO + A1[2]**ZGEO;
        *YGSM = A2[0]**XGEO + A2[1]**YGEO + A2[2]**ZGEO;
        *ZGSM = A3[0]**XGEO + A3[1]**YGEO + A3[2]**ZGEO;
    } else {
        *XGEO = A1[0]**XGSM + A2[0]**YGSM + A3[0]**ZGSM;
        *YGEO = A1[1]**XGSM + A2[1]**YGSM + A3[1]**ZGSM;
        *ZGEO = A1[2]**XGSM + A2[2]**YGSM + A3[2]**ZGSM;
    }

    #else

    if(J > 0){
        *XGSM = DEF_BX;
        *YGSM = DEF_BY;
        *ZGSM = DEF_BZ;
    } else {
        *XGEO = DEF_BX;
        *YGEO = DEF_BY;
        *ZGEO = DEF_BZ;
    }

    #endif

}

void calculate_trajectory_point(FILE* outfil, const double *rig, const float *th, const float *f,
                                const float *vx, const float *vy, const float *vz, const float *rad,
                                const float *vv, const float *c, const float *r, const float *td, const float *fei,
                                const float *alength, const float *time){

    float vr, vt, vf, at, af, ast, asf;

    carbsp(th, f, vx, vy, vz, &vr, &vt, &vf);
    asymp(th, f, &vr, &vt, &vf, &at, &af);
    ast = at / *rad;
    asf = af / *rad;
    if(asf > 360.f){
        asf -= 360.f;
    }

    #pragma omp ordered
    fprintf(outfil, "  %f\t%.10lf\t%.6f\t%.3f\t%.3f\t%.3f\t%.3f\t%f\t%.2f\n",
            *rig, *vv / *c, *r, *td, *fei, ast, asf, *time, *alength/1000.f);

}

void asymp(const float *th, const float *f, const float *vr, const float *vt, const float *vf, float *at, float *af){

    const float ct = cosf(*th);
    const float st = sinf(*th);
    const float va = *vt*ct + *vr*st;
    const float dp = -*vt*st + *vr*ct;
    const float dm = sqrtf(*vt**vt + va*va);
    *at = atan2f(dp, dm);
    *af = *f + atan2f(*vf, va);

}

void carbsp(const float *th, const float *f, const float *bx, const float *by, const float *bz,
            float *br, float *bt, float *bf){

    const float s = sinf(*th);
    const float c = cosf(*th);
    const float sf = sinf(*f);
    const float cf = cosf(*f);
    const float be = *bx*cf + *by*sf;
    *br = be*s + *bz*c;
    *bf = -*bx*sf + *by*cf;
    *bt = be*c - *bz*s;

}

void calculate_cutoff(FILE* outfil, const float *rmx1, const float *rmx2, const int *nza,
                      const float *rmi, const float *rni, const double *del){

    float rmx = fmaxf(*rmx1, *rmx2);
    const float zan = *nza - (*rmi - *rni) / *del;
    float rms = *rmi + zan**del;
    rmx = rmx + *del;

    fprintf(outfil, "  CUTOFF with rigidities P(S),P(C),P(M) are:\n\t%.5f\t\t%.5f\t\t%.5f\n", *rmi, rmx, rms);

}
