#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "utility.h"
#include "positions.h"
#include "recalc.h"
#include "btot.h"
#include "test.h"

int main(int argc, char *argv[]) {

    const double start_time = omp_get_wtime();

#ifdef TRAJ_TEST
    print_test_notif();
#endif

    unsigned int step_limit = 0, igrf_coefs_version = 0;

#ifdef TRAJ_TEST
    float epsilon = .01f;
    if(argc != 6) {
        throwError("Error: Invalid arguments\nUsage: executable <infile> <outfile> <igrf_ver> <number_of_steps> <epsilon>\n");
    } else {
        char *p;
        igrf_coefs_version = strtol(argv[3], &p, 10);
        if(errno != 0 || *p != '\0' || igrf_coefs_version < 9 || igrf_coefs_version > 13) {
            throwError("Error: Invalid version of IGRF coefficients.\n");
        }

        step_limit = strtol(argv[4], &p, 10);
        if (errno != 0 || *p != '\0') {
            throwError("Error: Value of the step limit was not valid.\n");
        }

        epsilon = strtof(argv[5], &p);
        if (errno != 0 || *p != '\0') {
            throwError("Error: Value of the epsilon was not valid.\n");
        }
    }
#else
    if(argc != 6 || (strncmp(argv[4], "seq", 3) != 0 && strncmp(argv[4], "par", 3) != 0)) {
        throwError("Error: Invalid arguments\nUsage: executable <infile> <outfile> <igrf_ver> <seq/par> <number_of_steps>\n");
    } else {
        char *p;
        igrf_coefs_version = strtol(argv[3], &p, 10);
        if(errno != 0 || *p != '\0' || igrf_coefs_version < 9 || igrf_coefs_version > 13) {
            throwError("Error: Invalid version of IGRF coefficients.\n");
        }

        step_limit = strtol(argv[5], &p, 10);
        if (errno != 0 || *p != '\0') {
            throwError("Error: Value of the step limit was not valid.\n");
        }
    }
#endif

    double rig = 0., del = 0.;
    float rni = 0.f, zn = 0.f, rk = 0.f, r0 = 0.f, the0 = 0.f, fi0 = 0.f, the1 = 0.f, fi1 = 0.f, PARMOD[10] = {0};
    int iy = 0, mes = 0, ide = 0, id = 0, ih = 0, min = 0, is = 0, nk1 = 0, iopt = 0, ist = 0;
    FILE *infil = fopen(argv[1], "r");
    if(!infil) {
        throwError("Error: Input file could not be opened\n");
    } else {
        if(fscanf(infil, "%lf %f %f\n"
                         "%f %f %f\n"
                         "%f %f\n"
                         "%d %d %d %d %d %d %d\n"
                         "%d %d %d %lf\n"
                         "%f %f %f %f\n"
                         "%f %f %f %f %f %f\n",
                         &rig, &zn, &rk,
                         &r0, &the0, &fi0,
                         &the1, &fi1,
                         &iy, &mes, &ide, &id, &ih, &min, &is,
                         &nk1, &iopt, &ist, &del,
                         &PARMOD[1], &PARMOD[0], &PARMOD[2], &PARMOD[3],
                         &PARMOD[4], &PARMOD[5], &PARMOD[6], &PARMOD[7], &PARMOD[8], &PARMOD[9]) != 29) {
            throwError("Error: Invalid number of variables or wrong data are present in the infil file.\n");
        }
        if(del <= 0.){
            throwError("Error: Step value in rigidity has to be higher than 0.0\n");
        }

#ifdef TRAJ_TEST
        if((float) rig != rk) {
            throwError("Error: START and END rigidity must be same for the test.\n");
        }
#endif

        fclose(infil);
        rni = rig;
    }

    // constants
    static const float pi = 3.141592654f;
    static const float rad = pi/180.f;
    static const float q = 1.6021e-19f;
    static const float c = 2.99725e+08f;
    static const float hm0 = 1.6725e-27f;
    static const float Re = 6.3712e+06f;
    //static const int nr = 1; // UNUSED

    FILE *outfil = fopen(argv[2], "w");
    if(!outfil) {
        throwError("Error: Output file could not be created\n");
    }

#ifdef TRAJ_TEST
    write_test_outfil_header(outfil, rig);
    const int run_in_parallel = 0;
    float max_rd = 0;
#else
    fprintf(outfil, "\n\n\n                ASYMPTOTIC COORDINATES\n"
                    "   calculated by model of exter.field T05   \n"
                    " Station with geo.latitude:%9.3f  & longitude:%9.3f & radius : %8.5f\n"
                    " Direction of trajectory with latitude: %9.3f & longitude: %9.3f\n"
                    " Date: %4d %2d %2d  time: %2d hod %2d min %2d sec\n"
                    " Starting rigidity : %8.4f GV Epsilon=%7.4f\n Limit of total number of steps : %d \n\n"
                    " rig : v : rad : eth : efi : ath : afi : time : length\n",
                    the0, fi0, r0, the1, fi1, iy, mes, ide, ih, min, is, rig, del, step_limit);

    const int run_in_parallel = strncmp(argv[4], "par", 3) == 0 ? 1 : 0;
#endif

    float G[105], H[105], REC[105], A1[3], A2[3], A3[3];
    recalc(igrf_coefs_version, iy, &id, &ih, &min, &is, G, H, REC, A1, A2, A3);

    const int rig_values_size = (int)round((rk-rig)/del);
    int nza = 0;
    float rmx1 = 0.f, rmx2 = 0.f, rmi = rk;
    #pragma omp parallel for ordered schedule(dynamic, 1) reduction(+:nza) if (run_in_parallel)
    for(int i = 0; i <= rig_values_size; i++) {
        double threadRig = rig + i*del;
        float alength = 0.f, time = 0.f;
        int nk0 = nk1;
        float a = threadRig*1.0e+09*zn*q / (hm0*c*c);
        float v = c * sqrtf(1.f - 1.f/(1.f + a*a));
        float r = r0;
        float td = 90.f - the0;
        float td1 = 90.f - the1;
        float fei = fi0;
        float fei1 = fi1;
        float th = td * rad;
        float th1 = td1 * rad;
        float f = fei * rad;
        float f1 = fei1 * rad;
        float rr = r * Re;

        float x, y, z;
        sphcar(&rr, &th, &f, &x, &y, &z, 1);
        float xs, ys, zs;
        geogsm(&x, &y, &z, &xs, &ys, &zs, 1, A1, A2, A3);

        float xx, yy, zz;
        set_double_coordinates(&xx, &yy, &zz, &xs, &ys, &zs, &Re);
        float bxs, bys, bzs;
        btot(&xx, &yy, &zz, G, H, REC, A1, A2, A3, PARMOD, &bxs, &bys, &bzs);
        float bx, by, bz;
        geogsm_b(&bx, &by, &bz, &bxs, &bys, &bzs, -1, A1, A2, A3);

        float b = 1.E-09f*sqrtf(bx*bx + by*by + bz*bz);
        float vx = v*sinf(th1)*cosf(f1);
        float vy = v*sinf(th1)*sinf(f1);
        float vz = v*cosf(th1);
        float vv = v;
        unsigned int jj = 0;
        do{
            jj++;

            float vc = vv/c;
            if(vc > 1.f){
                set_rmx_value(&nza, &rmx1, &threadRig);
                break;
            }

            const float hm = hm0 / sqrtf(1.f - vc*vc);
            const float t = (2.f*pi*hm) / (fabsf(zn)*q*b);
            const float h = t/nk0;

            v = vv;
            float d = (1.E-09f*zn*q) / hm;
            const float xp = x;
            const float yp = y;
            const float zp = z;
            const float vxp = vx;
            const float vyp = vy;
            const float vzp = vz;
            calculate_position_and_speed(&h, &d, &bx, &by, &bz, &vx, &vy, &vz, &x, &y, &z, &vv);

            const float gc = (vx*vxp + vy*vyp + vz*vzp)/(vv*v);
            if(gc < 1.f){
                const float gs = sqrtf((1.f-gc) * (1.f+gc));
                const float gu = atan2f(gs, gc);
#ifdef TRAJ_TEST
                if(gu >= epsilon && nk0 <= 50000){
#else
                if(gu >= DEF_TU && nk0 <= 50000){
#endif
                    reassign_previous_coordinates(&nk0, &x, &y, &z, &xp, &yp, &zp, &vx, &vy, &vz, &vxp, &vyp, &vzp);
                    jj--;
                    continue;
                }
            }

            sphcar(&rr, &th, &f, &x, &y, &z, -1);

            float xsp, ysp, zsp;
            geogsm(&x, &y, &z, &xsp, &ysp, &zsp, 1, A1, A2, A3);

            const float pax = (xsp/Re + 25.3f)/36.08f;
            const float paz = (ysp*ysp + zsp*zsp)/(Re*Re);

            r = rr/Re;
            const float paus = pax*pax + paz/459.6736f;
            const float teh = th/rad;
            td = 90.f - teh;
            fei = f/rad;

#ifdef TRAJ_TEST
            float new_rd = fabsf(r - r0) * Re;
            if(new_rd > max_rd) {
                max_rd = new_rd;
            }

#ifdef WRITE_TRAJ
            if(jj % DEF_WRITE_STEP == 0) {
                write_test_outfil_line(outfil, x, y, z, rr, th, f, time, max_rd);
            }
#endif
#ifdef CHECK_R
            if(new_rd > DEF_MAX_RD) {
                break;
            }
#endif
#else // TRAJ_TEST
            if(r > 25.f || paus > 1.f){
                if((r > 25.f && fabsf(r - 25.f) < 0.002f) || (r <= 25.f && paus > 1.f && fabsf(paus - 1.f) < 0.002f)){
                    #pragma omp atomic write
                    rmi = fminf(rmi, threadRig);
                    #pragma omp ordered
                    calculate_trajectory_point(outfil, &threadRig, &th, &f, &vx, &vy, &vz, &rad,
                            &vv, &c, &r, &td, &fei, &alength, &time);
                    break;
                }

                reassign_previous_coordinates(&nk0, &x, &y, &z, &xp, &yp, &zp, &vx, &vy, &vz, &vxp, &vyp, &vzp);
                jj--;
                continue;
            }

            if(r < 1.f){
                set_rmx_value(&nza, &rmx1, &threadRig);
                break;
            }

            if(f > 31.4f){
                set_rmx_value(&nza, &rmx2, &threadRig);
                break;
            }
#endif // TRAJ_TEST

            geogsm(&x, &y, &z, &xs, &ys, &zs, 1, A1, A2, A3);

            set_double_coordinates(&xx, &yy, &zz, &xs, &ys, &zs, &Re);

            btot(&xx, &yy, &zz, G, H, REC, A1, A2, A3, PARMOD, &bxs, &bys, &bzs);
            geogsm_b(&bx, &by, &bz, &bxs, &bys, &bzs, -1, A1, A2, A3);

            b = 1.E-09f*sqrtf(bx*bx + by*by + bz*bz);

            if(jj > step_limit){
                set_rmx_value(&nza, &rmx2, &threadRig);
                break;
            }

            time = time + t/nk0;
            const float a1length = sqrtf((x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp));
            alength = alength + a1length;

            rr = sqrtf(x*x + y*y + z*z);
            nk0 = nk1;
        }while(jj <= step_limit);

#if defined (TRAJ_TEST) && defined (CHECK_R)
        if(jj >= step_limit) {
            printf("TEST PASSED\n");
        } else {
            printf("TEST FAILED\nRadius exceeded allowed deviation of %d.\nFailed on step: %d/%d\n", DEF_MAX_RD, jj, step_limit);
        }
#endif
    }

#ifndef TRAJ_TEST
    calculate_cutoff(outfil, &rmx1, &rmx2, &nza, &rmi, &rni, &del);
#endif
    
    fclose(outfil);

#ifndef TRAJ_TEST
    printf("\nINFO: Successfully finished...\n");
#endif

    printf("The simulation took %lf seconds to finish...\n",
            omp_get_wtime() - start_time);

    return 0;
}
