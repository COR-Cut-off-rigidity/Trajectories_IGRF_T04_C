#include <stdio.h>
#include <stdarg.h>

#include "test.h"

#ifdef TRAJ_TEST

void print_test_notif() {
    printf("\n==========================================\n\n");
    printf("This executable is intended for TESTING!\n");
    printf("\n==========================================\n\n");
}

void write_test_outfil_header(FILE* outfil, double rig) {
    fprintf(outfil, "Particle's rigidity: %g\nIntensity of magnetic field: X = %lf, Y = %lf, Z = %lf\n",
        rig, (double) DEF_BX, (double) DEF_BY, (double) DEF_BZ);

#if defined (WRITE_TRAJ) && defined (DEF_WRITE_STEP)
    fprintf(outfil, "Print step: %d\n", DEF_WRITE_STEP);
#endif

#ifdef WRITE_TRAJ
    fprintf(outfil, "\n");
#ifdef WRITE_VAL_X
    fprintf(outfil, "x ");
#endif
#ifdef WRITE_VAL_Y
    fprintf(outfil, "y ");
#endif
#ifdef WRITE_VAL_Z
    fprintf(outfil, "z ");
#endif
#ifdef WRITE_VAL_R
    fprintf(outfil, "r ");
#endif
#ifdef WRITE_VAL_TH
    fprintf(outfil, "th ");
#endif
#ifdef WRITE_VAL_F
    fprintf(outfil, "f ");
#endif
#ifdef WRITE_VAL_TIME
    fprintf(outfil, "time ");
#endif
#ifdef WRITE_VAL_MAX_RD
    fprintf(outfil, "max_rd ");
#endif
    fprintf(outfil, "\n");
#endif // WRITE_TRAJ
}

void write_test_outfil_line(FILE* outfil, float x, float y, float z, float r, float th, float f, float time, float max_rd) {
#ifdef WRITE_VAL_X
    fprintf(outfil, "%f ", x);
#endif
#ifdef WRITE_VAL_Y
    fprintf(outfil, "%f ", y);
#endif
#ifdef WRITE_VAL_Z
    fprintf(outfil, "%f ", z);
#endif
#ifdef WRITE_VAL_R
    fprintf(outfil, "%f ", r);
#endif
#ifdef WRITE_VAL_TH
    fprintf(outfil, "%f ", th);
#endif
#ifdef WRITE_VAL_F
    fprintf(outfil, "%f ", f);
#endif
#ifdef WRITE_VAL_TIME
    fprintf(outfil, "%f ", time);
#endif
#ifdef WRITE_VAL_MAX_RD
    fprintf(outfil, "%f ", max_rd);
#endif
    fprintf(outfil, "\n");
}

#endif // TRAJ_TEST