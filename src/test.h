#ifdef TRAJ_TEST

#ifndef DEF_BX
#define DEF_BX 0
#endif

#ifndef DEF_BY
#define DEF_BY 86.9263319175807497238
#endif

#ifndef DEF_BZ
#define DEF_BZ 0
#endif

#if defined (CHECK_R) && !defined (DEF_MAX_RD)
#define DEF_MAX_RD 655000
#endif

#if defined (WRITE_TRAJ) && !defined (DEF_WRITE_STEP)
#define DEF_WRITE_STEP 100
#endif

#ifndef TEST_H
#define TEST_H

void print_test_notif();
void write_test_outfil_header(FILE*, double);
void write_test_outfil_line(FILE*, float, float, float, float, float, float, float, float);

#endif // TEST_H

#endif // TRAJ_TEST