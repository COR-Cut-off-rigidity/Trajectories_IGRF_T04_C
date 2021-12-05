//
// Created by pc on 7.4.2020.
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "utility.h"
#include "btot.h"
#include "deformed.h"
#include "birk_tot.h"
#include "full_rc.h"
#include "test.h"

void btot(const float *X, const float *Y, const float *Z, const float *G, const float *H, const float *REC,
          const float *A1, const float *A2, const float *A3, const float *PARMOD, float *BX, float *BY, float *BZ){

    #ifndef TRAJ_TEST

    float PS = 0.f, BXEXT, BYEXT, BZEXT;
    T04_s(PARMOD, &PS, X, Y, Z, &BXEXT, &BYEXT, &BZEXT);
    float HX, HY, HZ;
    igrf_gsm(*X, *Y, *Z, G, H, REC, A1, A2, A3, &HX, &HY, &HZ);
    *BX = HX + BXEXT;
    *BY = HY + BYEXT;
    *BZ = HZ + BZEXT;

    #else

    *BX = DEF_BX;
    *BY = DEF_BY;
    *BZ = DEF_BZ;
    
    #endif

    // IGRF print
    /* float F = sqrtf(HX*HX + HY*HY + HZ*HZ);
    printf("%f\n", F);
    exit(0); */
}

void T04_s(const float *PARMOD, const float *PS, const float *X, const float *Y, const float *Z,
           float *BX, float *BY, float *BZ){

    /*IOPGEN - GENERAL OPTION FLAG:   IOPGEN=0 - CALCULATE TOTAL FIELD
                                      IOPGEN=1 - DIPOLE SHIELDING ONLY
                                      IOPGEN=2 - TAIL FIELD ONLY
                                      IOPGEN=3 - BIRKELAND FIELD ONLY
                                      IOPGEN=4 - RING CURRENT FIELD ONLY
                                      IOPGEN=5 - INTERCONNECTION FIELD ONLY

       IOPT -  TAIL FIELD FLAG:       IOPT=0  -  BOTH MODES
                                      IOPT=1  -  MODE 1 ONLY
                                      IOPT=2  -  MODE 2 ONLY

       IOPB -  BIRKELAND FIELD FLAG:  IOPB=0  -  ALL 4 TERMS
                                      IOPB=1  -  REGION 1, MODES 1 AND 2
                                      IOPB=2  -  REGION 2, MODES 1 AND 2

       IOPR -  RING CURRENT FLAG:     IOPR=0  -  BOTH SRC AND PRC
                                      IOPR=1  -  SRC ONLY
                                      IOPR=2  -  PRC ONLY*/

    const int IOPGEN = 0;
    const int IOPTT = 0;
    const int IOPB = 0;
    const int IOPR = 0;

    double BBX, BBY, BBZ;
    const double XX = *X;
    const double YY = *Y;
    const double ZZ = *Z;
    externfunc(PARMOD, &IOPGEN, &IOPTT, &IOPB, &IOPR, PS, &XX, &YY, &ZZ, &BBX, &BBY, &BBZ);

    *BX = BBX;
    *BY = BBY;
    *BZ = BBZ;

}

void externfunc(const float *PARMOD, const int *IOPGEN, const int *IOPT, const int *IOPB, const int *IOPR,
                const float *PS, const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ) {

    static const double A[69] = {1.00000f,5.44118f,0.891995f,9.09684f,0.00000f,-7.18972f,12.2700f,
                          -4.89408f,0.00000f,0.870536f,1.36081f,0.00000f,0.688650f,0.602330f,
                          0.00000f,0.316346f,1.22728f,-0.363620E-01f,-0.405821f,0.452536f,
                          0.755831f,0.215662f,0.152759f,5.96235f,23.2036f,11.2994f,69.9596f,
                          0.989596f,-0.132131E-01f,0.985681f,0.344212E-01f,1.02389f,0.207867f,
                          1.51220f,0.682715E-01f,1.84714f,1.76977f,1.37690f,0.696350f,0.343280f,
                          3.28846f,111.293f,5.82287f,4.39664f,0.383403f,0.648176f,0.318752E-01f,
                          0.581168f,1.15070f,0.843004f,0.394732f,0.846509f,0.916555f,0.550920f,
                          0.180725f,0.898772f,0.387365f,2.26596f,1.29123f,0.436819f,1.28211f,
                          1.33199f,.405553f,1.6229f,.699074f,1.26131f,2.42297f,.537116f,.619441f};

    //const double BXIMF = 0.; // CURRENTLY UNUSED
    const double PDYN = PARMOD[0];
    const double DST = PARMOD[1]*0.8f - 13.f*sqrt(PDYN);
    const double BYIMF = PARMOD[2];
    const double BZIMF = PARMOD[3];
    const double W1 = PARMOD[4];
    const double W2 = PARMOD[5];
    const double W3 = PARMOD[6];
    const double W4 = PARMOD[7];
    const double W5 = PARMOD[8];
    const double W6 = PARMOD[9];

    const double XAPPA = pow((PDYN / 2.f),A[22]);
    const double XAPPA3 = custom_pow(XAPPA, 3);
    const double XX = *X * XAPPA;
    const double YY = *Y * XAPPA;
    const double ZZ = *Z * XAPPA;

    const double A0_A = 34.586;
    const double A0_S0 = 1.1960;
    const double A0_X0 = 3.4397;
    const double DSIG = 0.005;
    const double SPS = sinf(*PS);
    const double S0 = A0_S0;
    const double AM = A0_A / XAPPA;
    const double X0 = A0_X0 / XAPPA;

    const double FACTIMF = A[19];
    const double OIMFX = 0.;
    const double OIMFY = BYIMF * FACTIMF;
    const double OIMFZ = BZIMF * FACTIMF;

    const double Y2 = custom_pow(*Y, 2);
    const double R = sqrt(custom_pow(*X, 2) + Y2 + custom_pow(*Z, 2));

    const double RH0 = 7.5;
    const double RH2 = -5.2;
    double XSS = *X;
    double ZSS = *Z;
    double RH, SINPSAS, COSPSAS, XSOLD, ZSOLD, DD;
    do {
        RH = RH0 + RH2 * custom_pow((ZSS / R), 2);
        SINPSAS = SPS / pow(1.f + custom_pow((R / RH), 3), 0.33333333);
        COSPSAS = sqrt(1.f - custom_pow(SINPSAS, 2));
        XSOLD = XSS;
        ZSOLD = ZSS;
        ZSS = *X * SINPSAS + *Z * COSPSAS;
        XSS = *X * COSPSAS - *Z * SINPSAS;
        DD = fabs(XSS - XSOLD) + fabs(ZSS - ZSOLD);
    } while (DD > 1E-6);

    double XMXM = AM + XSS - X0;
    if (XMXM < 0) {
        XMXM = 0.;
    }

    const double AXX0 = custom_pow(XMXM, 2);
    const double ASQ = custom_pow(AM, 2);
    const double RHO2 = custom_pow(*Y, 2) + custom_pow(ZSS, 2);
    const double ARO = ASQ + RHO2;
    const double SIGMA = sqrt((ARO + AXX0 + sqrt(custom_pow((ARO + AXX0), 2) - 4.f * ASQ * AXX0)) / (2.f * ASQ));

    double QX = 0., QY = 0., QZ = 0.;
    if(SIGMA < S0+DSIG){
        double BXCF, BYCF, BZCF;
        if(*IOPGEN <= 1){
            double CFX, CFY, CFZ;
            shlcar3x3(&XX, &YY, &ZZ, PS, &CFX, &CFY, &CFZ);
            BXCF = CFX*XAPPA3;
            BYCF = CFY*XAPPA3;
            BZCF = CFZ*XAPPA3;
        } else {
            BXCF = 0.;
            BYCF = 0.;
            BZCF = 0.;
        }

        double ZNAM, BXT1, BYT1, BZT1, BXT2, BYT2, BZT2;
        if(*IOPGEN == 0 || *IOPGEN == 2){
            double DSTT = -20.;
            if(DST < DSTT){
                DSTT = DST;
            }
            ZNAM = pow(fabs(DSTT),0.37f);
            const double DXSHIFT1 = A[23] - A[24]/ZNAM;
            const double DXSHIFT2 = A[25] - A[26]/ZNAM;
            const double D = A[35]*exp(-W1/A[36]) + A[68];
            const float DELTADY = 4.7f;

            deformed(IOPT, *PS, &XX, &YY, &ZZ, &DXSHIFT1, &DXSHIFT2, &D, &DELTADY, &BXT1, &BYT1, &BZT1, &BXT2, &BYT2, &BZT2);
        } else {
            BXT1 = 0.;
            BYT1 = 0.;
            BZT1 = 0.;
            BXT2 = 0.;
            BYT2 = 0.;
            BZT2 = 0.;
        }

        double BXR11, BYR11, BZR11, BXR12, BYR12, BZR12, BXR21, BYR21, BZR21, BXR22, BYR22, BZR22;
        if(*IOPGEN == 0 || *IOPGEN == 3){
            ZNAM = fabs(DST);
            if(DST >= -20.){
                ZNAM = 20.;
            }
            const double XKAPPA1 = A[31]*pow(ZNAM/20.,A[32]);
            const double XKAPPA2 = A[33]*pow(ZNAM/20.,A[34]);

            birk_tot(IOPB, PS, &XX, &YY, &ZZ, &XKAPPA1, &XKAPPA2,
                     &BXR11, &BYR11, &BZR11, &BXR12, &BYR12, &BZR12, &BXR21, &BYR21, &BZR21, &BXR22, &BYR22, &BZR22);
        } else {
            BXR11 = 0.;
            BYR11 = 0.;
            BZR11 = 0.;
            BXR21 = 0.;
            BYR21 = 0.;
            BZR21 = 0.;
        }

        double BXSRC, BYSRC, BZSRC, BXPRC, BYPRC, BZPRC;
        if(*IOPGEN == 0 || *IOPGEN == 4){
            const double PHI = A[37];
            ZNAM = fabs(DST);
            if(DST >= -20.){
                ZNAM = 20.;
            }
            const double SC_SY = A[27]*pow((20./ZNAM),A[28])*XAPPA;
            const double SC_PR = A[29]*pow((20./ZNAM),A[30])*XAPPA;
            full_rc(IOPR, PS, &XX, &YY, &ZZ, &PHI, &SC_SY, &SC_PR, &BXSRC, &BYSRC, &BZSRC, &BXPRC, &BYPRC, &BZPRC);
        } else {
            BXSRC = 0.;
            BYSRC = 0.;
            BZSRC = 0.;
            BXPRC = 0.;
            BYPRC = 0.;
            BZPRC = 0.;
        }

        double HXIMF, HYIMF, HZIMF;
        if(*IOPGEN == 0 || *IOPGEN == 5){
            HXIMF = 0.;
            HYIMF = BYIMF;
            HZIMF = BZIMF;
        } else {
            HXIMF = 0.;
            HYIMF = 0.;
            HZIMF = 0.;
        }

        const double DLP1 = pow((PDYN/2.), A[20]);
        const double DLP2 = pow((PDYN/2.), A[21]);

        const double TAMP1 = A[1] + A[2]*DLP1 + A[3]*(A[38]*W1)/sqrt(custom_pow(W1, 2) + custom_pow(A[38], 2)) + A[4]*DST;
        const double TAMP2 = A[5] + A[6]*DLP2 + A[7]*(A[39]*W2)/sqrt(custom_pow(W2, 2) + custom_pow(A[39], 2)) + A[8]*DST;

        const double A_SRC = A[9] + A[10]*(A[40]*W3)/sqrt(custom_pow(W3, 2) + custom_pow(A[40], 2)) + A[11]*DST;
        const double A_PRC = A[12] + A[13]*(A[41]*W4)/sqrt(custom_pow(W4, 2) + custom_pow(A[41], 2)) + A[14]*DST;

        const double A_R11 = A[15] + A[16]*(A[42]*W5)/sqrt(custom_pow(W5, 2) + custom_pow(A[42], 2));
        const double A_R21 = A[17] + A[18]*(A[43]*W6)/sqrt(custom_pow(W6, 2) + custom_pow(A[43], 2));

        const double BBX = A[0]*BXCF + TAMP1*BXT1 + TAMP2*BXT2 + A_SRC*BXSRC + A_PRC*BXPRC + A_R11*BXR11 + A_R21*BXR21 + A[19]*HXIMF;
        const double BBY = A[0]*BYCF + TAMP1*BYT1 + TAMP2*BYT2 + A_SRC*BYSRC + A_PRC*BYPRC + A_R11*BYR11 + A_R21*BYR21 + A[19]*HYIMF;
        const double BBZ = A[0]*BZCF + TAMP1*BZT1 + TAMP2*BZT2 + A_SRC*BZSRC + A_PRC*BZPRC + A_R11*BZR11 + A_R21*BZR21 + A[19]*HZIMF;

        if(SIGMA < S0-DSIG){
            *BX = BBX;
            *BY = BBY;
            *BZ = BBZ;
        } else {
            dipole(PS, X, Y, Z, &QX, &QY, &QZ);
            const double FINT = 0.5*(1. - (SIGMA-S0)/DSIG);
            const double FEXT = 0.5*(1. + (SIGMA-S0)/DSIG);
            *BX = (BBX + QX)*FINT+OIMFX*FEXT - QX;
            *BY = (BBY + QY)*FINT+OIMFY*FEXT - QY;
            *BZ = (BBZ + QZ)*FINT+OIMFZ*FEXT - QZ;
        }
    } else {
        dipole(PS, X, Y, Z, &QX, &QY, &QZ);
        *BX = OIMFX - QX;
        *BY = OIMFY - QY;
        *BZ = OIMFZ - QZ;
    }

}

void shlcar3x3(const double *X, const double *Y, const double *Z, const float *PS, double *BX, double *BY, double *BZ){

    static const float A[50] = {-901.2327248,895.8011176,817.6208321,-845.5880889,-83.73539535,86.58542841,336.8781402,
                         -329.3619944,-311.2947120,308.6011161,31.94469304,-31.30824526,125.8739681,-372.3384278,
                         -235.4720434,286.7594095,21.86305585,-27.42344605,-150.4874688,2.669338538,1.395023949,
                         -.5540427503,-56.85224007,3.681827033,-43.48705106,5.103131905,1.073551279,-.6673083508,
                         12.21404266,4.177465543,5.799964188,-.3977802319,-1.044652977,.5703560010,3.536082962,
                         -3.222069852,9.620648151,6.082014949,27.75216226,12.44199571,5.122226936,6.982039615,20.12149582,
                         6.150973118,4.663639687,15.73319647,2.303504968,5.840511214,.8385953499E-01,.3477844929};

    const double T1 = A[48];
    const double T2 = A[49];
    const double ST1 = sin(*PS*T1);
    const double CT1 = cos(*PS*T1);
    const double ST2 = sin(*PS*T2);
    const double CT2 = cos(*PS*T2);

    const double X1 = *X*CT1 - *Z*ST1;
    const double Z1 = *X*ST1 + *Z*CT1;
    const double X2 = *X*CT2 - *Z*ST2;
    const double Z2 = *X*ST2 + *Z*CT2;
    const double P1 = A[36];
    const double P12 = 1./custom_pow(P1,2);
    const double R1 = A[39];
    const double R12 = 1./custom_pow(R1,2);

    double CZR = cos(Z1/R1);
    double SZR = sin(Z1/R1);
    double CYP = cos(*Y/P1);
    double SYP = sin(*Y/P1);
    double SQPR = sqrt(P12 + R12);
    double EXPR = exp(SQPR*X1);
    double FX1 = -SQPR*EXPR*CYP*SZR;
    double HY1 = EXPR/P1*SYP*SZR;
    double FZ1 = -EXPR*CYP/R1*CZR;
    double HX1 = FX1*CT1 + FZ1*ST1;
    double HZ1 = -FX1*ST1 + FZ1*CT1;

    const double R2 = A[40];
    const double R22 = 1.f/custom_pow(R2,2);
    SQPR = sqrt(P12 + R22);
    CYP = cos(*Y/P1);
    SYP = sin(*Y/P1);
    CZR = cos(Z1/R2);
    SZR = sin(Z1/R2);
    EXPR = exp(SQPR*X1);
    double FX2 = -SQPR*EXPR*CYP*SZR;
    double HY2 = EXPR/P1*SYP*SZR;
    double FZ2 = -EXPR*CYP/R2*CZR;
    double HX2 = FX2*CT1 + FZ2*ST1;
    double HZ2 = -FX2*ST1 + FZ2*CT1;

    const double R3 = A[41];
    const double R32 = custom_pow(R3,2);
    const double R32_1 = 1./R32;
    SQPR = sqrt(P12 + R32_1);
    CYP = cos(*Y/P1);
    SYP = sin(*Y/P1);
    CZR = cos(Z1/R3);
    SZR = sin(Z1/R3);
    EXPR = exp(SQPR*X1);
    double FX3 = -EXPR*CYP*(SQPR*Z1*CZR + SZR/R3*(X1 + 1./SQPR));
    double HY3 = EXPR/P1*SYP*(Z1*CZR + X1/R3*SZR/SQPR);
    double FZ3 = -EXPR*CYP*(CZR*(1. + X1/R32/SQPR) - Z1/R3*SZR);
    double HX3 = FX3*CT1 + FZ3*ST1;
    double HZ3 = -FX3*ST1 + FZ3*CT1;

    const double P2 = A[37];
    const double P22 = 1./custom_pow(P2,2);
    SQPR = sqrt(P22 + R12);
    CYP = cos(*Y/P2);
    SYP = sin(*Y/P2);
    CZR = cos(Z1/R1);
    SZR = sin(Z1/R1);
    EXPR = exp(SQPR*X1);
    double FX4 = -SQPR*EXPR*CYP*SZR;
    double HY4 = EXPR/P2*SYP*SZR;
    double FZ4 = -EXPR*CYP/R1*CZR;
    double HX4 = FX4*CT1 + FZ4*ST1;
    double HZ4 = -FX4*ST1 + FZ4*CT1;

    SQPR = sqrt(P22 + R22);
    CYP = cos(*Y/P2);
    SYP = sin(*Y/P2);
    CZR = cos(Z1/R2);
    SZR = sin(Z1/R2);
    EXPR = exp(SQPR*X1);
    double FX5 = -SQPR*EXPR*CYP*SZR;
    double HY5 = EXPR/P2*SYP*SZR;
    double FZ5 = -EXPR*CYP/R2*CZR;
    double HX5 = FX5*CT1 + FZ5*ST1;
    double HZ5 = -FX5*ST1 + FZ5*CT1;

    SQPR = sqrt(P22 + R32_1);
    CYP = cos(*Y/P2);
    SYP = sin(*Y/P2);
    CZR = cos(Z1/R3);
    SZR = sin(Z1/R3);
    EXPR = exp(SQPR*X1);
    double FX6 = -EXPR*CYP*(SQPR*Z1*CZR + SZR/R3*(X1 + 1./SQPR));
    double HY6 = EXPR/P2*SYP*(Z1*CZR + X1/R3*SZR/SQPR);
    double FZ6 = -EXPR*CYP*(CZR*(1. + X1/R32/SQPR)-Z1/R3*SZR);
    double HX6 = FX6*CT1 + FZ6*ST1;
    double HZ6 = -FX6*ST1 + FZ6*CT1;

    const double P3 = A[38];
    const double P32 = 1./custom_pow(P3,2);
    SQPR = sqrt(P32 + R12);
    CYP = cos(*Y/P3);
    SYP = sin(*Y/P3);
    CZR = cos(Z1/R1);
    SZR = sin(Z1/R1);
    EXPR = exp(SQPR*X1);
    double FX7 = -SQPR*EXPR*CYP*SZR;
    double HY7 = EXPR/P3*SYP*SZR;
    double FZ7 = -EXPR*CYP/R1*CZR;
    double HX7 = FX7*CT1 + FZ7*ST1;
    double HZ7 = -FX7*ST1 + FZ7*CT1;

    SQPR = sqrt(P32 + R22);
    CYP = cos(*Y/P3);
    SYP = sin(*Y/P3);
    CZR = cos(Z1/R2);
    SZR = sin(Z1/R2);
    EXPR = exp(SQPR*X1);
    double FX8 = -SQPR*EXPR*CYP*SZR;
    double HY8 = EXPR/P3*SYP*SZR;
    double FZ8 = -EXPR*CYP/R2*CZR;
    double HX8 = FX8*CT1 + FZ8*ST1;
    double HZ8 = -FX8*ST1 + FZ8*CT1;

    SQPR = sqrt(P32 + R32_1);
    CYP = cos(*Y/P3);
    SYP = sin(*Y/P3);
    CZR = cos(Z1/R3);
    SZR = sin(Z1/R3);
    EXPR = exp(SQPR*X1);
    double FX9 = -EXPR*CYP*(SQPR*Z1*CZR + SZR/R3*(X1 + 1./SQPR));
    double HY9 = EXPR/P3*SYP*(Z1*CZR + X1/R3*SZR/SQPR);
    double FZ9 = -EXPR*CYP*(CZR*(1. + X1/R32/SQPR) - Z1/R3*SZR);
    double HX9 = FX9*CT1 + FZ9*ST1;
    double HZ9 = -FX9*ST1 + FZ9*CT1;

    const double CPS = cosf(*PS);
    double A1 = A[0] + A[1]*CPS;
    double A2 = A[2] + A[3]*CPS;
    double A3 = A[4] + A[5]*CPS;
    double A4 = A[6] + A[7]*CPS;
    double A5 = A[8] + A[9]*CPS;
    double A6 = A[10] + A[11]*CPS;
    double A7 = A[12] + A[13]*CPS;
    double A8 = A[14] + A[15]*CPS;
    double A9 = A[16] + A[17]*CPS;
    *BX = A1*HX1 + A2*HX2 + A3*HX3 + A4*HX4 + A5*HX5 + A6*HX6 + A7*HX7 + A8*HX8 + A9*HX9;
    *BY = A1*HY1 + A2*HY2 + A3*HY3 + A4*HY4 + A5*HY5 + A6*HY6 + A7*HY7 + A8*HY8 + A9*HY9;
    *BZ = A1*HZ1 + A2*HZ2 + A3*HZ3 + A4*HZ4 + A5*HZ5 + A6*HZ6 + A7*HZ7 + A8*HZ8 + A9*HZ9;

    const double Q1 = A[42];
    const double Q12 = 1./custom_pow(Q1,2);
    const double S1 = A[45];
    const double S12 = 1./custom_pow(S1,2);
    double SQQS = sqrt(Q12 + S12);
    double CYQ = cos(*Y/Q1);
    double SYQ = sin(*Y/Q1);
    double CZS = cos(Z2/S1);
    double SZS = sin(Z2/S1);
    double EXQS = exp(SQQS*X2);

    const double SPS = sinf(*PS);
    FX1 = -SQQS*EXQS*CYQ*CZS*SPS;
    HY1 = EXQS/Q1*SYQ*CZS*SPS;
    FZ1 = EXQS*CYQ/S1*SZS*SPS;
    HX1 = FX1*CT2 + FZ1*ST2;
    HZ1 = -FX1*ST2 + FZ1*CT2;

    const double S2 = A[46];
    const double S22 = 1./custom_pow(S2,2);
    SQQS = sqrt(Q12 + S22);
    //CYQ = cos(*Y/Q1);
    //SYQ = sin(*Y/Q1);
    CZS = cos(Z2/S2);
    SZS = sin(Z2/S2);
    EXQS = exp(SQQS*X2);
    FX2 = -SQQS*EXQS*CYQ*CZS*SPS;
    HY2 = EXQS/Q1*SYQ*CZS*SPS;
    FZ2 = EXQS*CYQ/S2*SZS*SPS;
    HX2 = FX2*CT2 + FZ2*ST2;
    HZ2 = -FX2*ST2 + FZ2*CT2;

    const double S3 = A[47];
    const double S32 = 1./custom_pow(S3,2);
    SQQS = sqrt(Q12+S32);
    //CYQ = cos(*Y/Q1);
    //SYQ = sin(*Y/Q1);
    CZS = cos(Z2/S3);
    SZS = sin(Z2/S3);
    EXQS = exp(SQQS*X2);
    FX3 = -SQQS*EXQS*CYQ*CZS*SPS;
    HY3 = EXQS/Q1*SYQ*CZS*SPS;
    FZ3 = EXQS*CYQ/S3*SZS*SPS;
    HX3 = FX3*CT2 + FZ3*ST2;
    HZ3 = -FX3*ST2 + FZ3*CT2;

    const double Q2 = A[43];
    const double Q22 = 1./custom_pow(Q2,2);
    SQQS = sqrt(Q22 + S12);
    CYQ = cos(*Y/Q2);
    SYQ = sin(*Y/Q2);
    CZS = cos(Z2/S1);
    SZS = sin(Z2/S1);
    EXQS = exp(SQQS*X2);
    FX4 = -SQQS*EXQS*CYQ*CZS*SPS;
    HY4 = EXQS/Q2*SYQ*CZS*SPS;
    FZ4 = EXQS*CYQ/S1*SZS*SPS;
    HX4 = FX4*CT2 + FZ4*ST2;
    HZ4 = -FX4*ST2 + FZ4*CT2;

    SQQS = sqrt(Q22 + S22);
    //CYQ = cos(*Y/Q2);
    //SYQ = sin(*Y/Q2);
    CZS = cos(Z2/S2);
    SZS = sin(Z2/S2);
    EXQS = exp(SQQS*X2);
    FX5 = -SQQS*EXQS*CYQ*CZS*SPS;
    HY5 = EXQS/Q2*SYQ*CZS*SPS;
    FZ5 = EXQS*CYQ/S2*SZS*SPS;
    HX5 = FX5*CT2 + FZ5*ST2;
    HZ5 = -FX5*ST2 + FZ5*CT2;

    SQQS = sqrt(Q22 + S32);
    //CYQ = cos(*Y/Q2);
    //SYQ = sin(*Y/Q2);
    CZS = cos(Z2/S3);
    SZS = sin(Z2/S3);
    EXQS = exp(SQQS*X2);
    FX6 = -SQQS*EXQS*CYQ*CZS*SPS;
    HY6 = EXQS/Q2*SYQ*CZS*SPS;
    FZ6 = EXQS*CYQ/S3*SZS*SPS;
    HX6 = FX6*CT2 + FZ6*ST2;
    HZ6 = -FX6*ST2 + FZ6*CT2;

    const double Q3 = A[44];
    const double Q32 = 1./custom_pow(Q3,2);
    SQQS = sqrt(Q32 + S12);
    CYQ = cos(*Y/Q3);
    SYQ = sin(*Y/Q3);
    CZS = cos(Z2/S1);
    SZS = sin(Z2/S1);
    EXQS = exp(SQQS*X2);
    FX7 = -SQQS*EXQS*CYQ*CZS*SPS;
    HY7 = EXQS/Q3*SYQ*CZS*SPS;
    FZ7 = EXQS*CYQ/S1*SZS*SPS;
    HX7 = FX7*CT2 + FZ7*ST2;
    HZ7 = -FX7*ST2 + FZ7*CT2;

    SQQS= sqrt(Q32 + S22);
    //CYQ = cos(*Y/Q3);
    //SYQ = sin(*Y/Q3);
    CZS = cos(Z2/S2);
    SZS = sin(Z2/S2);
    EXQS = exp(SQQS*X2);
    FX8 = -SQQS*EXQS*CYQ*CZS*SPS;
    HY8 = EXQS/Q3*SYQ*CZS*SPS;
    FZ8 = EXQS*CYQ/S2*SZS*SPS;
    HX8 = FX8*CT2 + FZ8*ST2;
    HZ8 = -FX8*ST2 + FZ8*CT2;

    SQQS = sqrt(Q32 + S32);
    //CYQ = cos(*Y/Q3);
    //SYQ = sin(*Y/Q3);
    CZS = cos(Z2/S3);
    SZS = sin(Z2/S3);
    EXQS = exp(SQQS*X2);
    FX9 = -SQQS*EXQS*CYQ*CZS*SPS;
    HY9 = EXQS/Q3*SYQ*CZS*SPS;
    FZ9 = EXQS*CYQ/S3*SZS*SPS;
    HX9 = FX9*CT2 + FZ9*ST2;
    HZ9 = -FX9*ST2 + FZ9*CT2;

    const double S2PS = 2.*CPS;
    A1 = A[18] + A[19]*S2PS;
    A2 = A[20] + A[21]*S2PS;
    A3 = A[22] + A[23]*S2PS;
    A4 = A[24] + A[25]*S2PS;
    A5 = A[26] + A[27]*S2PS;
    A6 = A[28] + A[29]*S2PS;
    A7 = A[30] + A[31]*S2PS;
    A8 = A[32] + A[33]*S2PS;
    A9 = A[34] + A[35]*S2PS;

    *BX = *BX + A1*HX1 + A2*HX2 + A3*HX3 + A4*HX4 + A5*HX5 + A6*HX6 + A7*HX7 + A8*HX8 + A9*HX9;
    *BY = *BY + A1*HY1 + A2*HY2 + A3*HY3 + A4*HY4 + A5*HY5 + A6*HY6 + A7*HY7 + A8*HY8 + A9*HY9;
    *BZ = *BZ + A1*HZ1 + A2*HZ2 + A3*HZ3 + A4*HZ4 + A5*HZ5 + A6*HZ6 + A7*HZ7 + A8*HZ8 + A9*HZ9;

}

void dipole(const float *PS, const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ){

    const double SPS = sinf(*PS);
    const double CPS = cosf(*PS);
    const double P = custom_pow(*X, 2);
    const double U = custom_pow(*Z, 2);
    const double V = 3.**Z**X;
    const double T = custom_pow(*Y, 2);
    const double Q = 30115./ custom_pow(sqrt(P + T + U), 5);
    *BX = Q*((T + U - 2.*P)*SPS - V*CPS);
    *BY = -3.**Y*Q*(*X*SPS + *Z*CPS);
    *BZ = Q*((P + T - 2.*U)*CPS - V*SPS);

}

void igrf_gsm(float XGSM, float YGSM, float ZGSM, const float *G, const float *H, const float *REC,
              const float *A1, const float *A2, const float *A3, float *HXGSM, float *HYGSM, float *HZGSM){

    float A[14], B[14];

    float XGEO, YGEO, ZGEO;
    geogsm(&XGEO, &YGEO, &ZGEO, &XGSM, &YGSM, &ZGSM, -1, A1, A2, A3);

    const float RHO2 = custom_powf(XGEO, 2) + custom_powf(YGEO, 2);
    const float R = sqrtf(RHO2 + custom_powf(ZGEO, 2));
    const float C = ZGEO/R;
    const float RHO = sqrtf(RHO2);
    const float S = RHO/R;
    float CF, SF;
    if(S < 1.E-5){
        CF = 1.;
        SF = 0.;
    } else {
        CF = XGEO/RHO;
        SF = YGEO/RHO;
    }

    const float PP = 1./R;
    float P = PP;

    const int IRP3 = (int) R + 2;
    int NM = 3 + 30/IRP3;
    if(NM > 13){
        NM = 13;
    }
    const int K = NM+1;
    for(int i = 1; i <= K; i++) {
        P = P * PP;
        A[i-1] = P;
        B[i-1] = P * i;
    }
    P = 1.;

    int MM = 0, MN;
    float D = 0.f, BBR = 0.f, BBT = 0.f, BBF = 0.f, X = 0., Y = 1., W, Q, Z, BI, P2, D2, AN, E, HH, QQ, XK, DP, PM;
    for(int i = 1; i <= K; i++){
        if(i != 1){
            MM = i - 1;
            W = X;
            X = W*CF + Y*SF;
            Y = Y*CF - W*SF;
        }
        Q = P;
        Z = D;
        BI = 0.;
        P2 = 0.;
        D2 = 0.;

        for(int j = i; j <= K; j++){
            AN = A[j-1];
            MN = j*(j - 1)/2 + i;
            E = G[MN-1];
            HH = H[MN-1];
            W = E*Y + HH*X;
            BBR = BBR + B[j-1]*W*Q;
            BBT = BBT - AN*W*Z;
            if(i != 1){
                QQ = Q;
                if(S < 1.E-5) {
                    QQ = Z;
                }
                BI = BI + AN*(E*X - HH*Y)*QQ;
            }
            XK = REC[MN-1];
            DP = C*Z - S*Q - XK*D2;
            PM = C*Q - XK*P2;
            D2 = Z;
            P2 = Q;
            Z = DP;
            Q = PM;
        }
        D = S*D + C*P;
        P = S*P;
        if(i != 1){
            BI = BI*MM;
            BBF = BBF + BI;
        }
    }

    const float BR = BBR;
    const float BT = BBT;
    float BF;
    if(S < 1.E-5){
        if(C < 0.){
            BBF = -BBF;
        }
        BF = BBF;
    } else {
        BF = BBF/S;
    }
    const float HE = BR*S + BT*C;
    float HXGEO = HE*CF - BF*SF;
    float HYGEO = HE*SF + BF*CF;
    float HZGEO = BR*C - BT*S;

    geogsm(&HXGEO, &HYGEO, &HZGEO, HXGSM, HYGSM, HZGSM, 1, A1, A2, A3);

}