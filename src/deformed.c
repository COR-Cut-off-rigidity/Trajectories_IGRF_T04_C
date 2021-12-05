//
// Created by pc on 7.4.2020.
//
#include <math.h>

#include "deformed.h"
#include "utility.h"

void deformed(const int *IOPT, float PS, const double *X, const double *Y, const double *Z,
              const double *DXSHIFT1, const double *DXSHIFT2, const double *D, const float *DELTADY,
              double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2){

    const double SPS = sinf(PS);
    //const double CPS = sqrt(1. - custom_pow(SPS, 2)); // UNUSED
    const double R2 = custom_pow(*X, 2) + custom_pow(*Y, 2) + custom_pow(*Z, 2);
    const double R = sqrt(R2);
    const double ZR = *Z/R;

    const double RH0 = 7.5;
    const double RH2 = -5.2;
    const double RH = RH0 + RH2* custom_pow(ZR, 2);
    const double DRHDR = -ZR/R*2.*RH2*ZR;
    const double DRHDZ = 2.*RH2*ZR/R;

    const int IEPS = 3;
    const double RRH = R/RH;
    const double F = 1./pow((1. + custom_pow(RRH, IEPS)),(1./IEPS));
    const double DFDR = custom_pow(-RRH, (IEPS - 1))* custom_pow(F, (IEPS + 1))/RH;
    const double DFDRH = -RRH*DFDR;

    const double SPSAS = SPS*F;
    const double CPSAS = sqrt(1. - custom_pow(SPSAS, 2));

    const double XAS = *X*CPSAS - *Z*SPSAS;
    const double ZAS = *X*SPSAS + *Z*CPSAS;

    const double FACPS = SPS/CPSAS*(DFDR + DFDRH*DRHDR)/R;
    const double PSASX = FACPS**X;
    const double PSASY = FACPS**Y;
    const double PSASZ = FACPS**Z + SPS/CPSAS*DFDRH*DRHDZ;

    const double DXASDX = CPSAS - ZAS*PSASX;
    const double DXASDY = -ZAS*PSASY;
    const double DXASDZ = -SPSAS - ZAS*PSASZ;
    const double DZASDX = SPSAS + XAS*PSASX;
    const double DZASDY = XAS*PSASY;
    const double DZASDZ = CPSAS + XAS*PSASZ;
    const double FAC1 = DXASDZ*DZASDY - DXASDY*DZASDZ;
    const double FAC2 = DXASDX*DZASDZ - DXASDZ*DZASDX;
    const double FAC3 = DZASDX*DXASDY - DXASDX*DZASDY;

    double BXAS1, BYAS1, BZAS1, BXAS2, BYAS2, BZAS2;
    warped(IOPT, &PS, &XAS, Y, &ZAS, DXSHIFT1, DXSHIFT2, D, DELTADY, &BXAS1, &BYAS1, &BZAS1, &BXAS2, &BYAS2, &BZAS2);

    *BX1 = BXAS1*DZASDZ - BZAS1*DXASDZ + BYAS1*FAC1;
    *BY1 = BYAS1*FAC2;
    *BZ1 = BZAS1*DXASDX - BXAS1*DZASDX + BYAS1*FAC3;

    *BX2 = BXAS2*DZASDZ - BZAS2*DXASDZ + BYAS2*FAC1;
    *BY2 = BYAS2*FAC2;
    *BZ2 = BZAS2*DXASDX - BXAS2*DZASDX + BYAS2*FAC3;

}

void warped(const int *IOPT, const float *PS, const double *X, const double *Y, const double *Z,
            const double *DXSHIFT1, const double *DXSHIFT2, const double *D, const float *DELTADY,
            double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2){

    const double G = 0.;
    const double DGDX = 0.;
    //const double XL = 20.;
    const double XL4 = 160000;
    const double XL3 = 8000;
    const double DXLDX = 0.;

    const double SPS = sinf(*PS);
    const double RHO2 = custom_pow(*Y, 2) + custom_pow(*Z, 2);
    const double RHO = sqrt(RHO2);

    double PHI, CPHI, SPHI;
    if(*Y == 0. && *Z == 0.){
        PHI = 0.;
        CPHI = 1.;
        SPHI = 0.;
    } else {
        PHI = atan2(*Z,*Y);
        CPHI = *Y/RHO;
        SPHI = *Z/RHO;
    }

    const double RHO2SQ = custom_pow(RHO2, 2);

    const double RR4L4 = RHO/(RHO2SQ + XL4);
    const double F = PHI+G*RHO2*RR4L4*CPHI*SPS;
    const double DFDPHI = 1. - G*RHO2*RR4L4*SPHI*SPS;
    const double DFDRHO = G* custom_pow(RR4L4, 2)*(3.* XL4- RHO2SQ)*CPHI*SPS;
    const double DFDX = RR4L4*CPHI*SPS*(DGDX*RHO2 - G*RHO*RR4L4*4.* XL3*DXLDX);

    const double CF = cos(F);
    const double SF = sin(F);
    const double YAS = RHO*CF;
    const double ZAS = RHO*SF;

    double BX_AS1, BY_AS1, BZ_AS1, BX_AS2, BY_AS2, BZ_AS2;
    unwarped(IOPT, X, &YAS, &ZAS, DXSHIFT1, DXSHIFT2, D, DELTADY, &BX_AS1, &BY_AS1, &BZ_AS1, &BX_AS2, &BY_AS2, &BZ_AS2);

    double BRHO_AS = BY_AS1*CF + BZ_AS1*SF;
    double BPHI_AS = -BY_AS1*SF + BZ_AS1*CF;
    double BRHO_S = BRHO_AS*DFDPHI;
    double BPHI_S = BPHI_AS - RHO*(BX_AS1*DFDX + BRHO_AS*DFDRHO);

    *BX1 = BX_AS1*DFDPHI;
    *BY1 = BRHO_S*CPHI - BPHI_S*SPHI;
    *BZ1 = BRHO_S*SPHI + BPHI_S*CPHI;

    BRHO_AS = BY_AS2*CF + BZ_AS2*SF;
    BPHI_AS = -BY_AS2*SF + BZ_AS2*CF;
    BRHO_S = BRHO_AS*DFDPHI;
    BPHI_S = BPHI_AS - RHO*(BX_AS2*DFDX + BRHO_AS*DFDRHO);

    *BX2 = BX_AS2*DFDPHI;
    *BY2 = BRHO_S*CPHI - BPHI_S*SPHI;
    *BZ2 = BRHO_S*SPHI + BPHI_S*CPHI;

}

void unwarped(const int *IOPT, const double *X, const double *Y, const double *Z,
              const double *DXSHIFT1, const double *DXSHIFT2, const double *D0, const float *DELTADY,
              double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2){

    const double DELTADX1 = 1.;
    const double ALPHA1 = 1.1;
    const double XSHIFT1 = 6.;
    const double DELTADX2 = 0.;
    const double ALPHA2 = 0.25;
    const double XSHIFT2 = 4.;

    static const float A1[60] = {-25.45869857,57.35899080,317.5501869,-2.626756717,
                          -93.38053698,-199.6467926,-858.8129729,34.09192395,845.4214929,
                          -29.07463068,47.10678547,-128.9797943,-781.7512093,6.165038619,
                          167.8905046,492.0680410,1654.724031,-46.77337920,-1635.922669,
                          40.86186772,-.1349775602,-.9661991179E-01,-.1662302354,
                          .002810467517,.2487355077,.1025565237,-14.41750229,-.8185333989,
                          11.07693629,.7569503173,-9.655264745,112.2446542,777.5948964,
                          -5.745008536,-83.03921993,-490.2278695,-1155.004209,39.08023320,
                          1172.780574,-39.44349797,-14.07211198,-40.41201127,-313.2277343,
                          2.203920979,8.232835341,197.7065115,391.2733948,-18.57424451,
                          -437.2779053,23.04976898,11.75673963,13.60497313,4.691927060,
                          18.20923547,27.59044809,6.677425469,1.398283308,2.839005878,
                          31.24817706,24.53577264};

    static const float A2[60] = {-287187.1962,4970.499233,410490.1952,-1347.839052,
                          -386370.3240,3317.983750,-143462.3895,5706.513767,171176.2904,
                          250.8882750,-506570.8891,5733.592632,397975.5842,9771.762168,
                          -941834.2436,7990.975260,54313.10318,447.5388060,528046.3449,
                          12751.04453,-21920.98301,-21.05075617,31971.07875,3012.641612,
                          -301822.9103,-3601.107387,1797.577552,-6.315855803,142578.8406,
                          13161.93640,804184.8410,-14168.99698,-851926.6360,-1890.885671,
                          972475.6869,-8571.862853,26432.49197,-2554.752298,-482308.3431,
                          -4391.473324,105155.9160,-1134.622050,-74353.53091,-5382.670711,
                          695055.0788,-916.3365144,-12111.06667,67.20923358,-367200.9285,
                          -21414.14421,14.75567902,20.75638190,59.78601609,16.86431444,
                          32.58482365,23.69472951,17.24977936,13.64902647,68.40989058,
                          11.67828167};

    const double XM1 = -12.;
    const double XM2 = -12.;

    if(*IOPT != 2) {
        const double XSC1 = (*X - XSHIFT1 - *DXSHIFT1) * ALPHA1 - XM1 * (ALPHA1 - 1.);
        const double YSC1 = *Y * ALPHA1;
        const double ZSC1 = *Z * ALPHA1;
        const double D0SC1 = *D0 * ALPHA1;

        double FX1, FY1, FZ1;
        taildisk(&D0SC1, &DELTADX1, DELTADY, &XSC1, &YSC1, &ZSC1, &FX1, &FY1, &FZ1);

        double HX1, HY1, HZ1;
        shlcar5x5(A1, X, Y, Z, DXSHIFT1, &HX1, &HY1, &HZ1);

        *BX1 = FX1 + HX1;
        *BY1 = FY1 + HY1;
        *BZ1 = FZ1 + HZ1;

        if(*IOPT == 1){
            *BX2 = 0.;
            *BY2 = 0.;
            *BZ2 = 0.;
        }
    }

    const double XSC2 = (*X - XSHIFT2 - *DXSHIFT2)*ALPHA2 - XM2*(ALPHA2 - 1.);
    const double YSC2 = *Y*ALPHA2;
    const double ZSC2 = *Z*ALPHA2;
    const double D0SC2 = *D0*ALPHA2;

    double FX2, FY2, FZ2;
    taildisk(&D0SC2, &DELTADX2, DELTADY, &XSC2, &YSC2, &ZSC2, &FX2, &FY2, &FZ2);

    double HX2, HY2, HZ2;
    shlcar5x5(A2, X, Y, Z, DXSHIFT2, &HX2, &HY2, &HZ2);

    *BX2 = FX2 + HX2;
    *BY2 = FY2 + HY2;
    *BZ2 = FZ2 + HZ2;

    if(*IOPT == 2){
        *BX1 = 0.;
        *BY1 = 0.;
        *BZ1 = 0.;
    }

}

void taildisk(const double *D0, const double *DELTADX, const float *DELTADY,
              const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ){

    static const double F[5] = {-71.09346626,-1014.308601,-1272.939359,-3224.935936,-44546.86232};
    static const double B[5] = {10.90101242,12.68393898,13.51791954,14.86775017,15.12306404};
    static const double C[5] = {.7954069972,.6716601849,1.174866319,2.565249920,10.01986790};

    const double RHO = sqrt(custom_pow(*X, 2) + custom_pow(*Y, 2));
    const double DRHODX = *X/RHO;
    const double DRHODY = *Y/RHO;

    const double DEX = exp(*X/7.);
    const double D = *D0 + *DELTADY* custom_pow((*Y / 20.), 2) + *DELTADX*DEX;
    const double DDDY = *DELTADY**Y*0.005;
    const double DDDX = *DELTADX/7.*DEX;

    const double DZETA = sqrt(custom_pow(*Z, 2) + custom_pow(D, 2));

    const double DDZETADX = D*DDDX/DZETA;
    const double DDZETADY = D*DDDY/DZETA;
    const double DDZETADZ = *Z/DZETA;

    double DBX = 0.;
    double DBY = 0.;
    double DBZ = 0.;
    double BI, CI, S1, S2, DS1DRHO, DS2DRHO, DS1DDZ, DS2DDZ, DS1DX, DS1DY, DS1DZ, DS2DX, DS2DY, DS2DZ,
            S1TS2, S1PS2, S1PS2SQ, FAC1, AS, DASDS1, DASDS2, DASDX, DASDY, DASDZ;
    for(int i=1; i <= 5; i++){
        BI = B[i-1];
        CI = C[i-1];
        const double DZETACI2 = custom_pow((DZETA + CI), 2);
        S1 = sqrt(custom_pow((RHO + BI), 2) + DZETACI2);
        S2 = sqrt(custom_pow((RHO - BI), 2) + DZETACI2);

        DS1DRHO = (RHO+BI) / S1;
        DS2DRHO = (RHO-BI) / S2;
        DS1DDZ = (DZETA+CI) / S1;
        DS2DDZ = (DZETA+CI) / S2;

        DS1DX = DS1DRHO*DRHODX + DS1DDZ*DDZETADX;
        DS1DY = DS1DRHO*DRHODY + DS1DDZ*DDZETADY;
        DS1DZ = DS1DDZ*DDZETADZ;
        DS2DX = DS2DRHO*DRHODX + DS2DDZ*DDZETADX;
        DS2DY = DS2DRHO*DRHODY + DS2DDZ*DDZETADY;
        DS2DZ = DS2DDZ*DDZETADZ;

        S1TS2 = S1*S2;
        S1PS2 = S1 + S2;
        S1PS2SQ = custom_pow(S1PS2, 2);

        FAC1 = sqrt(S1PS2SQ - custom_pow((2. * BI), 2));
        AS = FAC1/(S1TS2*S1PS2SQ);
        DASDS1 = (1./(FAC1*S2) - AS/S1PS2*(S2*S2 + S1*(3.*S1 + 4.*S2))) / (S1*S1PS2);
        DASDS2 = (1./(FAC1*S1) - AS/S1PS2*(S1*S1 + S2*(3.*S2 + 4.*S1))) / (S2*S1PS2);

        DASDX = DASDS1*DS1DX + DASDS2*DS2DX;
        DASDY = DASDS1*DS1DY + DASDS2*DS2DY;
        DASDZ = DASDS1*DS1DZ + DASDS2*DS2DZ;

        DBX = DBX - F[i-1]**X*DASDZ;
        DBY = DBY - F[i-1]**Y*DASDZ;
        DBZ = DBZ + F[i-1]*(2.*AS + *X*DASDX + *Y*DASDY);
    }

    *BX = DBX;
    *BY = DBY;
    *BZ = DBZ;

}

void shlcar5x5(const float *A, const double *X, const double *Y, const double *Z, const double *DSHIFT,
               double *HX, double *HY, double *HZ){

    double RP, CYPI, SYPI, RR, SZRK, CZRK, SQPR, EPR, DBX, DBY, DBZ, COEF;
    int L = 0;
    double DHX = 0.;
    double DHY = 0.;
    double DHZ = 0.;
    for(int i = 1; i <= 5; i++){
        RP = 1./A[49+i];
        CYPI = cos(*Y*RP);
        SYPI = sin(*Y*RP);

        double SQPR_PR = custom_pow(RP, 2);

        for(int k = 1; k <= 5; k++) {
            RR = 1. / A[54+k];
            SQPR = sqrt(SQPR_PR + custom_pow(RR, 2));
            EPR = exp(*X * SQPR);
            SZRK = sin(*Z * RR);
            CZRK = cos(*Z * RR);

            DBX = -SQPR * EPR * CYPI * SZRK;
            DBY = RP * EPR * SYPI * SZRK;
            DBZ = -RR * EPR * CYPI * CZRK;

            L += 2;
            COEF = A[L-2] + A[L-1] * *DSHIFT;

            DHX = DHX + COEF * DBX;
            DHY = DHY + COEF * DBY;
            DHZ = DHZ + COEF * DBZ;
        }
    }

    *HX = DHX;
    *HY = DHY;
    *HZ = DHZ;

}
