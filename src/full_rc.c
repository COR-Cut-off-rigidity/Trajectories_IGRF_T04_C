//
// Created by pc on 7.4.2020.
//
#include <math.h>

#include "full_rc.h"
#include "utility.h"

void full_rc(const int *IOPR, const float *PS, const double *X, const double *Y, const double *Z, const double *PHI,
             const double *SC_SY, const double *SC_PR, double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC,
             double *BYPRC, double *BZPRC){

    static const double C_SY[86] = {-957.2534900f,-817.5450246f,583.2991249f,758.8568270f,
                             13.17029064f,68.94173502f,-15.29764089f,-53.43151590f,27.34311724f,
                             149.5252826f,-11.00696044f,-179.7031814f,953.0914774f,817.2340042f,
                             -581.0791366f,-757.5387665f,-13.10602697f,-68.58155678f,15.22447386f,
                             53.15535633f,-27.07982637f,-149.1413391f,10.91433279f,179.3251739f,
                             -6.028703251f,1.303196101f,-1.345909343f,-1.138296330f,-0.06642634348f,
                             -0.3795246458f,.07487833559f,.2891156371f,-.5506314391f,-.4443105812f,
                             0.2273682152f,0.01086886655f,-9.130025352f,1.118684840f,1.110838825f,
                             .1219761512f,-.06263009645f,-.1896093743f,.03434321042f,.01523060688f,
                             -.4913171541f,-.2264814165f,-.04791374574f,.1981955976f,-68.32678140f,
                             -48.72036263f,14.03247808f,16.56233733f,2.369921099f,6.200577111f,
                             -1.415841250f,-0.8184867835f,-3.401307527f,-8.490692287f,3.217860767f,
                             -9.037752107f,66.09298105f,48.23198578f,-13.67277141f,-16.27028909f,
                             -2.309299411f,-6.016572391f,1.381468849f,0.7935312553f,3.436934845f,
                             8.260038635f,-3.136213782f,8.833214943f,8.041075485f,8.024818618f,
                             35.54861873f,12.55415215f,1.738167799f,3.721685353f,23.06768025f,
                             6.871230562f,6.806229878f,21.35990364f,1.687412298f,3.500885177f,
                             0.3498952546f,0.6595919814f};

    static const double C_PR[86] = {-64820.58481f,-63965.62048f,66267.93413f,135049.7504f,
                             -36.56316878f,124.6614669f,56.75637955f,-87.56841077f,5848.631425f,
                             4981.097722f,-6233.712207f,-10986.40188f,68716.52057f,65682.69473f,
                             -69673.32198f,-138829.3568f,43.45817708f,-117.9565488f,-62.14836263f,
                             79.83651604f,-6211.451069f,-5151.633113f,6544.481271f,11353.03491f,
                             23.72352603f,-256.4846331f,25.77629189f,145.2377187f,-4.472639098f,
                             -3.554312754f,2.936973114f,2.682302576f,2.728979958f,26.43396781f,
                             -9.312348296f,-29.65427726f,-247.5855336f,-206.9111326f,74.25277664f,
                             106.4069993f,15.45391072f,16.35943569f,-5.965177750f,-6.079451700f,
                             115.6748385f,-35.27377307f,-32.28763497f,-32.53122151f,93.74409310f,
                             84.25677504f,-29.23010465f,-43.79485175f,-6.434679514f,-6.620247951f,
                             2.443524317f,2.266538956f,-43.82903825f,6.904117876f,12.24289401f,
                             17.62014361f,152.3078796f,124.5505289f,-44.58690290f,-63.02382410f,
                             -8.999368955f,-9.693774119f,3.510930306f,3.770949738f,-77.96705716f,
                             22.07730961f,20.46491655f,18.67728847f,9.451290614f,9.313661792f,
                             644.7620970f,418.2515954f,7.183754387f,35.62128817f,19.43180682f,
                             39.57218411f,15.69384715f,7.123215241f,2.300635346f,21.90881131f,
                             -.01775839370f,.3996346710f};

    double HXSRC, HYSRC, HZSRC, HXPRC, HYPRC, HZPRC;
    src_prc(IOPR, SC_SY, SC_PR, PHI, PS, X, Y, Z, &HXSRC, &HYSRC, &HZSRC, &HXPRC, &HYPRC, &HZPRC);

    double X_SC = *SC_SY - 1.;
    double FSX, FSY, FSZ;
    if(*IOPR == 0 || *IOPR == 1){
        rc_shield(C_SY, PS, &X_SC, X, Y, Z, &FSX, &FSY, &FSZ);
    } else {
        FSX = 0.;
        FSY = 0.;
        FSZ = 0.;
    }

    X_SC = *SC_PR - 1.;
    double FPX, FPY, FPZ;
    if(*IOPR == 0 || *IOPR == 2){
        rc_shield(C_PR, PS, &X_SC, X, Y, Z, &FPX, &FPY, &FPZ);
    } else {
        FPX = 0.;
        FPY = 0.;
        FPZ = 0.;
    }

    *BXSRC = HXSRC + FSX;
    *BYSRC = HYSRC + FSY;
    *BZSRC = HZSRC + FSZ;

    *BXPRC = HXPRC + FPX;
    *BYPRC = HYPRC + FPY;
    *BZPRC = HZPRC + FPZ;

}

void src_prc(const int *IOPR, const double *SC_SY, const double *SC_PR, const double *PHI, const float *PS,
             const double *X, const double *Y, const double *Z, double *BXSRC, double *BYSRC, double *BZSRC,
             double *BXPRC, double *BYPRC, double *BZPRC){

    const double CPS = cosf(*PS);
    const double SPS = sinf(*PS);

    const double XT = *X*CPS - *Z*SPS;
    const double ZT = *Z*CPS + *X*SPS;

    const double XTS = XT / *SC_SY;
    const double YTS = *Y / *SC_SY;
    const double ZTS = ZT / *SC_SY;

    const double XTA = XT / *SC_PR;
    const double YTA = *Y / *SC_PR;
    const double ZTA = ZT / *SC_PR;

    double BXS = 0.;
    double BYS = 0.;
    double BZS = 0.;
    if(*IOPR <= 1){
        rc_symm(&XTS, &YTS, &ZTS, &BXS, &BYS, &BZS);
    }

    double BXA_S = 0.;
    double BYA_S = 0.;
    double BZA_S = 0.;
    if(*IOPR == 0 || *IOPR == 2){
        prc_symm(&XTA, &YTA, &ZTA, &BXA_S, &BYA_S, &BZA_S);
    }

    const double CP = cos(*PHI);
    const double SP = sin(*PHI);
    const double XR = XTA*CP - YTA*SP;
    const double YR = XTA*SP + YTA*CP;

    double BXA_QR = 0.;
    double BYA_QR = 0.;
    double BZA_Q = 0.;
    if(*IOPR == 0 || *IOPR == 2){
        prc_quad(&XR, &YR, &ZTA, &BXA_QR, &BYA_QR, &BZA_Q);
    }

    const double BXA_Q = BXA_QR*CP + BYA_QR*SP;
    const double BYA_Q = -BXA_QR*SP + BYA_QR*CP;
    const double BXP = BXA_S + BXA_Q;
    const double BYP = BYA_S + BYA_Q;
    const double BZP = BZA_S + BZA_Q;

    *BXSRC = BXS*CPS + BZS*SPS;
    *BYSRC = BYS;
    *BZSRC = BZS*CPS - BXS*SPS;

    *BXPRC = BXP*CPS + BZP*SPS;
    *BYPRC = BYP;
    *BZPRC = BZP*CPS - BXP*SPS;

}

void rc_symm(const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ){

    const double DS = 1.E-2;
    const double DC = 0.99994999875;
    const double D = 1.E-4;
    const double DRD = 5E3;

    const double RHO2 = custom_pow(*X, 2) + custom_pow(*Y, 2);
    const double R2 = RHO2 + custom_pow(*Z, 2);
    const double R = sqrt(R2);
    const double RP = R + D;
    const double RM = R - D;
    const double SINT = sqrt(RHO2)/R;
    const double COST = *Z/R;

    if(SINT < DS){
        const double A = ap(&R,&DS,&DC)/DS;
        const double DARDR = (RP*ap(&RP,&DS,&DC) - RM*ap(&RM,&DS,&DC))*DRD;
        const double FXY = *Z*(2.*A - DARDR)/(R*R2);
        *BX = FXY**X;
        *BY = FXY**Y;
        *BZ = (2.*A* custom_pow(COST, 2)+DARDR* custom_pow(SINT, 2))/R;
    } else {
        const double THETA = atan2(SINT,COST);
        const double TP = THETA + D;
        const double TM = THETA - D;
        const double SINTP = sin(TP);
        const double SINTM = sin(TM);
        const double COSTP = cos(TP);
        const double COSTM = cos(TM);
        const double BR = (SINTP*ap(&R,&SINTP,&COSTP) - SINTM*ap(&R,&SINTM,&COSTM))/(R*SINT)*DRD;
        const double BT = (RM*ap(&RM,&SINT,&COST) - RP*ap(&RP,&SINT,&COST))/R*DRD;
        const double FXY = (BR + BT*COST/SINT)/R;
        *BX = FXY**X;
        *BY = FXY**Y;
        *BZ = BR*COST - BT*SINT;
    }

}

double ap(const double *R, const double *SINT, const double *COST){

    const double A1 = -456.5289941f;
    const double A2 = 375.9055332f;
    const double RRC1 = 4.274684950f;
    const double DD1 = 2.439528329f;
    const double RRC2 = 3.367557287f;
    const double DD2 = 3.146382545f;
    const double P1 = -0.2291904607f;
    const double R1 = 3.746064740f;
    const double DR1 = 1.508802177f;
    const double DLA1 = 0.5873525737f;
    const double P2 = 0.1556236119f;
    const double R2 = 4.993638842f;
    const double DR2 = 3.324180497f;
    const float DLA2 = 0.4368407663f;
    const double P3 = 0.1855957207f;
    const double R3 = 2.969226745f;
    const double DR3 = 2.243367377f;

    int PROX = 0;
    double SINT1 = *SINT;
    double COST1 = *COST;
    if(SINT1 < 1.E-2){
        SINT1 = 1.E-2;
        COST1 = .99994999875f;
        PROX = 1;
    }

    const double ARG1 = -custom_pow(((*R - R1) / DR1), 2)- custom_pow((COST1 / DLA1), 2);
    const double ARG2 = -custom_pow(((*R - R2) / DR2), 2)- custom_pow((COST1 / DLA2), 2);
    const double ARG3 = -custom_pow(((*R - R3) / DR3), 2);

    double DEXP1;
    if(ARG1 < -500.){
        DEXP1 = 0.;
    } else {
        DEXP1 = exp(ARG1);
    }

    double DEXP2;
    if(ARG2 < -500.){
        DEXP2 = 0.;
    } else {
        DEXP2 = exp(ARG2);
    }

    double DEXP3;
    if(ARG3 < -500.){
        DEXP3 = 0.;
    } else {
        DEXP3 = exp(ARG3);
    }

    const double ALPHA = custom_pow(SINT1, 2) / *R;
    const double ALPHA_S = ALPHA*(1. + P1*DEXP1 + P2*DEXP2 + P3*DEXP3);

    const double GAMMA = COST1/ custom_pow(*R, 2);
    const double GAMMA_S = GAMMA;
    const double GAMMAS2 = custom_pow(GAMMA_S, 2);

    const double ALSQH = custom_pow(ALPHA_S, 2)/2.;
    const double F = 64./27.*GAMMAS2 + custom_pow(ALSQH, 2);
    const double Q = pow((sqrt(F) + ALSQH),(1./3.));
    double C = Q - 4.*pow(GAMMAS2,(1./3.))/(3.*Q);
    if(C < 0.){
        C = 0.;
    }
    const double G = sqrt(custom_pow(C, 2) + 4.*pow(GAMMAS2,(1./3.)));
    const double RS = 4./((sqrt(2.*G - C)+sqrt(C))*(G + C));
    const double COSTS = GAMMA_S* custom_pow(RS, 2);
    const double SINTS = sqrt(1. - custom_pow(COSTS, 2));
    const double RHOS = RS*SINTS;
    //const double RHOS2 = custom_pow(RHOS, 2); // UNUSED
    const double ZS = RS*COSTS;
    const double ZS2 = custom_pow(ZS, 2);

    double P = custom_pow((RRC1 + RHOS), 2) + ZS2 + custom_pow(DD1, 2);
    double XK2 = 4.*RRC1*RHOS/P;
    double XK = sqrt(XK2);
    double XKRHO12 = XK*sqrt(RHOS);

    double XK2S = 1. - XK2;
    double DL = log(1./XK2S);
    double ELK = 1.38629436112 + XK2S*(0.09666344259 + XK2S*(0.03590092383f + XK2S*(0.03742563713f + XK2S*0.01451196212f))) +
                 DL*(0.5 + XK2S*(0.12498593597 + XK2S*(0.06880248576 + XK2S*(0.03328355346 + XK2S*0.00441787012))));
    double ELE = 1. + XK2S*(0.44325141463 + XK2S*(0.0626060122 + XK2S*(0.04757383546 + XK2S*0.01736506451))) +
                 DL*XK2S*(0.2499836831 + XK2S*(0.09200180037 + XK2S*(0.04069697526 + XK2S*0.00526449639)));

    const double APHI1 = ((1. - XK2*0.5)*ELK - ELE)/XKRHO12;

    P = custom_pow((RRC2 + RHOS), 2) + ZS2 + custom_pow(DD2, 2);
    XK2 = 4.*RRC2*RHOS/P;
    XK = sqrt(XK2);
    XKRHO12 = XK*sqrt(RHOS);
    XK2S = 1. - XK2;
    DL = log(1./XK2S);
    ELK = 1.38629436112 + XK2S*(0.09666344259 + XK2S*(0.03590092383f + XK2S*(0.03742563713f + XK2S*0.01451196212f))) +
          DL*(0.5 + XK2S*(0.12498593597 + XK2S*(0.06880248576 + XK2S*(0.03328355346 + XK2S*0.00441787012))));
    ELE = 1. + XK2S*(0.44325141463 + XK2S*(0.0626060122 + XK2S*(0.04757383546 + XK2S*0.01736506451))) +
          DL*XK2S*(0.2499836831 + XK2S*(0.09200180037 + XK2S*(0.04069697526 + XK2S*0.00526449639)));

    const double APHI2 = ((1. - XK2*0.5)*ELK-ELE)/XKRHO12;

    double AP = A1*APHI1 + A2*APHI2;
    if(PROX){
        AP = AP**SINT/SINT1;
    }

    return AP;
}

void prc_symm(const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ){

    const double DS = 1.E-2;
    const double DC = 0.99994999875;
    const double D = 1.E-4;
    const double DRD = 5.E3;

    const double RHO2 = custom_pow(*X, 2) + custom_pow(*Y, 2);
    const double R2 = RHO2 + custom_pow(*Z, 2);
    const double R = sqrt(R2);
    const double RP = R + D;
    const double RM = R - D;
    const double SINT = sqrt(RHO2)/R;
    const double COST = *Z/R;

    double FXY;
    if(SINT < DS){
        const double A = apprc(&R,&DS,&DC)/DS;
        const double DARDR = (RP*apprc(&RP,&DS,&DC) - RM*apprc(&RM,&DS,&DC))*DRD;
        FXY = *Z*(2.*A - DARDR)/(R*R2);
        *BZ = (2.*A* custom_pow(COST, 2) + DARDR* custom_pow(SINT, 2))/R;
    } else {
        const double THETA = atan2(SINT,COST);
        const double TP = THETA + D;
        const double TM = THETA - D;
        const double SINTP = sin(TP);
        const double SINTM = sin(TM);
        const double COSTP = cos(TP);
        const double COSTM = cos(TM);
        const double BR = (SINTP*apprc(&R,&SINTP,&COSTP) - SINTM*apprc(&R,&SINTM,&COSTM))/(R*SINT)*DRD;
        const double BT = (RM*apprc(&RM,&SINT,&COST) - RP*apprc(&RP,&SINT,&COST))/R*DRD;
        FXY = (BR + BT*COST/SINT)/R;
        *BZ = BR*COST - BT*SINT;
    }
    *BX = FXY**X;
    *BY = FXY**Y;

}

double apprc(const double *R, const double *SINT, const double *COST){

    const double A1 = -80.11202281f;
    const double A2 = 12.58246758f;
    const double RRC1 = 6.560486035f;
    const double DD1 = 1.930711037f;
    const double RRC2 = 3.827208119f;
    const double DD2 = .7789990504f;
    const double P1 = .3058309043f;
    const double ALPHA1 = .1817139853f;
    const double DAL1 = .1257532909f;
    const double BETA1 = 3.422509402f;
    const double DG1 = .04742939676f;
    const double P2 = -4.800458958f;
    const double ALPHA2 = -.02845643596f;
    const double DAL2 = .2188114228f;
    const double BETA2 = 2.545944574f;
    const double DG2 = .00813272793f;
    const double BETA3 = .35868244f;
    const double P3 = 103.1601001f;
    const double ALPHA3 = -.00764731187f;
    const double DAL3 = .1046487459f;
    const double BETA4 = 2.958863546f;
    const double DG3 = .01172314188f;
    const double BETA5 = .4382872938f;
    const double Q0 = .01134908150f;
    const double Q1 = 14.51339943f;
    const double ALPHA4 = .2647095287f;
    const double DAL4 = .07091230197f;
    const double DG4 = .01512963586f;
    const double Q2 = 6.861329631f;
    const double ALPHA5 = .1677400816f;
    const double DAL5 = .04433648846f;
    const double DG5 = .05553741389f;
    const double BETA6 = .7665599464f;
    const double BETA7 = .7277854652f;

    int PROX = 0;
    double SINT1 = *SINT;
    double COST1 = *COST;
    if(SINT1 < 1.E-2){
        SINT1 = 1.E-2;
        COST1 = .99994999875;
        PROX = 1;
    }

    const double ALPHA = custom_pow(SINT1, 2)/ *R;
    const double GAMMA = COST1/pow(*R,2);

    const double ARG1 = -custom_pow((GAMMA / DG1), 2);
    const double ARG2 = -custom_pow(((ALPHA - ALPHA4) / DAL4), 2)- custom_pow((GAMMA / DG4), 2);

    double DEXP1;
    if(ARG1 < -500.){
        DEXP1 = 0.;
    } else {
        DEXP1 = exp(ARG1);
    }

    double DEXP2;
    if(ARG2 < -500.){
        DEXP2 = 0.;
    } else {
        DEXP2 = exp(ARG2);
    }

    const double ALPHA_S = ALPHA*(1. + P1/pow((1. + pow(((ALPHA-ALPHA1)/DAL1),2)),BETA1)*DEXP1 +
                                  P2*(ALPHA - ALPHA2)/pow((1. + pow(((ALPHA-ALPHA2)/DAL2),2)),BETA2)/pow((1. + pow((GAMMA/DG2),2)),BETA3) +
                                  P3*pow((ALPHA - ALPHA3),2)/pow((1. + pow(((ALPHA-ALPHA3)/DAL3),2)),BETA4)/pow((1. + pow((GAMMA/DG3),2)),BETA5));

    const double GAMMA_S = GAMMA*(1. + Q0+Q1*(ALPHA - ALPHA4)*DEXP2 +
                                  Q2*(ALPHA - ALPHA5)/pow((1. + pow(((ALPHA - ALPHA5)/DAL5),2)),BETA6)/pow((1. + pow((GAMMA/DG5),2)),BETA7));

    const double GAMMAS2 = custom_pow(GAMMA_S, 2);

    const double ALSQH = custom_pow(ALPHA_S, 2)/2.;
    const double F = 64./27.*GAMMAS2 + custom_pow(ALSQH, 2);
    const double Q = pow((sqrt(F) + ALSQH),(1./3.));
    double C = Q - 4.*pow(GAMMAS2,(1./3.))/(3.*Q);
    if(C < 0.){
        C = 0.;
    }
    const double G = sqrt(custom_pow(C, 2) + 4.*pow(GAMMAS2,(1./3.)));
    const double RS = 4./((sqrt(2.*G - C) + sqrt(C))*(G + C));
    const double COSTS = GAMMA_S* custom_pow(RS, 2);
    const double SINTS = sqrt(1. - custom_pow(COSTS, 2));
    const double RHOS = RS*SINTS;
    //const double RHOS2 = custom_pow(RHOS, 2); // UNUSED
    const double ZS = RS*COSTS;
    const double ZS2 = custom_pow(ZS, 2);

    double P = custom_pow((RRC1 + RHOS), 2) + ZS2 + custom_pow(DD1, 2);
    double XK2 = 4.*RRC1*RHOS/P;
    double XK = sqrt(XK2);
    double XKRHO12 = XK*sqrt(RHOS);

    double XK2S = 1. - XK2;
    double DL = log(1./XK2S);
    double ELK = 1.38629436112 + XK2S*(0.09666344259 + XK2S*(0.03590092383f + XK2S*(0.03742563713f + XK2S*0.01451196212f))) +
                 DL*(0.5 + XK2S*(0.12498593597 + XK2S*(0.06880248576 + XK2S*(0.03328355346 + XK2S*0.00441787012))));
    double ELE = 1. + XK2S*(0.44325141463 + XK2S*(0.0626060122 + XK2S*(0.04757383546 + XK2S*0.01736506451))) +
                 DL*XK2S*(0.2499836831 + XK2S*(0.09200180037 + XK2S*(0.04069697526 + XK2S*0.00526449639)));

    const double APHI1 = ((1. - XK2*0.5)*ELK - ELE)/XKRHO12;

    P = custom_pow((RRC2 + RHOS), 2) + ZS2 + custom_pow(DD2, 2);
    XK2 = 4.*RRC2*RHOS/P;
    XK = sqrt(XK2);
    XKRHO12 = XK*sqrt(RHOS);

    XK2S = 1. - XK2;
    DL = log(1./XK2S);
    ELK = 1.38629436112 + XK2S*(0.09666344259 + XK2S*(0.03590092383f + XK2S*(0.03742563713f + XK2S*0.01451196212f))) +
          DL*(0.5 + XK2S*(0.12498593597 + XK2S*(0.06880248576 + XK2S*(0.03328355346 + XK2S*0.00441787012))));
    ELE = 1. + XK2S*(0.44325141463 + XK2S*(0.0626060122 + XK2S*(0.04757383546 + XK2S*0.01736506451))) +
          DL*XK2S*(0.2499836831 + XK2S*(0.09200180037 + XK2S*(0.04069697526 + XK2S*0.00526449639)));

    const double APHI2 = ((1. - XK2*0.5)*ELK - ELE)/XKRHO12;
    double APPRC = A1*APHI1 + A2*APHI2;
    if(PROX) {
        APPRC = APPRC**SINT/SINT1;
    }
    return APPRC;
}

void prc_quad(const double *X, const double *Y, const double *Z, double *BX, double *BY, double *BZ){

    const double D = 1.E-4;
    const double DD = 2.E-4;
    const double DS = 1.E-2;
    const double DC = 0.99994999875;

    const double RHO2 = custom_pow(*X, 2) + custom_pow(*Y, 2);
    const double R = sqrt(RHO2 + custom_pow(*Z, 2));
    const double RHO = sqrt(RHO2);
    const double SINT = RHO/R;
    const double COST = *Z/R;
    const double RP = R + D;
    const double RM = R - D;

    if(SINT > DS){
        const double CPHI = *X/RHO;
        const double SPHI = *Y/RHO;
        const double BR = br_prc_q(&R,&SINT,&COST);
        const double BT = bt_prc_q(&R,&SINT,&COST);
        const double DBRR = (br_prc_q(&RP,&SINT,&COST) - br_prc_q(&RM,&SINT,&COST))/DD;
        const double THETA = atan2(SINT,COST);
        const double TP = THETA+D;
        const double TM = THETA-D;
        const double SINTP = sin(TP);
        const double COSTP = cos(TP);
        const double SINTM = sin(TM);
        const double COSTM = cos(TM);
        const double DBTT = (bt_prc_q(&R,&SINTP,&COSTP) - bt_prc_q(&R,&SINTM,&COSTM))/DD;
        *BX = SINT*(BR + (BR + R*DBRR + DBTT)* custom_pow(SPHI, 2))+COST*BT;
        *BY = -SINT*SPHI*CPHI*(BR + R*DBRR + DBTT);
        *BZ = (BR*COST - BT*SINT)*CPHI;
    } else {
        const double ST = DS;
        double CT = DC;
        if(*Z < 0.){
            CT = -DC;
        }
        const double THETA = atan2(ST,CT);
        const double TP = THETA+D;
        const double TM = THETA-D;
        const double SINTP = sin(TP);
        const double COSTP = cos(TP);
        const double SINTM = sin(TM);
        const double COSTM = cos(TM);
        const double BR = br_prc_q(&R,&ST,&CT);
        const double BT = bt_prc_q(&R,&ST,&CT);
        const double DBRR = (br_prc_q(&RP,&ST,&CT) - br_prc_q(&RM,&ST,&CT))/DD;
        const double DBTT = (bt_prc_q(&R,&SINTP,&COSTP) - bt_prc_q(&R,&SINTM,&COSTM))/DD;
        const double FCXY = R*DBRR+DBTT;
        const double RST2 = custom_pow((R * ST), 2);
        *BX = (BR*(custom_pow(*X, 2) + 2.* custom_pow(*Y, 2)) + FCXY* custom_pow(*Y, 2))/ RST2 + BT*COST;
        *BY = -(BR + FCXY)**X**Y/ RST2;
        *BZ = (BR*COST/ST - BT)**X/R;
    }

}

double br_prc_q(const double *R, const double *SINT, const double *COST){

    const double A1 = -21.2666329f;
    const double A2 = 32.24527521f;
    const double A3 = -6.062894078f;
    const double A4 = 7.515660734f;
    const double A5 = 233.7341288f;
    const double A6 = -227.1195714f;
    const double A7 = 8.483233889f;
    const double A8 = 16.80642754f;
    const double A9 = -24.63534184f;
    const double A10 = 9.067120578f;
    const double A11 = -1.052686913f;
    const double A12 = -12.08384538f;
    const double A13 = 18.61969572f;
    const double A14 = -12.71686069f;
    const double A15 = 47017.35679f;
    const double A16 = -50646.71204f;
    const double A17 = 7746.058231f;
    const double A18 = 1.531069371f;
    const double XK1 = 2.318824273f;
    const double AL1 = .1417519429f;
    const double DAL1 = .6388013110E-02f;
    const double B1 = 5.303934488f;
    const double BE1 = 4.213397467f;
    const double XK2 = .7955534018f;
    const double AL2 = .1401142771f;
    const double DAL2 = .2306094179E-01f;
    const double B2 = 3.462235072f;
    const double BE2 = 2.568743010f;
    const double XK3 = 3.477425908f;
    const double XK4 = 1.922155110f;
    const double AL3 = .1485233485f;
    const double DAL3 = .2319676273E-01f;
    const double B3 = 7.830223587f;
    const double BE3 = 8.492933868f;
    const double AL4 = .1295221828f;
    const double DAL4 = .01753008801f;
    const double DG1 = .01125504083f;
    const double AL5 = .1811846095f;
    const double DAL5 = .04841237481f;
    const double DG2 = .01981805097f;
    const double C1 = 6.557801891f;
    const double C2 = 6.348576071f;
    const double C3 = 5.744436687f;
    const double AL6 = .2265212965f;
    const double DAL6 = .1301957209f;
    const double DRM = .5654023158f;

    const double SINT2 = custom_pow(*SINT, 2);
    const double COST2 = custom_pow(*COST, 2);
    const double SC = *SINT**COST;
    const double ALPHA = SINT2 / *R;
    const double GAMMA = *COST/ custom_pow(*R, 2);

    double F = 0.;
    double FA = 0.;
    double FS = 0.;
    calc_ffs(&ALPHA, &AL1, &DAL1, &F, &FA, &FS);
    const double D1 = SC*pow(F,XK1)/(pow((*R/B1),BE1) + 1.);
    const double D2 = D1*COST2;

    calc_ffs(&ALPHA, &AL2, &DAL2, &F, &FA, &FS);
    const double D3 = SC*pow(FS,XK2)/(pow((*R/B2),BE2) + 1.);
    const double D4 = D3*COST2;

    calc_ffs(&ALPHA, &AL3, &DAL3, &F, &FA, &FS);
    const double D5 = SC*(pow(ALPHA,XK3))*(pow(FS,XK4))/(pow((*R/B3),BE3) + 1.);
    const double D6 = D5*COST2;

    double ARGA = custom_pow(((ALPHA - AL4) / DAL4), 2) + 1.;
    double ARGG = 1. + custom_pow((GAMMA / DG1), 2);

    const double D7 = SC/ARGA/ARGG;
    const double D8 = D7/ARGA;
    const double D9 = D8/ARGA;
    const double D10 = D9/ARGA;

    ARGA = custom_pow(((ALPHA - AL5) / DAL5), 2) + 1.;
    ARGG = 1. + custom_pow((GAMMA / DG2), 2);

    const double D11 = SC/ARGA/ARGG;
    const double D12 = D11/ARGA;
    const double D13 = D12/ARGA;
    const double D14 = D13/ARGA;

    const double R4 = custom_pow(*R, 4);

    const double D15 = SC/(R4 + custom_pow(C1, 4));
    const double D16 = SC/(R4 + custom_pow(C2, 4))*COST2;
    const double D17 = SC/(R4 + custom_pow(C3, 4))* custom_pow(COST2, 2);

    calc_ffs(&ALPHA, &AL6, &DAL6, &F, &FA, &FS);
    const double D18 = SC*FS/(1. + custom_pow(((*R - 1.2) / DRM), 2));

    const double BR_PRC_Q = A1*D1 + A2*D2 + A3*D3 + A4*D4 + A5*D5 + A6*D6 + A7*D7 + A8*D8 + A9*D9 +
                            A10*D10 + A11*D11 + A12*D12 + A13*D13 + A14*D14 + A15*D15 + A16*D16 + A17*D17 + A18*D18;

    return BR_PRC_Q;

}

double bt_prc_q(const double *R, const double *SINT, const double *COST){

    const double A1 = 12.74640393f;
    const double A2 = -7.516393516f;
    const double A3 = -5.476233865f;
    const double A4 = 3.212704645f;
    const double A5 = -59.10926169f;
    const double A6 = 46.62198189f;
    const double A7 = -.01644280062f;
    const double A8 = .1234229112f;
    const double A9 = -.08579198697f;
    const double A10 = .01321366966f;
    const double A11 = .8970494003f;
    const double A12 = 9.136186247f;
    const double A13 = -38.19301215f;
    const double A14 = 21.73775846f;
    const double A15 = -410.0783424f;
    const double A16 = -69.90832690f;
    const double A17 = -848.8543440f;
    const double XK1 = 1.243288286f;
    const double AL1 = .2071721360f;
    const double DAL1 = .05030555417f;
    const double B1 = 7.471332374f;
    const double BE1 = 3.180533613f;
    const double XK2 = 1.376743507f;
    const double AL2 = .1568504222f;
    const double DAL2 = .02092910682f;
    const double BE2 = 1.985148197f;
    const double XK3 = .3157139940f;
    const double XK4 = 1.056309517f;
    const double AL3 = .1701395257f;
    const double DAL3 = .1019870070f;
    const double B3 = 6.293740981f;
    const double BE3 = 5.671824276f;
    const double AL4 = .1280772299f;
    const double DAL4 = .02189060799f;
    const double DG1 = .01040696080f;
    const double AL5 = .1648265607f;
    const double DAL5 = .04701592613f;
    const double DG2 = .01526400086f;
    const double C1 = 12.88384229f;
    const double C2 = 3.361775101f;
    const double C3 = 23.44173897f;

    const double SINT2 = custom_pow(*SINT, 2);
    const double COST2 = custom_pow(*COST, 2);
    //const double SC = *SINT**COST; // UNUSED
    const double ALPHA = SINT2 / *R;
    const double GAMMA = *COST/ custom_pow(*R, 2);

    double F = 0.;
    double FA = 0.;
    double FS = 0.;
    calc_ffs(&ALPHA, &AL1, &DAL1, &F, &FA, &FS);
    const double D1 = pow(F,XK1)/(pow((*R/B1),BE1) + 1.);
    const double D2 = D1*COST2;

    calc_ffs(&ALPHA, &AL2, &DAL2, &F, &FA, &FS);
    const double D3 = pow(FA,XK2)/pow(*R,BE2);
    const double D4 = D3*COST2;

    calc_ffs(&ALPHA, &AL3, &DAL3, &F, &FA, &FS);
    const double D5 = pow(FS,XK3)*pow(ALPHA,XK4)/(pow((*R/B3),BE3) + 1.);
    const double D6 = D5*COST2;

    const double A0_ZERO = 0.f;
    calc_ffs(&GAMMA, &A0_ZERO, &DG1, &F, &FA, &FS);
    const double FCC = (1. + custom_pow(((ALPHA - AL4) / DAL4), 2));
    const double D7 = 1./FCC*FS;
    const double D8 = D7/FCC;
    const double D9 = D8/FCC;
    const double D10 = D9/FCC;

    const double ARG = 1. + custom_pow(((ALPHA - AL5) / DAL5), 2);
    const double D11 = 1./ARG/(1. + custom_pow((GAMMA / DG2), 2));
    const double D12 = D11/ARG;
    const double D13 = D12/ARG;
    const double D14 = D13/ARG;

    const double R4 = custom_pow(*R, 4);

    const double D15 = 1./(R4 + custom_pow(C1, 2));
    const double D16 = COST2/(R4 + custom_pow(C2, 2));
    const double D17 = custom_pow(COST2, 2)/(R4 + custom_pow(C3, 2));

    const double BT_PRC_Q = A1*D1 + A2*D2 + A3*D3 + A4*D4 + A5*D5 + A6*D6 + A7*D7 + A8*D8 + A9*D9 +
                            A10*D10 + A11*D11 + A12*D12 + A13*D13 + A14*D14 + A15*D15 + A16*D16 + A17*D17;

    return BT_PRC_Q;
}

void calc_ffs(const double *A, const double *A0, const double *DA, double *F, double *FA, double *FS){

    const double DA2 = custom_pow(*DA, 2);
    const double SQ1 = sqrt(custom_pow((*A + *A0), 2) + DA2);
    const double SQ2 = sqrt(custom_pow((*A - *A0), 2) + DA2);
    *FA = 2./(SQ1+SQ2);
    *F = *FA * *A;
    *FS = 0.5*(SQ1 + SQ2)/(SQ1*SQ2)*(1. - *F * *F);

}

void rc_shield(const double *A, const float *PS, const double *X_SC, const double *X, const double *Y, const double *Z,
               double *BX, double *BY, double *BZ){

    const double FAC_SC = custom_pow((*X_SC + 1.), 3);

    const double CPS = cosf(*PS);
    const double SPS = sinf(*PS);
    const double S3PS = 2.*CPS;

    const double PST1 = *PS*A[84];
    const double PST2 = *PS*A[85];
    const double ST1 = sin(PST1);
    const double CT1 = cos(PST1);
    const double ST2 = sin(PST2);
    const double CT2 = cos(PST2);

    const double X1 = *X*CT1 - *Z*ST1;
    const double Z1 = *X*ST1 + *Z*CT1;
    const double X2 = *X*CT2 - *Z*ST2;
    const double Z2 = *X*ST2 + *Z*CT2;

    int L = 0;
    double GX = 0.;
    double GY = 0.;
    double GZ = 0.;
    double P, Q, CYPI, CYQI, SYPI, SYQI, R, S, SZRK, CZSK, CZRK, SZSK, SQPR, SQQS,
            EPR, EQS, FX, FY, FZ, HX, HY, HZ, HXR, HZR;
    for(int i = 1; i <= 2; i++){
        for(int j = 1; j <= 3; j++){
            P = A[71+j];
            Q = A[77+j];
            CYPI = cos(*Y/P);
            CYQI = cos(*Y/Q);
            SYPI = sin(*Y/P);
            SYQI = sin(*Y/Q);

            double P_SQPR = 1./ custom_pow(P, 2);
            double Q_SQQS = 1./ custom_pow(Q, 2);

            for(int k = 1; k <= 3; k++){
                R = A[74+k];
                S = A[80+k];
                SZRK = sin(Z1/R);
                CZSK = cos(Z2/S);
                CZRK = cos(Z1/R);
                SZSK = sin(Z2/S);
                SQPR = sqrt(P_SQPR + 1./ custom_pow(R, 2));
                SQQS = sqrt(Q_SQQS + 1./ custom_pow(S, 2));
                EPR = exp(X1*SQPR);
                EQS = exp(X2*SQQS);
                for(int l = 1; l <= 2; l++){
                    for(int m = 1; m <= 2; m++){
                        if(i == 1){
                            FX = -SQPR*EPR*CYPI*SZRK*FAC_SC;
                            FY = EPR*SYPI*SZRK/P*FAC_SC;
                            FZ = -EPR*CYPI*CZRK/R*FAC_SC;
                            if(l == 1){
                                if(m == 1) {
                                    HX = FX;
                                    HY = FY;
                                    HZ = FZ;
                                } else {
                                    HX = FX**X_SC;
                                    HY = FY**X_SC;
                                    HZ = FZ**X_SC;
                                }
                            } else {
                                if(m == 1){
                                    HX = FX*CPS;
                                    HY = FY*CPS;
                                    HZ = FZ*CPS;
                                } else {
                                    HX = FX*CPS**X_SC;
                                    HY = FY*CPS**X_SC;
                                    HZ = FZ*CPS**X_SC;
                                }
                            }
                        } else {
                            FX = -SPS*SQQS*EQS*CYQI*CZSK*FAC_SC;
                            FY = SPS/Q*EQS*SYQI*CZSK*FAC_SC;
                            FZ = SPS/S*EQS*CYQI*SZSK*FAC_SC;
                            if(l == 1){
                                if(m == 1){
                                    HX = FX;
                                    HY = FY;
                                    HZ = FZ;
                                } else {
                                    HX = FX**X_SC;
                                    HY = FY**X_SC;
                                    HZ = FZ**X_SC;
                                }
                            } else {
                                if(m == 1){
                                    HX = FX*S3PS;
                                    HY = FY*S3PS;
                                    HZ = FZ*S3PS;
                                } else {
                                    HX = FX*S3PS**X_SC;
                                    HY = FY*S3PS**X_SC;
                                    HZ = FZ*S3PS**X_SC;
                                }
                            }
                        }

                        L = L + 1;
                        if(i == 1){
                            HXR = HX*CT1 + HZ*ST1;
                            HZR = -HX*ST1 + HZ*CT1;
                        } else {
                            HXR = HX*CT2 + HZ*ST2;
                            HZR = -HX*ST2 + HZ*CT2;
                        }

                        GX = GX + HXR*A[L-1];
                        GY = GY + HY*A[L-1];
                        GZ = GZ + HZR*A[L-1];
                    }
                }
            }
        }
    }

    *BX = GX;
    *BY = GY;
    *BZ = GZ;

}