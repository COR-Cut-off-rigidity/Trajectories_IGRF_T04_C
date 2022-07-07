//
// Created by pc on 7.4.2020.
//
#include <math.h>

#include "birk_tot.h"
#include "utility.h"

void birk_tot(const int *IOPB, const float *PS,
              const double *X, const double *Y, const double *Z, const double *XKAPPA1, const double *XKAPPA2,
              double *BX11, double *BY11, double *BZ11, double *BX12, double *BY12, double *BZ12,
              double *BX21, double *BY21, double *BZ21, double *BX22, double *BY22, double *BZ22){

    static const float SH11[86] = {46488.84663,-15541.95244,-23210.09824,-32625.03856,
                            -109894.4551,-71415.32808,58168.94612,55564.87578,-22890.60626,
                            -6056.763968,5091.368100,239.7001538,-13899.49253,4648.016991,
                            6971.310672,9699.351891,32633.34599,21028.48811,-17395.96190,
                            -16461.11037,7447.621471,2528.844345,-1934.094784,-588.3108359,
                            -32588.88216,10894.11453,16238.25044,22925.60557,77251.11274,
                            50375.97787,-40763.78048,-39088.60660,15546.53559,3559.617561,
                            -3187.730438,309.1487975,88.22153914,-243.0721938,-63.63543051,
                            191.1109142,69.94451996,-187.9539415,-49.89923833,104.0902848,
                            -120.2459738,253.5572433,89.25456949,-205.6516252,-44.93654156,
                            124.7026309,32.53005523,-98.85321751,-36.51904756,98.88241690,
                            24.88493459,-55.04058524,61.14493565,-128.4224895,-45.35023460,
                            105.0548704,-43.66748755,119.3284161,31.38442798,-92.87946767,
                            -33.52716686,89.98992001,25.87341323,-48.86305045,59.69362881,
                            -126.5353789,-44.39474251,101.5196856,59.41537992,41.18892281,
                            80.86101200,3.066809418,7.893523804,30.56212082,10.36861082,
                            8.222335945,19.97575641,2.050148531,4.992657093,2.300564232,
                            .2256245602,-.05841594319};

    static const float SH12[86] = {210260.4816,-1443587.401,-1468919.281,281939.2993,
                            -1131124.839,729331.7943,2573541.307,304616.7457,468887.5847,
                            181554.7517,-1300722.650,-257012.8601,645888.8041,-2048126.412,
                            -2529093.041,571093.7972,-2115508.353,1122035.951,4489168.802,
                            75234.22743,823905.6909,147926.6121,-2276322.876,-155528.5992,
                            -858076.2979,3474422.388,3986279.931,-834613.9747,3250625.781,
                            -1818680.377,-7040468.986,-414359.6073,-1295117.666,-346320.6487,
                            3565527.409,430091.9496,-.1565573462,7.377619826,.4115646037,
                            -6.146078880,3.808028815,-.5232034932,1.454841807,-12.32274869,
                            -4.466974237,-2.941184626,-.6172620658,12.64613490,1.494922012,
                            -21.35489898,-1.652256960,16.81799898,-1.404079922,-24.09369677,
                            -10.99900839,45.94237820,2.248579894,31.91234041,7.575026816,
                            -45.80833339,-1.507664976,14.60016998,1.348516288,-11.05980247,
                            -5.402866968,31.69094514,12.28261196,-37.55354174,4.155626879,
                            -33.70159657,-8.437907434,36.22672602,145.0262164,70.73187036,
                            85.51110098,21.47490989,24.34554406,31.34405345,4.655207476,
                            5.747889264,7.802304187,1.844169801,4.867254550,2.941393119,
                            .1379899178,.06607020029};

    static const float SH21[86] = {162294.6224,503885.1125,-27057.67122,-531450.1339,
                            84747.05678,-237142.1712,84133.61490,259530.0402,69196.05160,
                            -189093.5264,-19278.55134,195724.5034,-263082.6367,-818899.6923,
                            43061.10073,863506.6932,-139707.9428,389984.8850,-135167.5555,
                            -426286.9206,-109504.0387,295258.3531,30415.07087,-305502.9405,
                            100785.3400,315010.9567,-15999.50673,-332052.2548,54964.34639,
                            -152808.3750,51024.67566,166720.0603,40389.67945,-106257.7272,
                            -11126.14442,109876.2047,2.978695024,558.6019011,2.685592939,
                            -338.0004730,-81.99724090,-444.1102659,89.44617716,212.0849592,
                            -32.58562625,-982.7336105,-35.10860935,567.8931751,-1.917212423,
                            -260.2023543,-1.023821735,157.5533477,23.00200055,232.0603673,
                            -36.79100036,-111.9110936,18.05429984,447.0481000,15.10187415,
                            -258.7297813,-1.032340149,-298.6402478,-1.676201415,180.5856487,
                            64.52313024,209.0160857,-53.85574010,-98.52164290,14.35891214,
                            536.7666279,20.09318806,-309.7349530,58.54144539,67.45226850,
                            97.92374406,4.752449760,10.46824379,32.91856110,12.05124381,
                            9.962933904,15.91258637,1.804233877,6.578149088,2.515223491,
                            .1930034238,-.02261109942};

    static const float SH22[86] = {-131287.8986,-631927.6885,-318797.4173,616785.8782,
                            -50027.36189,863099.9833,47680.20240,-1053367.944,-501120.3811,
                            -174400.9476,222328.6873,333551.7374,-389338.7841,-1995527.467,
                            -982971.3024,1960434.268,297239.7137,2676525.168,-147113.4775,
                            -3358059.979,-2106979.191,-462827.1322,1017607.960,1039018.475,
                            520266.9296,2627427.473,1301981.763,-2577171.706,-238071.9956,
                            -3539781.111,94628.16420,4411304.724,2598205.733,637504.9351,
                            -1234794.298,-1372562.403,-2.646186796,-31.10055575,2.295799273,
                            19.20203279,30.01931202,-302.1028550,-14.78310655,162.1561899,
                            .4943938056,176.8089129,-.2444921680,-100.6148929,9.172262228,
                            137.4303440,-8.451613443,-84.20684224,-167.3354083,1321.830393,
                            76.89928813,-705.7586223,18.28186732,-770.1665162,-9.084224422,
                            436.3368157,-6.374255638,-107.2730177,6.080451222,65.53843753,
                            143.2872994,-1028.009017,-64.22739330,547.8536586,-20.58928632,
                            597.3893669,10.17964133,-337.7800252,159.3532209,76.34445954,
                            84.74398828,12.76722651,27.63870691,32.69873634,5.145153451,
                            6.310949163,6.996159733,1.971629939,4.436299219,2.904964304,
                            .1486276863,.06859991529};

    double X_SC = *XKAPPA1 - 1.1;

    if(*IOPB == 0 || *IOPB == 1){
        int NUMB = 1, MODE = 1;
        double FX11, FY11, FZ11, HX11, HY11, HZ11;
        birk_1n2(&NUMB, &MODE, PS, X, Y, Z, XKAPPA1, &FX11, &FY11, &FZ11);
        birk_shl(SH11, PS, &X_SC, X, Y, Z, &HX11, &HY11, &HZ11);
        *BX11 = FX11 + HX11;
        *BY11 = FY11 + HY11;
        *BZ11 = FZ11 + HZ11;

        MODE = 2;
        double FX12, FY12, FZ12, HX12, HY12, HZ12;
        birk_1n2(&NUMB, &MODE, PS, X, Y, Z, XKAPPA1, &FX12, &FY12, &FZ12);
        birk_shl(SH12, PS, &X_SC, X, Y, Z, &HX12, &HY12, &HZ12);
        *BX12 = FX12 + HX12;
        *BY12 = FY12 + HY12;
        *BZ12 = FZ12 + HZ12;
    }

    X_SC = *XKAPPA2 - 1.0;

    if(*IOPB == 0 || *IOPB == 2){
        int NUMB = 2, MODE = 1;
        double FX21, FY21, FZ21, HX21, HY21, HZ21;
        birk_1n2(&NUMB, &MODE, PS, X, Y, Z, XKAPPA2, &FX21, &FY21, &FZ21);
        birk_shl(SH21, PS, &X_SC, X, Y, Z, &HX21, &HY21, &HZ21);
        *BX21 = FX21 + HX21;
        *BY21 = FY21 + HY21;
        *BZ21 = FZ21 + HZ21;

        MODE = 2;
        double FX22, FY22, FZ22, HX22, HY22, HZ22;
        birk_1n2(&NUMB, &MODE, PS, X, Y, Z, XKAPPA2, &FX22, &FY22, &FZ22);
        birk_shl(SH22, PS, &X_SC, X, Y, Z, &HX22, &HY22, &HZ22);
        *BX22 = FX22 + HX22;
        *BY22 = FY22 + HY22;
        *BZ22 = FZ22 + HZ22;
    }

}

void birk_1n2(const int *NUMB, const int *MODE, const float *PS, const double *X, const double *Y, const double *Z,
              const double *XKAPPA, double *BX, double *BY, double *BZ){

    const double BETA = 0.9;
    const double RH = 10.;
    const double EPS = 3.;

    static const float A11[31] = {.1618068350,-.1797957553,2.999642482,-.9322708978,
                           -.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
                           -16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
                           1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
                           -1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
                           .09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
                           2.492118385,.7113544659};

    static const float A12[31] = {.7058026940,-.2845938535,5.715471266,-2.472820880,
                           -.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
                           -212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
                           2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
                           -1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
                           .1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
                           1.212634762,.5567714182};

    static const float A21[31] = {.1278764024,-.2320034273,1.805623266,-32.37241440,
                           -.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
                           -6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
                           1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
                           -1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
                           .1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
                           1.102649543,.8867880020};

    static const float A22[31] = {.4036015198,-.3302974212,2.827730930,-45.44405830,
                           -1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
                           -233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
                           .7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
                           -1.460805289,.7719653528,-.6658988668,.2515179349E-05,
                           .02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
                           2.503482679,1.071587299,.7247997430};

    double DPHI = 0., DTHETA = 0.;
    if(*NUMB == 1){
        DPHI = 0.055;
        DTHETA = 0.06;
    }

    if(*NUMB == 2){
        DPHI = 0.030;
        DTHETA = 0.09;
    }

    const double Xsc = *X**XKAPPA;
    const double Ysc = *Y**XKAPPA;
    const double Zsc = *Z**XKAPPA;
    const double Xsc2 = custom_pow(Xsc, 2);
    const double Zsc2 = custom_pow(Zsc, 2);
    
    const double RHO = sqrt(Xsc2 + Zsc2);
    const double RHO2 = custom_pow(RHO, 2);
    const double Rsc = sqrt(Xsc2 + custom_pow(Ysc, 2) + Zsc2);
    const double RHO_0 = 7.0;
    const double RHO_02 = custom_pow(RHO_0, 2);

    double PHI;
    if(Xsc == 0. && Zsc == 0.){
        PHI = 0.;
    } else {
        PHI = atan2(-Zsc,Xsc);
    }

    const double SPHIC = sin(PHI);
    const double CPHIC = cos(PHI);

    const double B = 0.5;
    const double BRACK = DPHI + B*RHO_02/(RHO_02 + 1.)*(RHO2 - 1.)/(RHO_02 + RHO2);
    const double R1RH = (Rsc - 1.)/RH;
    const double PSIAS = BETA**PS/pow((1. + pow(R1RH,EPS)),(1./EPS));

    const double PHIS = PHI - BRACK*sin(PHI) - PSIAS;
    const double DPHISPHI = 1. - BRACK*cos(PHI);
    const double DPHIS_part = BETA**PS*pow(R1RH,(EPS - 1.));
    const double DPHIS_part_2 = (RH*Rsc*pow((1. + pow(R1RH,EPS)),(1./EPS + 1.)));
    const double DPHISRHO = -2.*B*RHO_02*RHO/ custom_pow((RHO_02 + RHO2), 2)*sin(PHI) + DPHIS_part*RHO/DPHIS_part_2;
    const double DPHISDY = DPHIS_part*Ysc/DPHIS_part_2;

    const double SPHICS = sin(PHIS);
    const double CPHICS = cos(PHIS);

    const double XS = RHO*CPHICS;
    const double ZS = -RHO*SPHICS;

    double BXS = 0.;
    double BYAS = 0.;
    double BZS = 0.;
    if(*NUMB == 1){
        if(*MODE == 1){
            twocones(A11, &XS, &Ysc, &ZS, MODE, &DTHETA, &BXS, &BYAS, &BZS);
        } else if(*MODE == 2){
            twocones(A12, &XS, &Ysc, &ZS, MODE, &DTHETA, &BXS, &BYAS, &BZS);
        }
    } else {
        if(*MODE == 1){
            twocones(A21, &XS, &Ysc, &ZS, MODE, &DTHETA, &BXS, &BYAS, &BZS);
        } else if(*MODE == 2){
            twocones(A22, &XS, &Ysc, &ZS, MODE, &DTHETA, &BXS, &BYAS, &BZS);
        }
    }

    const double BRHOAS = BXS*CPHICS - BZS*SPHICS;
    const double BPHIAS = -BXS*SPHICS - BZS*CPHICS;

    const double BRHO_S = BRHOAS*DPHISPHI**XKAPPA;
    const double BPHI_S = (BPHIAS - RHO*(BYAS*DPHISDY + BRHOAS*DPHISRHO))**XKAPPA;
    const double BY_S = BYAS*DPHISPHI**XKAPPA;

    *BX = BRHO_S*CPHIC - BPHI_S*SPHIC;
    *BY = BY_S;
    *BZ = -BRHO_S*SPHIC - BPHI_S*CPHIC;

}

void twocones(const float *A, const double *X, const double *Y, const double *Z, const int *MODE, const double *DTHETA,
              double *BX, double *BY, double *BZ){

    double BXN, BYN, BZN;
    one_cone(A, X, Y, Z, MODE, DTHETA, &BXN, &BYN, &BZN);

    const double YY = -*Y;
    const double ZZ = -*Z;
    double BXS, BYS, BZS;
    one_cone(A, X, &YY, &ZZ, MODE, DTHETA, &BXS, &BYS, &BZS);

    *BX = BXN - BXS;
    *BY = BYN + BYS;
    *BZ = BZN + BZS;

}

void one_cone(const float *A, const double *X, const double *Y, const double *Z, const int *MODE, const double *DTHETA,
              double *BX, double *BY, double *BZ){

    const double DR = 1E-6;
    const double DT = 1E-6;

    const double THETA0 = A[30];

    const double RHO2 = custom_pow(*X, 2) + custom_pow(*Y, 2);
    const double RHO = sqrt(RHO2);
    const double R = sqrt(RHO2 + custom_pow(*Z, 2));
    const double THETA = atan2(RHO,*Z);
    const double PHI = atan2(*Y,*X);

    const double RS = r_s(A,&R,&THETA);
    const double THETAS = theta_s(A,&R,&THETA);
    const double PHIS = PHI;

    double BTAST, BFAST;
    fialcos(&RS, &THETAS, &PHIS, MODE, &THETA0, DTHETA, &BTAST, &BFAST);

    const double RSR = RS/R;
    const double STSST = sin(THETAS)/sin(THETA);

    const double R_DR_PLUS = R + DR;
    const double R_DR_MINUS = R - DR;
    const double DRSDR = (r_s(A,&R_DR_PLUS,&THETA) - r_s(A,&R_DR_MINUS,&THETA))/(2.*DR);

    const double THETA_DT_PLUS = THETA + DT;
    const double THETA_DT_MINUS = THETA - DT;
    const double DRSDT = (r_s(A,&R,&THETA_DT_PLUS) - r_s(A,&R,&THETA_DT_MINUS))/(2.*DT);

    const double DTSDR = (theta_s(A,&R_DR_PLUS,&THETA) - theta_s(A,&R_DR_MINUS,&THETA))/(2.*DR);
    const double DTSDT = (theta_s(A,&R,&THETA_DT_PLUS) - theta_s(A,&R,&THETA_DT_MINUS))/(2.*DT);

    const double BR = -RSR/R*STSST*BTAST*DRSDT;
    const double BTHETA = RSR*STSST*BTAST*DRSDR;
    const double BPHI = RSR*BFAST*(DRSDR*DTSDT - DRSDT*DTSDR);

    const  double S = RHO/R;
    const double C = *Z/R;
    const double SF = *Y/RHO;
    const double CF = *X/RHO;
    const double BE = BR*S + BTHETA*C;

    *BX = A[0]*(BE*CF - BPHI*SF);
    *BY = A[0]*(BE*SF + BPHI*CF);
    *BZ = A[0]*(BR*C - BTHETA*S);

}

double r_s(const float *A, const double *R, const double *THETA){
    const double R2 = custom_pow(*R, 2);
    return *R +
           A[1] / *R + A[2]**R/ sqrt(R2 + custom_pow(A[10], 2)) + A[3]**R/(R2 + custom_pow(A[11], 2)) +
           (A[4] + A[5] / *R + A[6]**R/sqrt(R2 + custom_pow(A[12],2)) + A[7]**R/(R2 + custom_pow(A[13], 2)))*cos(*THETA) +
           (A[8]**R/sqrt(R2 + custom_pow(A[14],2)) + A[9]**R/custom_pow((R2 + custom_pow(A[15],2)),2))*cos(2.**THETA);
}

double theta_s(const float *A, const double *R, const double *THETA){
    const double R2 = custom_pow(*R, 2);
    return *THETA +
           (A[16] + A[17] / *R + A[18]/ R2 + A[19]**R/sqrt(R2 + custom_pow(A[26], 2)))*sin(*THETA) +
           (A[20] + A[21]**R/sqrt(R2 + custom_pow(A[27],2)) + A[22]**R/(R2 + custom_pow(A[28], 2)))*sin(2.**THETA) +
           (A[23] + A[24] / *R + A[25]**R/(R2 + custom_pow(A[29], 2)))*sin(3.**THETA);
}

void fialcos(const double *R, const double *THETA, const double *PHI, const int *MODE, const double *THETA0,
             const double *DT, double *BTHETA, double *BPHI){

    double BTN[10], BPN[10], CCOS[10], SSIN[10];

    const double SINTE = sin(*THETA);
    const double RO = *R*SINTE;
    const double COSTE = cos(*THETA);
    const double SINFI = sin(*PHI);
    const double COSFI = cos(*PHI);
    const double TG = SINTE/(1. + COSTE);
    const double CTG = SINTE/(1. - COSTE);

    const double TETANP = *THETA0 + *DT;
    const double TETANM = *THETA0 - *DT;
    double TGP = 0.;
    double TGM = 0.;
    double TGM2 = 0.;
    double TGP2 = 0.;
    if(*THETA >= TETANM){
        TGP = tan(TETANP*0.5);
        TGM = tan(TETANM*0.5);
        TGM2 = TGM*TGM;
        TGP2 = TGP*TGP;
    }

    double COSM1 = 1.;
    double SINM1 = 0.;
    double TM = 1.;
    double TGM2M = 1.;
    double TGP2M = 1.;
    for(int i=1; i <= *MODE; i++){
        TM = TM*TG;
        CCOS[i] = COSM1*COSFI - SINM1*SINFI;
        SSIN[i] = SINM1*COSFI + COSM1*SINFI;
        COSM1 = CCOS[i];
        SINM1 = SSIN[i];
        double T, DTT, FC, FC1, TGM2M1, TG21;
        //double DTT0; // UNUSED
        if(*THETA < TETANM){
            T = TM;
            DTT = 0.5*i*TM*(TG + CTG);
            //DTT0 = 0.; // UNUSED
        } else if(*THETA < TETANP){
            TGM2M = TGM2M*TGM2;
            FC = 1./(TGP - TGM);
            FC1 = 1./(2*i +1);
            TGM2M1 = TGM2M*TGM;
            TG21 = 1. + TG*TG;
            T = FC*(TM*(TGP - TG) + FC1*(TM*TG - TGM2M1/TM));
            DTT = 0.5*i*FC*TG21*(TM/TG*(TGP - TG) - FC1*(TM - TGM2M1/(TM*TG)));
            //DTT0 = 0.5*FC*((TGP + TGM)*(TM*TG - FC1*(TM*TG - TGM2M1/TM)) + TM*(1. - TGP*TGM) - (1. + TGM2)*TGM2M/TM); // UNUSED
        } else {
            TGP2M = TGP2M*TGP2;
            TGM2M = TGM2M*TGM2;
            FC = 1./(TGP - TGM);
            FC1 = 1./(2*i + 1);
            T = FC*FC1*(TGP2M*TGP - TGM2M*TGM)/TM;
            DTT = -T*i*0.5*(TG + CTG);
        }

        BTN[i] = *MODE*T*CCOS[i]/RO;
        BPN[i] = -DTT*SSIN[i] / *R;
    }

    *BTHETA = BTN[*MODE]*800.;
    *BPHI = BPN[*MODE]*800.;

}

void birk_shl(const float *A, const float *PS, const double *X_SC, const double *X, const double *Y, const double *Z,
              double *BX, double *BY, double *BZ){

    const double SPS = sinf(*PS);
    const double CPS = cosf(*PS);
    const double S3PS = 2.*CPS;

    const double PST1 = *PS*A[84];
    const double PST2 = *PS*A[85];
    
    const double ST1 = sin(PST1);
    const double CT1 = cos(PST1);
    const double ST2 = sin(PST2);
    const double CT2 = cos(PST2);

    const double X2 = *X*CT2 - *Z*ST2;    
    const double X1 = *X*CT1 - *Z*ST1;
    const double Z1 = *X*ST1 + *Z*CT1;
    const double Z2 = *X*ST2 + *Z*CT2;

    int L = 3;
    double GX = 0., GY = 0., GZ = 0.;

    for(int j = 1; j <= 3; j++){
        const double P = A[71+j];
        const double YP = *Y / P;
        const double SYPI = sin(YP);
        const double CYPI = cos(YP);
        const double P_SQPR = 1./ (P * P);

        for(int k = 1; k <= 3; k++){
            const double R = A[74+k];
            const double Z1R = Z1 / R;
            const double SZRK = sin(Z1R);
            const double CZRK = cos(Z1R);
            const double SQPR = sqrt(P_SQPR + 1./ (R * R));
            const double EPR = exp(X1*SQPR);

            const double FX = -SQPR*EPR*CYPI*SZRK;
            const double FY = EPR*SYPI*SZRK/P;
            const double FZ = -EPR*CYPI*CZRK/R;

            double HX = FX;
            double HY = FY;
            double HZ = FZ;

            double HXR = HX*CT1 + HZ*ST1;
            double HZR = -HX*ST1 + HZ*CT1;
            GX += HXR * A[L - 3];
            GY += HY * A[L - 3];
            GZ += HZR * A[L - 3];

            HX = FX * *X_SC;
            HY = FY * *X_SC;
            HZ = FZ * *X_SC;

            HXR = HX*CT1 + HZ*ST1;
            HZR = -HX*ST1 + HZ*CT1;
            GX += HXR * A[L - 2];
            GY += HY * A[L - 2];
            GZ += HZR * A[L - 2];

            HX = FX * CPS;
            HY = FY * CPS;
            HZ = FZ * CPS;

            HXR = HX*CT1 + HZ*ST1;
            HZR = -HX*ST1 + HZ*CT1;
            GX += HXR * A[L - 1];
            GY += HY * A[L - 1];
            GZ += HZR * A[L - 1];

            const double CPS_X_SC = CPS * *X_SC;
            HX = FX * CPS_X_SC;
            HY = FY * CPS_X_SC;
            HZ = FZ * CPS_X_SC;

            HXR = HX*CT1 + HZ*ST1;
            HZR = -HX*ST1 + HZ*CT1;
            GX += HXR * A[L];
            GY += HY * A[L];
            GZ += HZR * A[L];

            L += 4;
        }
    }

    for(int j = 1; j <= 3; j++){
        const double Q = A[77+j];
        const double YQ = *Y / Q;
        const double SYQI = sin(YQ);
        const double CYQI = cos(YQ);
        const double Q_SQQS = 1./ (Q * Q);

        for(int k = 1; k <= 3; k++){
            const double S = A[80+k];
            const double Z2S = Z2 / S;
            const double SZSK = sin(Z2S);
            const double CZSK = cos(Z2S);
            const double SQQS = sqrt(Q_SQQS + 1./ (S * S));
            const double EQS = exp(X2*SQQS);

            const double FX = -SPS*SQQS*EQS*CYQI*CZSK;
            const double FY = SPS/Q*EQS*SYQI*CZSK;
            const double FZ = SPS/S*EQS*CYQI*SZSK;

            double HX = FX;
            double HY = FY;
            double HZ = FZ;

            double HXR = HX*CT2 + HZ*ST2;
            double HZR = -HX*ST2 + HZ*CT2;
            GX += HXR * A[L - 3];
            GY += HY * A[L - 3];
            GZ += HZR * A[L - 3];

            HX = FX * *X_SC;
            HY = FY * *X_SC;
            HZ = FZ * *X_SC;

            HXR = HX*CT2 + HZ*ST2;
            HZR = -HX*ST2 + HZ*CT2;
            GX += HXR * A[L - 2];
            GY += HY * A[L - 2];
            GZ += HZR * A[L - 2];

            HX = FX * S3PS;
            HY = FY * S3PS;
            HZ = FZ * S3PS;

            HXR = HX*CT2 + HZ*ST2;
            HZR = -HX*ST2 + HZ*CT2;
            GX += HXR * A[L - 1];
            GY += HY * A[L - 1];
            GZ += HZR * A[L - 1];

            const double S3PS_X_SC = S3PS * *X_SC;
            HX = FX * S3PS_X_SC;
            HY = FY * S3PS_X_SC;
            HZ = FZ * S3PS_X_SC;
            
            HXR = HX*CT2 + HZ*ST2;
            HZR = -HX*ST2 + HZ*CT2;
            GX += HXR * A[L];
            GY += HY * A[L];
            GZ += HZR * A[L];

            L += 4;
        }
    }

    *BX = GX;
    *BY = GY;
    *BZ = GZ;
}