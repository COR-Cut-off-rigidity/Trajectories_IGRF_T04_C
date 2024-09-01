//
// Created by pc on 7.4.2020.
//
#include <math.h>

#include "recalc.h"
#include "utility.h"
#include "coefs/igrf_coefs_9.h"
#include "coefs/igrf_coefs_10.h"
#include "coefs/igrf_coefs_11.h"
#include "coefs/igrf_coefs_12.h"
#include "coefs/igrf_coefs_13.h"

void recalc(int igrfCoefsVersion, int iyear, const int *iday, const int *ihour, const int *min, const int *isec,
            float *G, float *H, float *REC, float *A1, float *A2, float *A3) {

    float N2;
    for (int i = 1; i < 15; ++i) {
        N2 = 2*i - 1;
        N2 = N2 * (N2-2);
        int MN = 0, j;
        for (j = 1; j < i+1; ++j) {
            MN = i*(i-1) / 2+j-1;
            REC[MN] = ((i-j)*(i+j-2)) / N2;
        }
    }

    interpolateOrExtrapolate(igrfCoefsVersion, &iyear, iday, G, H);

    float S = 1.f, P, AA;
    int MN = 0, MNN;
    for (int i = 2; i <= 14; ++i) {
        MN = (i*(i-1)/2+1)-1;
        S = S * (2*i-3) / (i-1);
        G[MN] = G[MN] * S;
        H[MN] = H[MN] * S;
        P = S;
        for (int j = 2; j <= i; j++){
            AA = 1.f;
            if (j == 2){
                AA = 2.f;
            }
            P = P * sqrtf(AA * (i-j+1) / (i+j-2));
            MNN = MN+j-1;
            G[MNN] = G[MNN] * P;
            H[MNN] = H[MNN] * P;
        }
    }

    const float G10 = -G[1];
    const float G11 = G[2];
    const float H11 = H[2];

    const float SQ = custom_pow(G11, 2) + custom_pow(H11, 2);
    const float SQQ = sqrtf(SQ);
    const float SQR = sqrtf(custom_powf(G10, 2) + SQ);

    float GST = 0.f, SLONG = 0.f, SRASN = 0.f, SDEC = 0.f;
    sun(&iyear, iday, ihour, min, isec, &GST, &SLONG, &SRASN, &SDEC);

    const float CT0 = G10/SQR;
    const float CL0 = -G11/SQQ;
    const float SL0 = -H11/SQQ;
    const float ST0 = SQQ/SQR;
    //const float CTCL = CT0*CL0; // UNUSED
    //const float CTSL = CT0*SL0; // UNUSED
    const float STCL = ST0*CL0;
    const float STSL = ST0*SL0;
    const float CGST = cosf(GST);
    const float SGST = sinf(GST);
    const float DIP1 = STCL*CGST - STSL*SGST;
    const float DIP2 = STCL*SGST + STSL*CGST;
    const float DIP3 = CT0;

    const float CSDEC = cosf(SDEC);

    const float S1 = cosf(SRASN) * CSDEC;
    const float S2 = sinf(SRASN) * CSDEC;
    const float S3 = sinf(SDEC);
    float Y1 = DIP2*S3 - DIP3*S2;
    float Y2 = DIP3*S1 - DIP1*S3;
    float Y3 = DIP1*S2 - DIP2*S1;
    const float Y = sqrtf(Y1*Y1 + Y2*Y2 + Y3*Y3);
    Y1 = Y1/Y;
    Y2 = Y2/Y;
    Y3 = Y3/Y;
    const float Z1 = S2*Y3 - S3*Y2;
    const float Z2 = S3*Y1 - S1*Y3;
    const float Z3 = S1*Y2 - S2*Y1;

    A1[0] = S1*CGST + S2*SGST;
    A1[1] = -S1*SGST + S2*CGST;
    A1[2] = S3;
    A2[0] = Y1*CGST + Y2*SGST;
    A2[1] = -Y1*SGST + Y2*CGST;
    A2[2] = Y3;
    A3[0] = Z1*CGST + Z2*SGST;
    A3[1] = -Z1*SGST + Z2*CGST;
    A3[2] = Z3;

    // UNUSED VARIABLES
    /*const float DJ = (365*(iyear-1900) + (iyear-1901)/4. + *iday)-0.5f + (*ihour*3600 + *min*60 + *isec)/86400.f;
    const float T = DJ / 36525.f;
    const float OBLIQ = (23.45229f - 0.0130125f*T) / 57.2957795f;
    const float DZ1 = 0.f;
    const float DZ2 = -sinf(OBLIQ);
    const float DZ3 = cosf(OBLIQ);
    const float DY1 = DZ2*S3 - DZ3*S2;
    const float DY2 = DZ3*S1 - DZ1*S3;
    const float DY3 = DZ1*S2 - DZ2*S1;

    const float CHI = Y1*DY1 + Y2*DY2 + Y3*DY3;
    const float SHI = Y1*DZ1 + Y2*DZ2 + Y3*DZ3;
    const float HI = asinf(SHI);

    const float SPS = DIP1*S1 + DIP2*S2 + DIP3*S3;
    const float CPS = sqrtf(1.f- custom_pow(SPS, 2));
    const float PSI = asinf(SPS);

    const float EXMAGX = CT0 * (CL0*CGST - SL0*SGST);
    const float EXMAGY = CT0 * (CL0*SGST + SL0*CGST);
    const float EXMAGZ = -ST0;
    const float EYMAGX = -(SL0*CGST + CL0*SGST);
    const float EYMAGY = -(SL0*SGST - CL0*CGST);
    const float CFI = Y1*EYMAGX + Y2*EYMAGY;
    const float SFI = Y1*EXMAGX + Y2*EXMAGY + Y3*EXMAGZ;

    const float XMUT = (atan2f(SFI,CFI) + 3.1415926536f) * 3.8197186342f;*/

}

void interpolateOrExtrapolate(int igrfCoefsVersion, int *year, const int *iday, float *G, float *H) {

    const float *GHs_9[] = { G1900_9, H1900_9, G1905_9, H1905_9, G1910_9, H1910_9, G1915_9, H1915_9, G1920_9, H1920_9, G1925_9, H1925_9, G1930_9, H1930_9, G1935_9, H1935_9, G1940_9, H1940_9, G1945_9, H1945_9, G1950_9, H1950_9, G1955_9, H1955_9, G1960_9, H1960_9, G1965_9, H1965_9, G1970_9, H1970_9, G1975_9, H1975_9, G1980_9, H1980_9, G1985_9, H1985_9, G1990_9, H1990_9, G1995_9, H1995_9, G2000_9, H2000_9, G2005_9, H2005_9 };
    const float *GHs_10[] = { G1900_9, H1900_9, G1905_9, H1905_9, G1910_9, H1910_9, G1915_9, H1915_9, G1920_9, H1920_9, G1925_9, H1925_9, G1930_9, H1930_9, G1935_9, H1935_9, G1940_9, H1940_9, G1945_9, H1945_9, G1950_9, H1950_9, G1955_9, H1955_9, G1960_9, H1960_9, G1965_9, H1965_9, G1970_9, H1970_9, G1975_9, H1975_9, G1980_9, H1980_9, G1985_9, H1985_9, G1990_9, H1990_9, G1995_9, H1995_9, G2000_9, H2000_9, G2005_10, H2005_10, G2010_10, H2010_10 };
    const float *GHs_11[] = { G1900_9, H1900_9, G1905_9, H1905_9, G1910_9, H1910_9, G1915_9, H1915_9, G1920_9, H1920_9, G1925_9, H1925_9, G1930_9, H1930_9, G1935_9, H1935_9, G1940_9, H1940_9, G1945_9, H1945_9, G1950_9, H1950_9, G1955_9, H1955_9, G1960_9, H1960_9, G1965_9, H1965_9, G1970_9, H1970_9, G1975_9, H1975_9, G1980_9, H1980_9, G1985_9, H1985_9, G1990_9, H1990_9, G1995_9, H1995_9, G2000_9, H2000_9, G2005_11, H2005_11, G2010_11, H2010_11, G2015_11, H2015_11 };
    const float *GHs_12[] = { G1900_9, H1900_9, G1905_9, H1905_9, G1910_9, H1910_9, G1915_9, H1915_9, G1920_9, H1920_9, G1925_9, H1925_9, G1930_9, H1930_9, G1935_9, H1935_9, G1940_9, H1940_9, G1945_9, H1945_9, G1950_9, H1950_9, G1955_9, H1955_9, G1960_9, H1960_9, G1965_9, H1965_9, G1970_9, H1970_9, G1975_9, H1975_9, G1980_9, H1980_9, G1985_9, H1985_9, G1990_9, H1990_9, G1995_9, H1995_9, G2000_9, H2000_9, G2005_11, H2005_11, G2010_12, H2010_12, G2015_12, H2015_12, G2020_12, H2020_12 };
    const float *GHs_13[] = { G1900_9, H1900_9, G1905_9, H1905_9, G1910_9, H1910_9, G1915_9, H1915_9, G1920_9, H1920_9, G1925_9, H1925_9, G1930_9, H1930_9, G1935_9, H1935_9, G1940_9, H1940_9, G1945_9, H1945_9, G1950_9, H1950_9, G1955_9, H1955_9, G1960_9, H1960_9, G1965_9, H1965_9, G1970_9, H1970_9, G1975_9, H1975_9, G1980_9, H1980_9, G1985_9, H1985_9, G1990_9, H1990_9, G1995_9, H1995_9, G2000_9, H2000_9, G2005_11, H2005_11, G2010_12, H2010_12, G2015_13, H2015_13, G2020_13, H2020_13, G2025_13, H2025_13 };
    
    const int GHsNumYears[] = { numYears_9, numYears_10, numYears_11, numYears_12, numYears_13 };
    const float **GHsVersions[] = { GHs_9, GHs_10, GHs_11, GHs_12, GHs_13 };

    const int GHsVersionIndex = igrfCoefsVersion - 9;
    const int numYears = GHsNumYears[GHsVersionIndex];
    const int GHsSize = numYears * 2;
    const float **GHs = GHsVersions[GHsVersionIndex];

    const int minYear = 1900;
    const int maxYear = 1895 + numYears * 5;

    if(*year < minYear) {
        *year = minYear;
    } else if (*year > maxYear){
        *year = maxYear;
    }

    if(maxYear-*year <= 5){
        extrapolate(G, H, GHs[GHsSize-4], GHs[GHsSize-3], GHs[GHsSize-2], GHs[GHsSize-1], maxYear, *year, *iday);
    } else {
        int GHYear;
        if(*year <= minYear){
            GHYear = maxYear / 5 - (minYear + 1) / 5;
        } else {
            GHYear = maxYear / 5 - *year / 5;
        }
        //const int GHPosition = ((GHsSize/2 - GHYear) * 2) - 1;
        const int GHPosition = GHsSize - 2*GHYear - 1;
        const int edgeYear = minYear + ((GHPosition - 1) / 2) * 5;
        interpolate(G, H, GHs[GHPosition-1], GHs[GHPosition], GHs[GHPosition+1], GHs[GHPosition+2], edgeYear, *year, *iday);
    }

}

void interpolate(float *G, float *H, const float *Gx, const float *Hx, const float *Gy, const float *Hy, int edgeYear, int year, int iday){

    const float F2 = (year + (iday-1) / 365.25f - edgeYear) / 5.0f;
    const float F1 = 1.f - F2;

    for (int i = 0; i < 105; ++i) {
        G[i] = Gx[i] * F1 + Gy[i] * F2;
        H[i] = Hx[i] * F1 + Hy[i] * F2;
    }

}

void extrapolate(float *G, float *H, const float *Gx, const float *Hx, const float *Gy, const float *Hy, int maxYear, int year, int iday){

    const float DT = year + (iday-1) / 365.25f - (maxYear-5.0f);

    for (int i = 0; i < 105; ++i) {
        G[i] = Gx[i];
        H[i] = Hx[i];
        if(i > 45) {
            continue;
        }
        G[i] = G[i] + Gy[i] * DT;
        H[i] = H[i] + Hy[i] * DT;
    }

}

void sun(int *iyear, const int *iday, const int *ihour, const int *min, const int *isec,
         float *GST, float *SLONG, float *SRASN, float *SDEC){

    if(*iyear < 1901 || *iyear > 2099){
        return;
    }

    const float RAD = 57.295779513f;
    const double FDAY = (*ihour*3600 + *min*60 + *isec) / 86400.;
    const double DJ = 365*(*iyear-1900) + (int)((*iyear-1901)/4) + *iday - 0.5 + FDAY;
    const float T = DJ/36525.f;
    *GST = fmod(279.690983f + 0.9856473354f * DJ + 360.f * FDAY + 180.f, 360.f) / RAD;

    const float VL = fmod(279.696678f + 0.9856473354f*DJ, 360.f);
    const float G = fmod(358.475845f + 0.985600267f * DJ, 360.f) / RAD;
    *SLONG = (VL + (1.91946f - 0.004789f*T)*sinf(G) + 0.020094f*sinf(2.f*G)) / RAD;
    if(*SLONG > 6.2831853f){
        *SLONG = *SLONG - 6.2831853f;
    } else if(*SLONG < 0){
        *SLONG = *SLONG + 6.2831853f;
    }

    const float OBLIQ = (23.45229f - 0.0130125f*T) / RAD;
    const float SLP = *SLONG-9.924E-5f;
    const float SOB = sinf(OBLIQ);

    const float SIND = SOB*sinf(SLP);
    const float COSD = sqrtf(1.f- custom_pow(SIND, 2));
    const float SC = SIND/COSD;
    *SDEC = atanf(SC);
    *SRASN = 3.141592654f - atan2f(cosf(OBLIQ)/SOB*SC,-cosf(SLP)/COSD);

}
