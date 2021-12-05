#include <math.h>

#include "positions.h"

void calculate_position_and_speed(const float *h, const float *d, const float *bx, const float *by, const float *bz,
                                  float *vx, float *vy, float *vz, float *x, float *y, float *z, float *vv){
    float a1 = *h*FNX(d,vy,vz,by,bz);
    float b1 = *h*FNY(d,vx,vz,bx,bz);
    float c1 = *h*FNZ(d,vx,vy,bx,by);

    const float y_a2 = *vy + b1/3.f;
    const float z_a2 = *vz + c1/3.f;
    const float a2 = *h*FNX(d,&y_a2,&z_a2,by,bz);
    const float y_b2 = *vx + a1/3.f;
    const float z_b2 = *vz + c1/3.f;
    const float b2 = *h*FNY(d,&y_b2,&z_b2,bx,bz);
    const float y_c2 = *vx + a1/3.f;
    const float z_c2 = *vy + b1/3.f;
    const float c2 = *h*FNZ(d,&y_c2,&z_c2,bx,by);

    const float y_a3 = *vy + (4.f*b1 + 6.f*b2)/25.f;
    const float z_a3 = *vz + (4.f*c1 + 6.f*c2)/25.f;
    const float a3 = *h*FNX(d,&y_a3,&z_a3,by,bz);
    const float y_b3 = *vx + (4.f*a1 + 6.f*a2)/25.f;
    const float z_b3 = *vz + (4.f*c1 + 6.f*c2)/25.f;
    const float b3 = *h*FNY(d,&y_b3,&z_b3,bx,bz);
    const float y_c3 = *vx + (4.f*a1 + 6.f*a2)/25.f;
    const float z_c3 = *vy + (4.f*b1 + 6.f*b2)/25.f;
    const float c3 = *h*FNZ(d,&y_c3,&z_c3,bx,by);

    const float y_a4 = *vy + (b1 - 12.f*b2+15.f*b3)/4.f;
    const float z_a4 = *vz + (c1 - 12.f*c2 + 15.f*c3)/4.f;
    const float a4 = *h*FNX(d,&y_a4,&z_a4,by,bz);
    const float y_b4 = *vx + (a1 - 12.f*a2+15.f*a3)/4.f;
    const float z_b4 = *vz + (c1 - 12.f*c2 + 15.f*c3)/4.f;
    const float b4 = *h*FNY(d,&y_b4,&z_b4,bx,bz);
    const float y_c4 = *vx + (a1 - 12.f*a2+15.f*a3)/4.f;
    const float z_c4 = *vy + (b1 - 12.f*b2 + 15.f*b3)/4.f;
    const float c4 = *h*FNZ(d,&y_c4,&z_c4,bx,by);

    const float y_a5 = *vy + (6.f*b1 + 90.f*b2 - 50.f*b3 + 8.f*b4)/81.f;
    const float z_a5 = *vz + (6.f*c1 + 90.f*c2 - 50.f*c3 + 8.f*c4)/81.f;
    const float a5 = *h*FNX(d,&y_a5,&z_a5,by,bz);
    const float y_b5 = *vx + (6.f*a1 + 90.f*a2 - 50.f*a3 + 8.f*a4)/81.f;
    const float z_b5 = *vz + (6.f*c1 + 90.f*c2 - 50.f*c3 + 8.f*c4)/81.f;
    const float b5 = *h*FNY(d,&y_b5,&z_b5,bx,bz);
    const float y_c5 = *vx + (6.f*a1 + 90.f*a2 - 50.f*a3 + 8.f*a4)/81.f;
    const float z_c5 = *vy + (6.f*b1 + 90.f*b2 - 50.f*b3 + 8.f*b4)/81.f;
    const float c5 = *h*FNZ(d,&y_c5,&z_c5,bx,by);

    const float y_a6 = *vy + (6.f*b1 + 36.f*b2 + 10.f*b3 + 8.f*b4)/75.f;
    const float z_a6 = *vz + (6.f*c1 + 36.f*c2 + 10.f*c3 + 8.f*c4)/75.f;
    const float a6 = *h*FNX(d,&y_a6,&z_a6,by,bz);
    const float y_b6 = *vx + (6.f*a1 + 36.f*a2 + 10.f*a3 + 8.f*a4)/75.f;
    const float z_b6 = *vz + (6.f*c1 + 36.f*c2 + 10.f*c3 + 8.f*c4)/75.f;
    const float b6 = *h*FNY(d,&y_b6,&z_b6,bx,bz);
    const float y_c6 = *vx + (6.f*a1 + 36.f*a2 + 10.f*a3 + 8.f*a4)/75.f;
    const float z_c6 = *vy + (6.f*b1 + 36.f*b2 + 10.f*b3 + 8.f*b4)/75.f;
    const float c6 = *h*FNZ(d,&y_c6,&z_c6,bx,by);

    const float aa1 = *h**vx;
    const float bb1 = *h**vy;
    const float cc1 = *h**vz;
    const float aa2 = *h*(*vx + aa1/3.f);
    const float bb2 = *h*(*vy + bb1/3.f);
    const float cc2 = *h*(*vz + cc1/3.f);
    const float aa3 = *h*(*vx + (4.f*aa1 + 6.f*aa2)/25.f);
    const float bb3 = *h*(*vy + (4.f*bb1 + 6.f*bb2)/25.f);
    const float cc3 = *h*(*vz + (4.f*cc1 + 6.f*cc2)/25.f);
    const float aa4 = *h*(*vx + (aa1 - 12.f*aa2 + 15.f*aa3)/4.f);
    const float bb4 = *h*(*vy + (bb1 - 12.f*bb2 + 15.f*bb3)/4.f);
    const float cc4 = *h*(*vz + (cc1 - 12.f*cc2 + 15.f*cc3)/4.f);
    const float aa5 = *h*(*vx + (6.f*aa1 + 90.f*aa2 - 50.f*aa3 + 8.f*aa4)/81.f);
    const float bb5 = *h*(*vy + (6.f*bb1 + 90.f*bb2 - 50.f*bb3 + 8.f*bb4)/81.f);
    const float cc5 = *h*(*vz + (6.f*cc1 + 90.f*cc2 - 50.f*cc3 + 8.f*cc4)/81.f);
    const float aa6 = *h*(*vx + (6.f*aa1 + 36.f*aa2 + 10.f*aa3 + 8.f*aa4)/75.f);
    const float bb6 = *h*(*vy + (6.f*bb1 + 36.f*bb2 + 10.f*bb3 + 8.f*bb4)/75.f);
    const float cc6 = *h*(*vz + (6.f*cc1 + 36.f*cc2 + 10.f*cc3 + 8.f*cc4)/75.f);

    *vx = *vx + (23.f*a1 + 125.f*a3 - 81.f*a5 + 125.f*a6)/192.f;
    *vy = *vy + (23.f*b1 + 125.f*b3 - 81.f*b5 + 125.f*b6)/192.f;
    *vz = *vz + (23.f*c1 + 125.f*c3 - 81.f*c5 + 125.f*c6)/192.f;
    *x = *x + (23.f*aa1 + 125.f*aa3 - 81.f*aa5 + 125.f*aa6)/192.f;
    *y = *y + (23.f*bb1 + 125.f*bb3 - 81.f*bb5 + 125.f*bb6)/192.f;
    *z = *z + (23.f*cc1 + 125.f*cc3 - 81.f*cc5 + 125.f*cc6)/192.f;
    *vv = sqrtf(*vx**vx + *vy**vy + *vz**vz);

}

float FNZ(const float *d, const float *x, const float *y, const float *bx, const float *by){
    return *d*(*x**by - *y**bx);
}

float FNY(const float *d, const float *x, const float *z, const float *bx, const float *bz){
    return *d*(*z**bx - *x**bz);
}

float FNX(const float *d, const float *y, const float *z, const float *by, const float *bz){
    return *d*(*y**bz - *z**by);
}
