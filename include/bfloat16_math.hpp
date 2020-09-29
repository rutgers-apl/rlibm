#include <stdint.h>
#include "bfloat16.hpp"
#include <math.h>
#include <cstdint>

typedef union {
    double d;
    uint64_t x;
} dx;

typedef union {
    float f;
    unsigned int x;
} fx;

bfloat16 myexpv2(bfloat16);
bfloat16 myexp2v2(bfloat16);
bfloat16 myexp10v2(bfloat16);
bfloat16 mylog(bfloat16);
bfloat16 mylog2(bfloat16);
bfloat16 mylog10(bfloat16);
bfloat16 mysinpi(bfloat16);
bfloat16 mycospi(bfloat16);
bfloat16 mysqrt(bfloat16);
bfloat16 mycbrt(bfloat16);

double myexpInternalv2(float);
double myexp2Internalv2(float);
double myexp10Internalv2(float);
double mylogInternal(float);
double mylog2Internal(float);
double mylog10Internal(float);
double mysinpiInternal(float);
double mycospiInternal(float);
double mysqrtInternal(float);
double mycbrtInternal(float);
