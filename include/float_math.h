#include <stdint.h>
#include "math.h"
#include <cstdint>

typedef union {
    double d;
    uint64_t x;
} doubleX;

typedef union {
    float f;
    unsigned int x;
} floatX;

float log2small(float);
float mylog2(float);
float mylog2v2(float);

double log2smallInternal(double);

