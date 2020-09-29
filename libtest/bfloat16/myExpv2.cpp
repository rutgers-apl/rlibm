#include "mpfr.h"
#include "bfloat16_math.hpp"
#include "bfloat16.hpp"
#include <math.h>
#include <x86intrin.h>

#define MPFR_PREC 2000
mpfr_t mval;

double MpfrCalculateExp(bfloat16 x) {
    mpfr_set_d(mval, (double)x, MPFR_RNDN);
    mpfr_exp(mval, mval, MPFR_RNDN);
    double retVal = mpfr_get_d(mval, MPFR_RNDN);

    if (retVal == 0) return 0.0;
    if (retVal != retVal) {
        return retVal;
    }
    
    if (mpfr_cmp_d(mval, pow(2, -134)) <= 0 && mpfr_cmp_d(mval, -pow(2, -134)) >= 0) {
        return 0;
    }

    long exp;
    double fr = mpfr_get_d_2exp(&exp, mval, MPFR_RNDN);
    fr *= 2;
    exp--;
    
    if (mpfr_cmp_d(mval, 0.0) > 0) {
        if (mpfr_cmp_d(mval, 1.5 * pow(2, -133)) < 0) return pow(2, -133);
        if (mpfr_cmp_d(mval, pow(2, -132)) < 0) return pow(2, -132);

    } else {
        if (mpfr_cmp_d(mval, -1.5 * pow(2, -133)) > 0) return -pow(2, -133);
        if (mpfr_cmp_d(mval, -pow(2, -132)) > 0) return -pow(2, -132);
    }
    
    if (exp >= -132 && exp <= -127) {
        int prec = 134 + exp;
        mpfr_t r;
        mpfr_init2(r, prec);
        mpfr_set(r, mval, MPFR_RNDN);
        retVal = mpfr_get_d(r, MPFR_RNDN);
        mpfr_clear(r);
        return retVal;
    } else {
        mpfr_t r;
        mpfr_init2(r, 8);
        mpfr_set(r, mval, MPFR_RNDN);
        retVal = mpfr_get_d(r, MPFR_RNDN);
        mpfr_clear(r);
        return retVal;
    }
}

bfloat16 myExpTest(bfloat16 x, unsigned long long& time) {
    unsigned int dummy;
    unsigned long long t1 = __rdtscp(&dummy);
    bfloat16 result = myexpv2(x);
    unsigned long long t2 = __rdtscp(&dummy);
    time += (t2 - t1);
    return result;
}

bfloat16 mlibExpTest(bfloat16 x, unsigned long long& time) {
    unsigned int dummy;
    unsigned long long t1 = __rdtscp(&dummy);
    bfloat16 result = expf((float)x);
    unsigned long long t2 = __rdtscp(&dummy);
    time += (t2 - t1);
    return result;
}

bfloat16 doubleExpTest(bfloat16 x, unsigned long long& time) {
    unsigned int dummy;
    unsigned long long t1 = __rdtscp(&dummy);
    bfloat16 result = exp((double)x);
    unsigned long long t2 = __rdtscp(&dummy);
    time += (t2 - t1);
    return result;
}

int main(int argc, char** argv) {
    mpfr_init2(mval, MPFR_PREC);
    int wrongBfloatCount = 0;
    int wrongFloatCount = 0;
    int wrongDoubleCount = 0;
    unsigned long long count = 0;
    unsigned long long myTime = 0;
    unsigned long long mlibTime = 0;
    unsigned long long doubleTime = 0;

    bfloat16 x = 0.0;
    for (; count < 0x10000; count++) {
        x.val = count;
        bfloat16 bres = myExpTest(x, myTime);
        bfloat16 bmy = MpfrCalculateExp(x);
        bfloat16 bfy = mlibExpTest(x, mlibTime);
        bfloat16 bdy = doubleExpTest(x, doubleTime);
        
        // if bres is nan and bmy is nan, continue
        if (bres != bres && bmy != bmy && bfy != bfy && bdy != bdy) continue;
        if (bres != bmy) {
            printf("x    = %.100e\n", (double)x);
            printf("bres = %.100e (%x)\n", (double)bres, bres.val);
            printf("bmy  = %.100e (%x)\n", (double)bmy, bmy.val);
            wrongBfloatCount++;
        }
        if (bfy != bmy) wrongFloatCount++;
        if (bdy != bmy) wrongDoubleCount++;
    }
    
    printf("Found %d/%llu values that did not calculate correctly\n", wrongBfloatCount, count);
    printf("Average time = %llu cycles\n", myTime / count);
    printf("Float computes %d/%llu values incorrectly\n", wrongFloatCount, count);
    printf("Average time = %llu cycles\n", mlibTime / count);
    printf("Double computes %d/%llu values incorrectly\n", wrongDoubleCount, count);
    printf("Average time = %llu cycles\n", doubleTime / count);
    mpfr_clear(mval);
}
