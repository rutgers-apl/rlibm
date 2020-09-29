#include "mpfr.h"
#include "posit16_math.h"
#include "softposit.h"

#define MPFR_PREC 2000

mpfr_t mval;

posit16_t MpfrCalculateExp10(posit16_t x) {
    mpfr_set_d(mval, convertP16ToDouble(x), MPFR_RNDN);
    mpfr_exp10(mval, mval, MPFR_RNDN);
    
    // Check for Nan
    if (mpfr_nan_p(mval) != 0) return castP16(0x8000);
    // Check for infinity
    if (mpfr_inf_p(mval) != 0) return castP16(0x8000);
    // Check for 0
    if (mpfr_cmp_d(mval, 0.0) == 0) return castP16(0x0);
    
    if (mpfr_cmp_d(mval, 0) > 0) {
        if (mpfr_cmp_d(mval, pow(2, 27)) > 0) return castP16(0x7fff);
        if (mpfr_cmp_d(mval, 1.5 * pow(2, 25)) >= 0) return castP16(0x7ffe);
        if (mpfr_cmp_d(mval, 1.5 * pow(2, 24)) > 0) return castP16(0x7ffd);
        if (mpfr_cmp_d(mval, pow(2, 24)) >= 0) return castP16(0x7ffc);
        
        
        if (mpfr_cmp_d(mval, pow(2, -27)) < 0) return castP16(0x0001);
        if (mpfr_cmp_d(mval, 1.5 * pow(2, -26)) <= 0) return castP16(0x0002);
        if (mpfr_cmp_d(mval, 1.5 * pow(2, -25)) < 0) return castP16(0x0003);
        if (mpfr_cmp_d(mval, pow(2, -24)) <= 0) return castP16(0x0004);
    } else {
        if (mpfr_cmp_d(mval, -pow(2, 27)) < 0) return castP16(0x8001);
        if (mpfr_cmp_d(mval, -1.5 * pow(2, 25)) <= 0) return castP16(0x8002);
        if (mpfr_cmp_d(mval, -1.5 * pow(2, 24)) < 0) return castP16(0x8003);
        if (mpfr_cmp_d(mval, -pow(2, 24)) <= 0) return castP16(0x8004);
        
        if (mpfr_cmp_d(mval, -pow(2, -27)) > 0) return castP16(0xffff);
        if (mpfr_cmp_d(mval, -1.5 * pow(2, -26)) >= 0) return castP16(0xfffe);
        if (mpfr_cmp_d(mval, -1.5 * pow(2, -25)) > 0) return castP16(0xfffd);
        if (mpfr_cmp_d(mval, -pow(2, -24)) >= 0) return castP16(0xfffc);
    }
    
    long exp;
    double fr = mpfr_get_d_2exp(&exp, mval, MPFR_RNDN);
    long origExp = exp;
    fr *= 2;
    exp--;
    if (exp < 0) {
        exp *= -1;
        exp--;;
    }
    exp >>= 1;
    long p = 13 - exp;
    mpfr_t r;
    mpfr_init2(r, p);
    mpfr_set(r, mval, MPFR_RNDN);
    double retVal = mpfr_get_d(r, MPFR_RNDN);
    mpfr_clear(r);
    return convertDoubleToP16(retVal);
}

int main(int argc, char** argv) {
    mpfr_init2(mval, MPFR_PREC);
    int wrongDoubleCount = 0;
    int count = 0;

    for (; count < 0x10000; count++) {
        posit16_t x = castP16(count);
        posit16_t bres = exp10(x);
        posit16_t bmy = MpfrCalculateExp10(x);
        
        if (!p16_eq(bres, bmy)) wrongDoubleCount++;
    }
    
    printf("Found %d/%d values that did not calculate correctly\n", wrongDoubleCount, count);
    
    mpfr_clear(mval);
}
