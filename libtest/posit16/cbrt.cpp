#import "mpfr.h"
#import "posit16_math.h"
#import "softposit.h"

#define MPFR_PREC 2000
mpfr_t mval;

double MpfrCalculateSqrt(posit16 x) {
    mpfr_set_d(mval, x.toDouble(), MPFR_RNDN);
    mpfr_sqrt(mval, mval, MPFR_RNDN);
    
    // Check for Nan
    if (mpfr_nan_p(mval) != 0) return castP16(0x8000);
    // Check for infinity
    if (mpfr_inf_p(mval) != 0) return castP16(0x8000);
    // Check for 0
    if (mpfr_cmp_d(mval, 0.0) == 0) return castP16(0x0);
    
    // Otherwise we need to run non-zero values.
    if (mpfr_cmp_d(mval, 0) > 0) {
        if (mpfr_cmp_d(mval, pow(2, 27)) > 0) return castP16(0x7fff);
        if (mpfr_cmp_d(mval, pow(2, 26)) >= 0) return castP16(0x7ffe);
        if (mpfr_cmp_d(mval, pow(2, -27)) < 0) return castP16(0x0001);
        if (mpfr_cmp_d(mval, pow(2, -26)) <= 0) return castP16(0x0002);
    } else {
        if (mpfr_cmp_d(mval, -pow(2, 27)) < 0) return castP16(0x8001);
        if (mpfr_cmp_d(mval, -pow(2, 26)) <= 0) return castP16(0x8002);
        if (mpfr_cmp_d(mval, -pow(2, -27)) > 0) return castP16(0xffff);
        if (mpfr_cmp_d(mval, -pow(2, -26)) >= 0) return castP16(0xfffe);
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

    posit16 x = 0.0;
    for (; count < 0x10000; count++) {
        x.value = count;
        posit16 bres = mysqrt(x);
        posit16 bmy = MpfrCalculateSqrt(x);
        
        // if bres is nan and bmy is nan, continue
        if (bres != bres && bmy != bmy) continue;
        if (bres != bmy) {
            printf("x = %.30e\n", x.toDouble());
            wrongDoubleCount++;
        }
    }
    
    printf("Found %d/%d values that did not calculate correctly\n", wrongDoubleCount, count);
    mpfr_clear(mval);
}
