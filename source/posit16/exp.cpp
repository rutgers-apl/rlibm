#include "posit16_math.h"

posit16_t exp(posit16_t x) {
    
    if (x.v == 0x8000) {
        // Take care of NaR
        return x;
    } else if (x.v > 0x8000 & x.v <= 0x8f52) {
        // Take care of when result is minpos. exp(x) for posit should never
        // return 0, because exp(x) is always > 0 as long as x != -infinity
        x.v = 0x1;
        return x;
    } else if (x.v >= 0x70ae && x.v < 0x8000) {
        // Take care of maxpos case.
        x.v = 0x7FFF;
        return x;
    } else if ((x.v >= 0xff80) || (x.v <= 0xbf)) {
        // The values in these range return 1.0.
        x.v = 0x4000;
        return x;
    }
    
    // Extract exponent and mantissa (where mantissa is between [1, 2)
    double xprime = convertP16ToDouble(x) *
    1.442695040888963387004650940070860087871551513671875;
    double modifier;
    double f = modf(xprime, &modifier);
    if (f < 0.0f) {
        f += 1.0;
        modifier -= 1.0;
    }
    
    // Now compute polynomial
    double y = 3.880795968797916273285153465621988289058208465576171875e-04 ;
    y *= f;
    y += 7.56416228767151587429606873769216690561734139919281005859375e-04;
    y *= f;
    y += 1.019728060416365440776775841413837042637169361114501953125e-02;
    y *= f;
    y += 5.522989323495941516029006379540078341960906982421875e-02;
    y *= f;
    y += 2.40286646236610668125877054990269243717193603515625e-01;
    y *= f;
    y += 6.93141820259879803955982424668036401271820068359375e-01;
    y *= f;
    y += 1.0000001298120355652798707524198107421398162841796875;

    y = ldexp(y, modifier);
    
    return convertDoubleToP16(y);
}

double expInternal(double dx) {
    if (dx <= -18.71875) {
        // Take care of when result is minpos. exp(x) for posit should never
        // return 0, because exp(x) is always > 0 as long as x != -infinity
        return 3.7252902984619140625e-09;
    } else if (dx >= 18.71875) {
        // Take care of maxpos case.
        return 2.68435456e+08;
    } else if (dx != dx) {
        // Take care of NaR
        return dx;
    } else if ((dx >= -6.103515625e-05) && (dx <= 1.2111663818359375e-04)) {
        // The values in these range return 1.0.
        return 1.0;
    }
    
    // Extract exponent and mantissa (where mantissa is between [1, 2)
    double xprime = dx *
    1.442695040888963387004650940070860087871551513671875;
    double modifier;
    double f = modf(xprime, &modifier);
    if (f < 0.0f) {
        f += 1.0;
        modifier -= 1.0;
    }
    
    // Now compute polynomial
    double y = 3.880795968797916273285153465621988289058208465576171875e-04 ;
    y *= f;
    y += 7.56416228767151587429606873769216690561734139919281005859375e-04;
    y *= f;
    y += 1.019728060416365440776775841413837042637169361114501953125e-02;
    y *= f;
    y += 5.522989323495941516029006379540078341960906982421875e-02;
    y *= f;
    y += 2.40286646236610668125877054990269243717193603515625e-01;
    y *= f;
    y += 6.93141820259879803955982424668036401271820068359375e-01;
    y *= f;
    y += 1.0000001298120355652798707524198107421398162841796875;

    y = ldexp(y, modifier);
    
    return y;
}

uint_fast64_t p16_exppoly( uint_fast64_t );

posit16_t p16_exp( posit16_t pA ){

    union ui16_p16 uA;
    uint_fast16_t uiA;
    uint_fast64_t bit, f, s = 0;

    uA.p = pA;
    uiA = uA.ui;
    f = uiA;

    // Calculate the exponential for given posit pA
    if ( uiA < 28846 ) {        // result does not round up to maxpos

        if ( uiA < 192 ) {    // small positive values that round to 1
            uA.ui = 0x4000;
            return uA.p;
        }

        if ( f & 0x4000 ) {    // decode regime
            s = 8;
            while (f & 0x2000) {
                f = f << 1;
                s += 2;
            }
        } else {
            s = 6;
            while ( !(f & 0x2000) ) {
                f = f << 1;
                s -= 2;
            }
        }

        if (f & 0x1000) s++;                    // decode exponent
        f = (f & 0x0FFF) | 0x1000;                // decode fraction
        f = (((int)s < 0 ? f >> -s : f << s) * 48408813) >> 20;
        s = f >> 25;                        // s now stores floor(x)
        f = p16_exppoly(f & 0x1FFFFFF);                // 37 fraction bits of exp(x)
        bit = (s & 1) << 37;                    // exponent bit of exp(x)
        s = s >> 1;                        // regime length of exp(x)
        f = ((0x10000000000 << s) - 0x8000000000) | bit | f;

        bit = (uint_fast64_t)1 << (24+s);            // location of bit n-plus-1
        if ( f & bit ) {
            if ((f & (bit - 1)) || (f & (bit << 1))) {
                f += bit;
            }
        }
        uA.ui = f >> (25 + s);
        return uA.p;        // return rounded exp(x) as posit

    } else if (uiA > 36690) {    // result does not round up to minpos

        if (uiA > 65407) {    // small negative values that round to 1
            uA.ui = 0x4000;
            return uA.p;
        }

        if (f & 0x4000) {    // decode regime
            s = 7;
            while (f & 0x2000) {
                f = f << 1;
                s -= 2;
            }
        } else {
            s = 9;
            while (!(f & 0x2000)) {
                f = f << 1;
                s += 2;
            }
        }

        if (f & 0x1000) s--;                // decode exponent
        f = (f & 0x0FFF) | 0x1FFE000;            // decode fraction
        f = ((int)s < 0) ? (f >> -s) | (0x2000000 - (1 << (13 + s))) : (f << s) & 0x1ffffff;
        f = (0x4000000000000 - ((0x2000000 - f) * 48408813)) >> 20;
        s = (f >> 25) - 32;                // s now stores floor(x)
        f = p16_exppoly(f & 0x1FFFFFF);            // 37 fraction bits of exp(x)
        bit = (s & 1) << 37;                // exponent bit of exp(x)
            s = (-1 - s) >> 1;
        f = 0x4000000000 | bit | f;            // Install regime end bit

        bit = (uint_fast64_t)1 << (24+s);        // location of bit n-plus-1
        if ( f & bit ) {
            if ((f & (bit - 1)) || (f & (bit << 1))) {
                f += bit;
            }
        }
        uA.ui = f >> (25 + s);
        return uA.p;    // return rounded exp(x) as posit

    }

    // Section for exception cases
    if (uiA < 0x8000) {
        uA.ui = 0x7FFF;
        return uA.p;        // return maxpos
    } else if (uiA > 0x8000) {
        uA.ui = 0x0001;
        return uA.p;        // return minpos
    } else {
        uA.ui = 0x8000;
        return uA.p;        // return NaR
    }

}


uint_fast64_t p16_exppoly( uint_fast64_t f ){

    uint_fast64_t s = 0;

    s = (f * 7529) >> 26;
    s = (f * (20487 + s)) >> 20;
    s = (f * (0x4F8300 + s)) >> 24;
    s = (f * (0x38CC980 + s)) >> 20;
    s = (f * (0x1EBFFC800 + s)) >> 26;
    s = ((f * (0x2C5C83600 + s)) >> 22) + 2048;

    return s;
}
