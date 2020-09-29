#include "posit16_math.h"

posit16_t exp2(posit16_t x) {
    
    if (x.v > 0x8000 & x.v <= 0x8d3f) {
        // Take care of when result is minpos. exp(x) for posit should never
        // return 0, because exp(x) is always > 0 as long as x != -infinity
        x.v = 0x1;
        return x;
    } else if (x.v >= 0x72c1 && x.v < 0x8000) {
        // Take care of maxpos case.
        x.v = 0x7FFF;
        return x;
    } else if (x.v == 0x8000) {
        // Take care of NaR
        return x;
    } else if ((x.v >= 0xff64) || (x.v <= 0xdc)) {
        // The values in these range return 1.0.
        x.v = 0x4000;
        return x;
    }
    
    // Extract exponent and mantissa (where mantissa is between [1, 2)
    double xprime = convertP16ToDouble(x);
    double modifier;
    double f = modf(xprime, &modifier);
    if (f < 0.0f) {
        f += 1.0;
        modifier -= 1.0;
    }
    
    // Now compute polynomial
    double y = 3.65720720951589702875372811519127935753203928470611572265625e-04 ;
    y *= f;
    y += 7.694531242965385771370723233530952711589634418487548828125e-04;
    y *= f;
    y += 1.02531888584330328761939910009459708817303180694580078125e-02;
    y *= f;
    y += 5.5157921920287643346991757198338746093213558197021484375e-02;
    y *= f;
    y += 2.40314607704837424062560558013501577079296112060546875e-01;
    y *= f;
    y += 6.9313911175081532878294865440693683922290802001953125e-01;
    y *= f;
    y += 1.0;

    y = ldexp(y, modifier);
    
    return convertDoubleToP16(y);
}

double exp2Internal(double dx) {
    if (dx <= -2.7015625e+01) {
        // Take care of when result is minpos. exp(x) for posit should never
        // return 0, because exp(x) is always > 0 as long as x != -infinity
        return 3.7252902984619140625e-09;
    } else if (dx >= 27.015625) {
        // Take care of maxpos case.
        return 2.68435456e+08;
    } else if (dx != dx) {
        // Take care of NaR
        return dx;
    } else if ((dx >= -8.7738037109375e-05) && (dx <= 1.7547607421875e-04)) {
        // The values in these range return 1.0.
        return 1.0;
    }
    
    // Extract exponent and mantissa (where mantissa is between [1, 2)
    double xprime = dx;
    double modifier;
    double f = modf(xprime, &modifier);
    if (f < 0.0f) {
        f += 1.0;
        modifier -= 1.0;
    }
    
    // Now compute polynomial
    double y = 3.65720720951589702875372811519127935753203928470611572265625e-04 ;
    y *= f;
    y += 7.694531242965385771370723233530952711589634418487548828125e-04;
    y *= f;
    y += 1.02531888584330328761939910009459708817303180694580078125e-02;
    y *= f;
    y += 5.5157921920287643346991757198338746093213558197021484375e-02;
    y *= f;
    y += 2.40314607704837424062560558013501577079296112060546875e-01;
    y *= f;
    y += 6.9313911175081532878294865440693683922290802001953125e-01;
    y *= f;
    y += 1.0;

    y = ldexp(y, modifier);
    
    return y;
}

uint_fast64_t p16_exp2poly( uint_fast64_t );

posit16_t p16_exp2( posit16_t pA ){

    union ui16_p16 uA;
    uint_fast16_t uiA;
    uint_fast64_t bit, f, s = 0;

    uA.p = pA;
    uiA = uA.ui;
    f = uiA;

    // Calculate the exponential for given posit pA
    if ( uiA < 29377 ) {                // result does not round up to maxpos

        if ( uiA < 221 ) {            // cases that round down to 1.
            uA.ui = 0x4000;
            return uA.p;
        }

        if ( f & 0x4000 ) {            // decode regime
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

        if (f & 0x1000) s++;            // decode exponent
        f = (f & 0x0FFF) | 0x1000;        // decode fraction
        f = ((int)s < 0 ? f >> -s : f << s);
        s = f >> 20;                // s now stores floor(x)
        f = p16_exp2poly(f & 0xFFFFF);        // fraction bits of exp2(x)
        bit = (s & 1) << 26;            // exponent bit of exp2(x)
        s = s >> 1;                // regime length of exp2(x)
        f = (((uint_fast64_t)0x20000000 << s) - 0x10000000) | bit | f;

        bit = (uint_fast64_t)1 << (13 + s);    // location of bit n-plus-1
        if ( f & bit ) {
            if ((f & (bit - 1)) || (f & (bit << 1))) {
                f += bit;
            }
        }
        uA.ui = f >> (14 + s);
        return uA.p;                // return rounded exp2(x) as posit

    } else if (uiA > 36159) {

        if (uiA > 65379 ) {            // cases that round up to 1.
            uA.ui = 0x4000;
            return uA.p;
        }
        if (f & 0x4000) {            // decode regime
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

        if (f & 0x1000) s--;            // decode exponent
        f = (f & 0xFFF) | 0x1FFE000;        // decode fraction
        f = ((int)s < 0) ? (f >> -s) | (0x2000000 - (1 << (13 + s))) : (f << s) & 0x1ffffff;
        s = (f >> 20) - 32;            // s now stores floor(x)
        f = p16_exp2poly(f & 0xFFFFF);        // fraction bits of exp2(x)
        bit = (s & 1) << 26;            // exponent bit of exp2(x)
        s = (-1 - s) >> 1;
        f = 0x8000000 | bit | f;        // Install regime end bit

        bit = (uint_fast64_t)1 << (13 + s);    // location of bit n-plus-1
        if ( f & bit ) {
            if ((f & (bit - 1)) || (f & (bit << 1))) {
                f += bit;
            }
        }
        uA.ui = f >> (14 + s);
        return uA.p;            // return rounded exp2(x) as posit

    }

    // Section for exception cases
    if ( uiA < 0x8000 ) {            // cases that round to maxpos
        uA.ui = 0x7FFF;
        return uA.p;
    } else if ( uiA > 0x8000 ) {        // cases that round to minpos
        uA.ui = 0x0001;
        return uA.p;
    } else {
        uA.ui = 0x8000;            // NaR case
        return uA.p;
    }

}


uint_fast64_t p16_exp2poly( uint_fast64_t f ){

    uint_fast64_t s = 0;

    s = (f * (0x9BA00000 + (f * 491))) >> 34;
    s = (f * (0x13F840 + s)) >> 20;
    s = (f * (0x718A80 + s)) >> 16;
    s = (f * (0x1EC04000 + s)) >> 21;
    s = ((f * (0x2C5C8000 + s)) >> 24);

    return s;
}
