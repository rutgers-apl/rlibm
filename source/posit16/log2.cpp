#include "posit16_math.h"

posit16_t log2(posit16_t x) {
    if (x.v == 0x0 || x.v >= 0x8000) {
        // If x == 0, NaR, or negative, then resutl should be NaR
        x.v = 0x8000;
        return x;
    }
    
    // Extract exponent and mantissa (where mantissa is between [1, 2))
    int m;
    double fx = frexp(convertP16ToDouble(x), &m);
    fx *= 2.0;
    m--;
    
    // Cody and Waite Transformation on input
    double dx = fx;
    double codyX = (dx - 1) / (dx + 1);

    double codyX2 = codyX * codyX;
    // Now compute polynomial
    double y = 6.8022527114824737903830964569351635873317718505859375e-01;
    y *= codyX2;
    y += 3.36729567454907396939489672149647958576679229736328125e-01;
    y *= codyX2;
    y += 5.827497609092706642996972732362337410449981689453125e-01;
    y *= codyX2;
    y += 9.616405555684151007511673014960251748561859130859375e-01;
    y *= codyX2;
    y += 2.88539115994917327867597123258747160434722900390625;
    y *= codyX;
    y += (double)m;
    
    return convertDoubleToP16(y);
}

double log2Internal(double dx) {
    dint ch;
    ch.d = dx;
    if (ch.x > 0x7FF0000000000000 || ch.x == 0) {
        // If x == 0, NaR, or negative, then resutl should be NaR
        return 0.0;
    }
    
    // Extract exponent and mantissa (where mantissa is between [1, 2))
    int m;
    double fx = frexp(dx, &m);
    fx *= 2.0;
    m--;
    
    // Cody and Waite Transformation on input
    double codyX = (fx - 1) / (fx + 1);

    double codyX2 = codyX * codyX;
    // Now compute polynomial
    double y = 6.8022527114824737903830964569351635873317718505859375e-01;
    y *= codyX2;
    y += 3.36729567454907396939489672149647958576679229736328125e-01;
    y *= codyX2;
    y += 5.827497609092706642996972732362337410449981689453125e-01;
    y *= codyX2;
    y += 9.616405555684151007511673014960251748561859130859375e-01;
    y *= codyX2;
    y += 2.88539115994917327867597123258747160434722900390625;
    y *= codyX;
    y += (double)m;
    
    return y;
}

uint_fast64_t p16_log2poly( uint_fast64_t );

posit16_t p16_log2( posit16_t pA ) {

    union ui16_p16 uA;
    uint_fast16_t uiA;
    uint_fast64_t bit, f, sign;
    int_fast64_t s;            // s can be negative
    uA.p = pA;
    uiA = uA.ui;
    f = uiA;

    if ( (f > 0x7FFF) || !f ) {     // if input is greater than maxpos or 0, return NaR
        uA.ui = 0x8000;
        return uA.p;
    }

    if ( f & 0x4000 ) {        // decode regime
        s = 0;
        while ( f & 0x2000 ) {
            f = f << 1;
            s += 2;
        }
    } else {
        s = -2;
        while  ( !(f & 0x2000) ) {
            f = f << 1;
            s -= 2;
        }
    }

    if ( f & 0x1000 ) s++;        // decode exponent
    f = f & 0xFFF;            // get 12-bit fraction, without hidden bit
    if (f) f = p16_log2poly(f);    // turn fraction into mantissa of logarithm
    f = (((s < 0) ? 64 + s : s) << 28) | f;
    sign = f & 0x200000000;
    if (sign) f = 0x400000000 - f;    // take absolute value of fixed-point result
    if (f < 0x10000000) {        // turn fixed-point into posit format
        if (f) {
            s = 30;
            while ( !(f & 0x8000000) ) {
                f = f << 1;
                s++;
            }
            f = (f ^ 0x18000000) | (( 1 ^ (s & 1)) << 27);
            s = s >> 1;
            bit = 1 << (s - 1);
            if ( f & bit ) {
                if ((f & (bit - 1)) || (f & (bit << 1))) f += bit;
            }
            f = f >> s;
        }
    } else {
        s = 0;
        while ( f > 0x1FFFFFFF ) {
            f = (f & 1) | (f >> 1);
            s++;
        }
        f = f & 0xFFFFFFF;
        if (s & 1) f = f | 0x10000000;
        s = s >> 1;
        f = (( (uint_fast64_t)0x80000000 << s) - 0x40000000) | f;
        bit = 0x8000 << s;
        if ( f & bit ) {
            if ((f & (bit - 1)) || (f & (bit << 1))) f += bit;
        }
        f = f >> (s + 16);
    }
    if (sign) f = 0x10000 - f;    // restore sign
    uA.ui = f;
    return uA.p;
}


uint_fast64_t p16_log2poly( uint_fast64_t f ) {

    uint_fast64_t s, z, zsq;

    z = (f << 29) / (f + 8192);        // fixed-point divide; discard remainder
    zsq = (( z*z ) >> 30);            // fixed-point squaring
    s = (zsq * 1661) >> 25;
    s = (zsq * (13209 + s)) >> 26;
    s = (zsq * (75694 + s)) >> 24;
    s = (zsq * (2017019 + s)) >> 24;
    s = (z   * (96817627 + s)) >> 26;

    return s;

}
