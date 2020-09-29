#include "posit16_math.h"

posit16_t log(posit16_t x) {
    if (x.v == 0x0 || x.v >= 0x8000) {
        // If x == 0, NaR, or negative, then resutl should be NaR
        x.v = 0x8000;
        return x;
    }
    
    // Extract exponent and mantissa (where mantissa is between [1, 2)
    int m;
    double fx = frexp(convertP16ToDouble(x), &m);
    fx *= 2.0;
    m--;
    
    // Cody and Waite Transformation on input
    double dx = fx;
    double codyX = (dx - 1) / (dx + 1);
    double codyX2 = codyX * codyX;
    
    // Now compute polynomial
    double y = 4.5254178489671204044242358577321283519268035888671875e-01;
    y *= codyX2;
    y += 3.9449243216490248453709455134230665862560272216796875e-01;
    y *= codyX2;
    y += 5.7802192858859535729010303839459083974361419677734375e-01;
    y *= codyX2;
    y += 9.6177728824005104257821585633791983127593994140625e-01;
    y *= codyX2;
    y += 2.8853901812623536926594169926829636096954345703125;
    y *= codyX;
    
    // Range propagation
    y = (m + y) / 1.442695040888963387004650940070860087871551513671875;
    
    return convertDoubleToP16(y);
}

double logInternal(double dx) {
    dint ch;
    ch.d = dx;
    if (ch.x > 0x7FF0000000000000 || ch.x == 0) {
        // If x == 0, NaR, or negative, then resutl should be NaR
        return 0.0;
    }
    
    // Extract exponent and mantissa (where mantissa is between [1, 2)
    int m;
    double fx = frexp(dx, &m);
    fx *= 2.0;
    m--;
    
    // Cody and Waite Transformation on input
    double codyX = (fx - 1) / (fx + 1);
    double codyX2 = codyX * codyX;
    
    // Now compute polynomial
    double y = 4.5254178489671204044242358577321283519268035888671875e-01;
    y *= codyX2;
    y += 3.9449243216490248453709455134230665862560272216796875e-01;
    y *= codyX2;
    y += 5.7802192858859535729010303839459083974361419677734375e-01;
    y *= codyX2;
    y += 9.6177728824005104257821585633791983127593994140625e-01;
    y *= codyX2;
    y += 2.8853901812623536926594169926829636096954345703125;
    y *= codyX;
    
    // Range propagation
    y = (m + y) / 1.442695040888963387004650940070860087871551513671875;
    
    return y;
}

uint_fast64_t p16_logpoly( uint_fast64_t );

posit16_t p16_log( posit16_t pA ) {

    union ui16_p16 uA;
    uint_fast16_t uiA;
    uint_fast64_t bit, f, sign;
    int_fast64_t s;            // s can be negative
    uA.p = pA;
    uiA = uA.ui;
    f = uiA;

    if ( (f > 0x7FFF) || !f ) {    // if input is 0, or greater than maxpos, return NaR
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
    if (f) f = p16_logpoly(f);    // turn fraction into mantissa of logarithm
    f = (((s < 0) ? 64 + s : s) << 30) | f;

    f = (s < 0) ? 0x1000000000 - (((0x1000000000 - f) * 186065280) >> 28) : ((f * 186065279) >> 28);

    sign = f & 0x800000000;
    if (sign) f = 0x1000000000 - f;    // take absolute value of fixed-point result
    if (f < 0x40000000) {        // turn fixed-point into posit format
        if (f) {
            s = 34;
            while ( !(f & 0x20000000) ) {
                f = f << 1;
                s++;
            }
            f = (f ^ 0x60000000) | (( 1 ^ (s & 1)) << 29);
            s = s >> 1;
            bit = 1 << (s - 1);
            if ( f & bit ) {
                if ((f & (bit - 1)) || (f & (bit << 1))) f += bit;
            }
            f = f >> s;
        }
    } else {
        s = 0;
        while ( f > 0x7FFFFFFF ) {
            f = (f & 1) | (f >> 1);
            s++;
        }
        f = f & 0x3FFFFFFF;
        if (s & 1) f = f | 0x40000000;
        s = s >> 1;
        f = (( (uint_fast64_t)0x200000000 << s) - 0x100000000) | f;
        bit = 0x20000 << s;
        if ( f & bit ) {
            if ((f & (bit - 1)) || (f & (bit << 1))) f += bit;
        }
        f = f >> (s + 18);
    }
    if (sign) f = 0x10000 - f;    // restore sign
    uA.ui = f;
    return uA.p;

}


uint_fast64_t p16_logpoly( uint_fast64_t f ) {

    uint_fast64_t s, z, zsq;

    z = ((f << 31) + 2) / (f + 8192);    // fixed-point divide; discard remainder
    zsq = (( z * z ) >> 30);        // fixed-point squaring
    s = (zsq * 1584) >> 28;
    s = (zsq * (26661 + s)) >> 29;
    s = (zsq * (302676 + s)) >> 27;
    s = (zsq * (16136153 + s)) >> 30;
    s = (z   * (193635259 + s)) >> 27;

    return s;

}
