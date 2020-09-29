#include "posit16_math.h"

posit16_t sinpi(posit16_t x) {
    if (x.v == 0x8000) {
        // Take care of NaN
        x.v = 0x8000;
        return x;
    }
    
    double dx = convertP16ToDouble(x);
    double modifier = 1;

    // sin(-x) = -sin(x)
    if (dx < 0.0f) {
        dx = -1.0 * dx;
        modifier *= -1;
    }
    
    // How do we reduce range of x?
    // Reduce x to [0, 1)
    double intPart;
    double frac = modf(dx, &intPart);
    int iIntPart = intPart;
    
    // if iIntPart is odd, then flip modifier
    if (iIntPart % 2 == 1) modifier *= -1;
    
    // sin(x) = sin(pi - x)
    if (frac >= 0.5) {
        frac = 1.0 - frac;
    }
    double fracPart = frac;
    double y = 0;
    if (fracPart <= 0.00252532958984375000) {
        y = 3.141577060931899811890843920991756021976470947265625 * fracPart;
    } else {
        double xSquared = fracPart * fracPart;
        y = 9.47599641221426869375221713198698125779628753662109375e-02;
        y *= xSquared;
        y += -6.0547119473342603246379667325527407228946685791015625e-01;
        y *= xSquared;
        y += 2.55098424541712009983029929571785032749176025390625;
        y *= xSquared;
        y += -5.1677486367595673044661452877335250377655029296875;
        y *= xSquared;
        y += 3.141593069399674309494230328709818422794342041015625;
        y *= fracPart;
    }
    
    y *= modifier;
    return convertDoubleToP16(y);
}

double sinpiInternal(double dx) {
    dint ch;
    ch.d = dx;
    if ((ch.x & 0x7FFFFFFFFFFFFFFF) > 0x7FF0000000000000) {
        return 0.0/0.0;
    }
    
    double modifier = 1;

    // sin(-x) = -sin(x)
    if (dx < 0.0f) {
        dx = -1.0 * dx;
        modifier *= -1;
    }
    
    // How do we reduce range of x?
    // Reduce x to [0, 1)
    double intPart;
    double frac = modf(dx, &intPart);
    int iIntPart = intPart;
    
    // if iIntPart is odd, then flip modifier
    if (iIntPart % 2 == 1) modifier *= -1;
    
    // sin(x) = sin(pi - x)
    if (frac >= 0.5) {
        frac = 1.0 - frac;
    }
    double fracPart = frac;
    double y = 0;
    if (fracPart <= 0.00252532958984375000) {
        y = 3.141577060931899811890843920991756021976470947265625 * fracPart;
    } else {
        double xSquared = fracPart * fracPart;
        y = 9.47599641221426869375221713198698125779628753662109375e-02;
        y *= xSquared;
        y += -6.0547119473342603246379667325527407228946685791015625e-01;
        y *= xSquared;
        y += 2.55098424541712009983029929571785032749176025390625;
        y *= xSquared;
        y += -5.1677486367595673044661452877335250377655029296875;
        y *= xSquared;
        y += 3.141593069399674309494230328709818422794342041015625;
        y *= fracPart;
    }
    
    y *= modifier;
    return y;
}

uint_fast64_t p16_sinpipoly( uint_fast64_t );

posit16_t p16_sinpi( posit16_t pA ) {

    union ui16_p16 uA;
    uint_fast16_t uiA;
    uint_fast64_t bit, f, sign;
    int_fast64_t s;                // s can be negative

    uA.p = pA;
    uiA = uA.ui;
    f = uiA;
    sign = f & 0x8000;
    if (sign) f = 0x10000 - f;        // 2's complement if negative
    if (f > 31743) {            // input value is an integer?
        if (f == 0x8000) {
            uA.ui = 0x8000;        // sinpi(NaR) is NaR
            return uA.p;
        } else {
            uA.ui = 0;        // sinpi of an integer is zero
            return uA.p;
        }
    }
    if (!f) {                // sinpi(0) = 0
        uA.ui = 0;
        return uA.p;
    }
    if ( f & 0x4000 ) {            // decode regime
        s = 16;
        while ( f & 0x2000 ) {
            f = f << 1;
            s += 2;
        }
    } else {
        s = 14;
        while  ( !(f & 0x2000) ) {
            f = f << 1;
            s -= 2;
        }
    }
    if ( f & 0x1000 ) s++;            // decode exponent
    f = (f & 0x0FFF) | 0x1000;        // get 12-bit fraction and restore hidden bit
    f = (s < 0) ? f >> -s : f << s;
    f = (f & 0x1FFFFFFF);            // fixed-point with 28-bit fraction
    s = f >> 27;                // the quadrant is the multiple of 1/2
    f = f & 0x7FFFFFF;            // input value modulo 1/2
    if (s & 2) sign = sign ^ 0x8000;    // quadrants 2 and 3 flip the sign
    if (!f) {
        uA.ui = (s & 1) ? (sign | 0x4000) : 0;
        return uA.p;
    }
    if (s & 1) f = 0x8000000 - f;
    f = p16_sinpipoly(f);
    s = 1;                    // convert 28-bit fixed-point to a posit
    while (!(f & 0x8000000)) {
        f = f << 1;
        s++;
    }
    bit = (s & 1);
    s = (s >> 1) + 14 + bit;
    if (!bit) f = (f & 0x7FFFFFF);        // encode exponent bit
    f = (f | 0x10000000);            // encode regime termination bit
    bit = ((uint_fast64_t)1 << (s - 1));
    if (f & bit) {                // round to nearest, tie to even
        if ( (f & (bit - 1)) || (f & (bit << 1)) ) f += bit;
    }
    f = (f >> s);
    uA.ui = (sign) ? (0x10000 - f) : f;
    return uA.p;

}


uint_fast64_t p16_sinpipoly( uint_fast64_t f) {

    uint_fast64_t fs, fsq, s;

    if (f < 0xA5801) return (f * 102943) >> 15;    // linear approximation suffices
    fs = (f >> 11);
    fsq = ((fs * fs) >> 8);
    s = ((fsq * 650) >> 25);
    s = ((fsq * (9813 - s)) >> 23);
    s = ((fsq * (334253 - s)) >> 23);
    s = ((fsq * (5418741 - s)) >> 22);
    return (fs * (52707180 - s)) >> 13;
}
