#include "posit16_math.h"

posit16_t cospi(posit16_t x) {
    if (x.v == 0x8000) {
        // Take care of NaN
        x.v = 0x8000;
        return x;
    }
    
    double dx = convertP16ToDouble(x);
    double modifier = 1;

    // cos(-x) = cos(x)
    if (dx < 0.0f) dx = -1 * dx;
    
    // How do we reduce range of x?
    // Reduce x to [0, 1)
    double intPart;
    double frac = modf(dx, &intPart);
    int iIntPart = intPart;
    
    // if iIntPart is odd, then flip modifier
    if (iIntPart % 2 == 1) modifier *= -1;
    
    // cos(pi - x) = -cos(x)
    if (frac >= 0.5) {
        frac = 1.0 - frac;
        modifier *= -1;
    }
    
    double y = 0;
    if (frac <= 0.003509521484375) {
        y = 1.0001220703125;
    } else if (frac < 0.5f) {
        double xSquared = frac * frac;
        y = 2.215338495769658688772096866159699857234954833984375e-01;
        y *= xSquared;
        y += -1.3327362938689424343152722940430976450443267822265625;
        y *= xSquared;
        y += 4.05853647916781223869975292473100125789642333984375;
        y *= xSquared;
        y += -4.93479863229652071510145106003619730472564697265625;
        y *= xSquared;
        y += 1.000000009410458634562246515997685492038726806640625;
    } else {
        y = 0.0f;
    }
    
    y *= modifier;
    return convertDoubleToP16(y);
}

double cospiInternal(double dx) {
    dint ch;
    ch.d = dx;
    if ((ch.x & 0x7FFFFFFFFFFFFFFF) > 0x7FF0000000000000) {
        return 0.0/0.0;
    }
    
    double modifier = 1;

    // cos(-x) = cos(x)
    if (dx < 0.0f) dx = -1 * dx;
    
    // How do we reduce range of x?
    // Reduce x to [0, 1)
    double intPart;
    double frac = modf(dx, &intPart);
    int iIntPart = intPart;
    
    // if iIntPart is odd, then flip modifier
    if (iIntPart % 2 == 1) modifier *= -1;
    
    // cos(pi - x) = -cos(x)
    if (frac >= 0.5) {
        frac = 1.0 - frac;
        modifier *= -1;
    }
    
    double y = 0;
    if (frac <= 0.003509521484375) {
        y = 1.0001220703125;
    } else if (frac < 0.5f) {
        double xSquared = frac * frac;
        y = 2.215338495769658688772096866159699857234954833984375e-01;
        y *= xSquared;
        y += -1.3327362938689424343152722940430976450443267822265625;
        y *= xSquared;
        y += 4.05853647916781223869975292473100125789642333984375;
        y *= xSquared;
        y += -4.93479863229652071510145106003619730472564697265625;
        y *= xSquared;
        y += 1.000000009410458634562246515997685492038726806640625;
    } else {
        y = 0.0f;
    }
    
    y *= modifier;
    return y;
}

uint_fast64_t p16_cospipoly(uint_fast64_t);

posit16_t p16_cospi(posit16_t pA) {

    union ui16_p16 uA;
    uint_fast16_t uiA;
    uint_fast64_t bit, f, sign = 0;
    int_fast64_t s;                     // s can be negative

    uA.p = pA;
    uiA = uA.ui;
    f = uiA;
    if (f == 0x8000)
        return uA.p;                // dispense with the NaR case

    if (f & 0x8000) f = 0x10000 - f;    // f = |f|

    if (f) {
        if (f & 0x4000) {           // decode regime
            s = 16;
            while (f & 0x2000) {
                f = f << 1;
                s += 2;
            }
        } else {
            s = 14;
            while (!(f & 0x2000)) {
                f = f << 1;
                s -= 2;
            }
        }
        if (f & 0x1000) s++;        // decode exponent
        f = (f & 0x0FFF) | 0x1000;  // get 12-bit fraction and restore hidden bit
        f = (s < 0) ? f >> -s : f << s;
    }
    s = f >> 27;                    // the quadrant is the multiple of 1/2
    f = f & 0x7FFFFFF;              // input value modulo 1/2
    if ((s + 1) & 2) sign = 0x8000; // cos is negative for quadrants 2 and 3
    if (!f) {
        uA.ui = (s & 1) ? 0 : (sign | 0x4000);
        return uA.p;
    }
    if (s & 1) f = 0x8000000 - f;
    f = p16_cospipoly(f);
    s = 1;                          // convert fixed-point to a posit
    while (!(f & 0x1000000)) {
        f = f << 1;
        s++;
    }
    bit = (s & 1);
    if (!bit) f = (f & 0xFFFFFF);   // encode exponent bit
    s = (s >> 1) + 12;
    if (!bit) s--;

    f = (f | 0x2000000);            // encode regime termination bit
    bit = ((uint_fast64_t)1 << (s - 1));
    if (f & bit) {                  // round to nearest, tie to even
        if ((f & (bit - 1)) || (f & (bit << 1))) f += bit;
    }
    f = (f >> s);
    uA.ui = (sign) ? (0x10000 - f) : f;
    return uA.p;

}

uint_fast64_t p16_cospipoly(uint_fast64_t f) {

    uint_fast64_t fsq, s;

    if (f < 0xE6001) return 0x1FFFFFF;  // this rounds up to 1.0
    fsq = (f >> 11);                    // convert to 17-bit fixed point
    fsq = ((fsq * fsq) >> 8);
    s = 349194 - ((fsq * 28875) >> 25);
    s = 4255560 - ((fsq * s) >> 24);
    s = 20698014 - ((fsq * s) >> 24);
    return 33554428 - ((fsq * s) >> 23);

}
