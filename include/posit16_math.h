#include "math.h"
#include "softposit.h"
#include <cstdint>

typedef union {
    double d;
    unsigned long long x;
} dint;

posit16_t cospi(posit16_t x);
posit16_t log(posit16_t x);
posit16_t log2(posit16_t x);
posit16_t log10(posit16_t x);
posit16_t sinpi(posit16_t x);
posit16_t mysqrt(posit16_t x);
posit16_t exp(posit16_t x);
posit16_t exp2(posit16_t x);
posit16_t exp10(posit16_t x);

double cospiInternal(double dx);
double logInternal(double dx);
double log2Internal(double dx);
double log10Internal(double dx);
double sinpiInternal(double dx);
double mysqrtInternal(double dx);
double expInternal(double dx);
double exp2Internal(double dx);
double exp10Internal(double dx);

uint_fast64_t p16_cospipoly(uint_fast64_t);
uint_fast64_t p16_logpoly( uint_fast64_t );
uint_fast64_t p16_log2poly( uint_fast64_t );
uint_fast64_t p16_sinpipoly( uint_fast64_t );
uint_fast64_t p16_exppoly( uint_fast64_t );
uint_fast64_t p16_exp2poly( uint_fast64_t );

posit16_t p16_cospi(posit16_t pA);
posit16_t p16_log( posit16_t pA );
posit16_t p16_log2( posit16_t pA );
posit16_t p16_sinpi( posit16_t pA );
posit16_t p16_exp( posit16_t pA );
posit16_t p16_exp2( posit16_t pA );
