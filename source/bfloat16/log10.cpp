#include "bfloat16_math.hpp"

bfloat16 mylog10(bfloat16 x) {
    // If x == 0, then it should be -inf
    if (x.val == 0x0 || x.val == 0x8000) {
        x.val = 0xFF80;
        return x;
    }

    // If x == inf, then it should be infinity
    if (x.val == 0x7f80) {
        return x;
    }
    
    // If x == NaN or negative, then it should be NaN
    if (x.val > 0x7F80) {
        x.val = 0xFFFF;
        return x;
    }

    float fInput = (float)x;
    
    // Extract exponent and mantissa (where mantissa is between [1, 2)
    int m;
    float fx = frexpf(fInput, &m);
    fx *= 2.0;
    m--;
    
    // Cody and Waite Transformation on input
    double dx = (double)fx;
    double codyX = (dx - 1) / (dx + 1);
    double codyX2 = codyX * codyX;
    
    // Now compute polynomial
    double y = 6.710954935542725596775426311069168150424957275390625e-01;
    y *= codyX2;
    y += 9.56484867363945223672772044665180146694183349609375e-01;
    y *= codyX2;
    y += 2.88545942229525831379532974096946418285369873046875;
    y *= codyX;
    
    // Range propagation
    y = (m + y) / 3.321928094887362181708567732130177319049835205078125;
    return y;
}

double mylog10Internal(float x) {
    fx input;
    input.f = x;
    
    // If x == 0, then it should be -inf
    if (input.x == 0x0 || input.x == 0x80000000) {
        input.x = 0xFF800000;
        return input.f;
    }

    // If x == inf, then it should be infinity
    if (input.x == 0x7f800000) {
        return x;
    }
    
    // If x == NaN or negative, then it should be NaN
    if (input.x > 0x7F800000) {
        input.x = 0xFFFFFFFF;
        return input.f;
    }
    
    
    // Extract exponent and mantissa (where mantissa is between [1, 2)
    int m;
    float fx = frexpf((float)x, &m);
    fx *= 2.0;
    m--;
    
    // Cody and Waite Transformation on input
    double dx = (double)fx;
    double codyX = (dx - 1) / (dx + 1);
    double codyX2 = codyX * codyX;
    
    // Now compute polynomial
    double y = 6.710954935542725596775426311069168150424957275390625e-01;
    y *= codyX2;
    y += 9.56484867363945223672772044665180146694183349609375e-01;
    y *= codyX2;
    y += 2.88545942229525831379532974096946418285369873046875;
    y *= codyX;
    
    // Range propagation
    y = (m + y) / 3.321928094887362181708567732130177319049835205078125;
    return y;
}
