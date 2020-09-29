#include "bfloat16_math.hpp"

bfloat16 myexpv2(bfloat16 x) {
    float fInput = (float)x;
    
    // Take care of 0 case
    if (fInput <= -93.0f) {
        return 0.0;
    }
    
    // Take care of infinity case
    if (fInput >= 89.0f) {
        x.val = 0x7F80;
        return x;
    }
    
    // Take care of NaN
    if ((x.val & 0x7FFF) > 0x7F80) {
        return x;
    }
    
    if (fInput <= 3.890991e-03 && fInput >= -1.953125e-03) {
        return 1.0;
    }
    
    double xprime = fInput * 1.442695040888963387004650940070860087871551513671875;
    double intPart;
    double fracPart = modf(xprime, &intPart);
    if (fracPart < 0.0f) {
        fracPart += 1.0;
        intPart -= 1.0;
    }

    double y = 1.45139853027161404297462610202273936010897159576416015625e-02;
    y *= fracPart;
    y += 5.014719237694532927296364732683287002146244049072265625e-02;
    y *= fracPart;
    y += 2.42560224581628236517616414857911877334117889404296875e-01;
    y *= fracPart;
    y += 6.9279247181322956006255253669223748147487640380859375e-01;
    y *= fracPart;
    y += 1.0000095976211798021182630691328085958957672119140625;
    
    return ldexp(y, intPart);

}

double myexpInternalv2(float x) {
    fx input;
    input.f = x;
    
    // Take care of 0 case
    if (input.f <= -134.0f) {
        return 0.0;
    }
    
    // Take care of infinity case
    if (input.f >= 128.0f) {
        input.x = 0x7F800000;
        return input.f;
    }
    
    // Take care of NaN
    if ((input.x & 0x7FFFFFFF) > 0x7F800000) {
        return x;
    }
    
    if (input.f <= 3.890991e-03 && input.f >= -1.953125e-03) {
        return 1.0;
    }
    
    double xprime = (double)x * 1.442695040888963387004650940070860087871551513671875;
    double intPart;
    double fracPart = modf(xprime, &intPart);
    if (fracPart < 0.0f) {
        fracPart += 1.0;
        intPart -= 1.0;
    }

    double y = 1.45139853027161404297462610202273936010897159576416015625e-02;
    y *= fracPart;
    y += 5.014719237694532927296364732683287002146244049072265625e-02;
    y *= fracPart;
    y += 2.42560224581628236517616414857911877334117889404296875e-01;
    y *= fracPart;
    y += 6.9279247181322956006255253669223748147487640380859375e-01;
    y *= fracPart;
    y += 1.0000095976211798021182630691328085958957672119140625;
    
    return ldexp(y, intPart);
}
