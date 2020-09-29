#include "bfloat16_math.hpp"

bfloat16 mycbrt(bfloat16 x) {
    
    // If x == 0, then it should be 0
    if (x.val == 0x0 || x.val == 0x8000) {
        return x;
    }

    // If x == inf, then it should be infinity
    if ((x.val & 0x7f80) == 0x7f80) {
        return x;
    }
    
    // If x == NaN, then it should be NaN
    if ((x.val & 0x7f80) > 0x7F80) {
        x.val = 0xFFFF;
        return x;
    }

    // Extract exponent and mantissa (where mantissa is between [1, 8))
    int m;
    float rx = frexpf((float)x, &m);
    rx *= 2.0;
    m--;
    
    int leftOver = m % 3;
    if (leftOver != 0) {
        if (leftOver < 0) leftOver += 3;
        
        fx fi;
        fi.x = 0x3f800000 + (leftOver << 23);
        rx *= fi.f;
        m -= leftOver;
    }

    m /= 3;
    if (rx < 0) rx *= -1;
    
    // Now compute polynomial
    double y = -1.7372029717703960593165601888898663673899136483669281005859375e-05;
    y *= rx;
    y += 5.241080546145838146843143334763226448558270931243896484375e-04;
    y *= rx;
    y += -6.5208421736825845915763721905022975988686084747314453125e-03;
    y *= rx;
    y += 4.3868412288261666998057108912689727731049060821533203125e-02;
    y *= rx;
    y += -1.80364291120356845521399691278929822146892547607421875e-01;
    y *= rx;
    y += 5.752913905623990853399618572439067065715789794921875e-01;
    y *= rx;
    y += 5.6860957346246798760347473944420926272869110107421875e-01;
    
    if (x < 0) y *= -1;
    return ldexp(y, m);
}

double mycbrtInternal(float x) {
    fx input;
    input.f = x;
    
    // If x == 0, then it should be -inf
    if (input.x == 0x0 || input.x == 0x80000000) {
        return x;
    }

    // If x == inf, then it should be infinity
    if ((input.x & 0x7f800000) == 0x7f800000) {
        return x;
    }
    
    // If x == NaN or negative, then it should be NaN
    if ((input.x & 0x7f800000) > 0x7F800000) {
        input.x = 0xFFFFFFFF;
        return input.f;
    }
    
    
    // Extract exponent and mantissa (where mantissa is between [1, 8))
    int m;
    float rx = frexpf((float)x, &m);
    rx *= 2.0;
    m--;
    
    int leftOver = m % 3;
    if (leftOver != 0) {
        if (leftOver < 0) leftOver += 3;
        
        fx fi;
        fi.x = 0x3f800000 + (leftOver << 23);
        rx *= fi.f;
        m -= leftOver;
    }

    m /= 3;
    if (rx < 0) rx *= -1;
    
    // Now compute polynomial
    double y = -1.7372029717724262278845392959425453227595426142215728759765625e-05;
    y *= rx;
    y += 5.24108054615050021618483100382945849560201168060302734375e-04;
    y *= rx;
    y += -6.520842173686691549405747281298317830078303813934326171875e-03;
    y *= rx;
    y += 4.3868412288279250155209609829398686997592449188232421875e-02;
    y *= rx;
    y += -1.80364291120395259238051721695228479802608489990234375e-01;
    y *= rx;
    y += 5.7529139056243938643575575042632408440113067626953125e-01;
    y *= rx;
    y += 5.6860957346245211141422259970568120479583740234375e-01;
    
    if (x < 0) y *= -1;
    return ldexp(y, m);
}
