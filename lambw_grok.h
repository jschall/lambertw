#ifndef LAMBERTW_H
#define LAMBERTW_H

#include <math.h>
#include <stdio.h>

/*
 * Computes the Lambert W function for the principal branch (k=0).
 * Input: x >= -1/e, where 1/e ≈ 0.36787944117
 * Output: W_0(x) such that W_0(x) * exp(W_0(x)) = x and W_0(x) >= -1
 */
float lambert_w0(float x) {
    const float INV_E = 0.36787944117f; /* 1/e with float precision */

    /* Check domain: x >= -1/e */
    if (x < -INV_E) {
        return NAN; /* Return NaN for invalid input */
    }

    /* Initial guess: x0 = x for all x >= -1/e */
    float w = x;

    /* Perform 17 iterations of Halley's method - determined experimentally*/
    for (int i = 0; i < 17; i++) {
        float ew = expf(w);         /* e^w */
        float f = w * ew - x;       /* f(w) = w * e^w - x */
        float fp = ew * (1.0f + w); /* f'(w) = e^w * (1 + w) */
        float fpp = ew * (2.0f + w);/* f''(w) = e^w * (2 + w) */
        float num = 2.0f * f * fp;  /* Numerator: 2 * f * f' */
        float den = 2.0f * fp * fp - f * fpp; /* Denominator: 2 * (f')^2 - f * f'' */

        /* Update w if denominator is non-zero */
        if (den != 0.0f) {
            w -= num / den;
        }
    }

    return w;
}

/*
 * Computes the Lambert W function for the k=-1 branch.
 * Input: y in [-1/e, 0), where 1/e ≈ 0.36787944117
 * Output: W_{-1}(y) such that W_{-1}(y) * exp(W_{-1}(y)) = y
 * Uses 32-bit floats, optimized for embedded systems with fixed iterations.
 */
float lambert_wm1(float y) {
    const float INV_E = 0.36787944117f; /* 1/e with float precision */

    /* Check domain: y must be in [-1/e, 0) */
    if (y < -INV_E || y >= 0.0f) {
        return NAN; /* Return NaN for invalid input */
    }

    /* Initial guess: x0 = ln(-y), which is real since -y > 0 */
    float x = logf(-y);

    /* Perform 7 iterations of Halley's method - determined experimentally */
    int i;
    for (i = 0; i < 7; i++) {
        float ex = expf(x);             /* Compute e^x once per iteration */
        if (ex == 0.0f) break;         /* Avoid underflow issues */
        float f = x * ex - y;           /* f(x) = x * e^x - y */
        float fp = ex * (1.0f + x);     /* f'(x) = e^x * (1 + x) */
        float fpp = ex * (2.0f + x);    /* f''(x) = e^x * (2 + x) */
        float num = 2.0f * f * fp;      /* Numerator for Halley's method */
        float den = 2.0f * fp * fp - f * fpp; /* Denominator */

        /* Update x if denominator is non-zero */
        if (den != 0.0f) {
            x -= num / den;
        }
    }

    return x;
}

#endif /* LAMBERTW_H */
