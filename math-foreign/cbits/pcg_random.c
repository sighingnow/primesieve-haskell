#if defined(__cplusplus)
extern "C" {
#endif  // __cplusplus

#include <assert.h>
#include <math.h>

#include "pcg_basic.h"
#include "pcg_random.h"

uint32_t pcg32_random_float() {
    return ldexpf((float)pcg32_random(), -32);
}

void pcg32_random_array(uint32_t nlen, uint32_t *r) {
    uint32_t i = 0;
    for (i = 0; i < nlen; ++i) {
        r[i] = pcg32_random();
    }
}

void pcg32_random_float_array(uint32_t nlen, float *r) {
    uint32_t i = 0;
    for (i = 0; i < nlen; ++i) {
        r[i] = ldexpf((float)pcg32_random(), -32);
    }
}

// Gaussian distribution generator.
//
// Ref: https://github.com/miloyip/normaldist-benchmark/blob/master/src/boxmuller.cpp
void boxmuller(uint32_t nlen, float* data) {
    assert(nlen % 2 == 0);

    uint32_t x = 1;
    static const float twopi = 2.0f * 3.14159265358979323846;

    union {
        uint32_t u;
        float f;
    } u;

    for (uint32_t i = 0; i < nlen; i += 2) {
        x = x * 1664525 + 1013904223;
        u.u = (x >> 9) | 0x3F800000;
        float u1 = 1.0f - (u.f - 1.0f);

        x = x * 1664525 + 1013904223;
        u.u = (x >> 9) | 0x3F800000;
        float u2 = u.f - 1.0f;

        float radius = sqrtf(-2 * logf(u1));
        float theta = twopi * u2;
        data[i] = radius * cosf(theta);
        data[i + 1] = radius * sinf(theta);
    }
}

#if defined(__cplusplus)
}
#endif  // __cplusplus
