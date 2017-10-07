#if defined(__cplusplus)
extern "C" {
#endif  // __cplusplus

#include <math.h>

#include "pcg_basic.h"
#include "pcg_random.h"

uint32_t pcg32_random_float() {
    return ldexp(pcg32_random(), -32);
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
        r[i] = ldexp(pcg32_random(), -32);
    }
}

#if defined(__cplusplus)
}
#endif  // __cplusplus
