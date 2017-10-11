#ifndef __PCG_RANDOM_H__
#define __PCG_RANDOM_H__

#if defined(__cplusplus)
extern "C" {
#endif  // __cplusplus

// extensions for pcg random generator.
uint32_t pcg32_random_float();
void pcg32_random_array(uint32_t nlen, uint32_t *r);
void pcg32_random_float_array(uint32_t nlen, float *r);

// Gaussian distribution generator.
void boxmuller(uint32_t nlen, float* data);

#if defined(__cplusplus)
}
#endif  // __cplusplus

#endif  // __PCG_RANDOM_H__
