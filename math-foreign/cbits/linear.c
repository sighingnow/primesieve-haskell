#include <assert.h>
#include <complex.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "linear.h"

#include <cblas.h>
#include <lapacke.h>

#if defined(__cplusplus)
extern "C" {
#endif  // __cplusplus

#if defined(DEBUG)
#include <stdio.h>
#define debug(fmt, ...) fprintf(stderr, "ccall: " fmt "\n", ##__VA_ARGS__)
#else
#define debug(...)
#endif  // DEBUG

/* max and min */
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

/* zero test */
static inline __attribute__((always_inline)) bool s_zero(float x) { return x > -1e-6 && x < 1e-6; }

/* zero test */
static inline __attribute__((always_inline)) bool d_zero(double x) { return x > -1e-6 && x < 1e-6; }

/* zero test */
static inline __attribute__((always_inline)) bool c_zero(float complex x) {
    return crealf(x) > -1e-6 && crealf(x) < 1e-6 && cimagf(x) > -1e-6 && cimagf(x) < 1e-6;
}

/* zero test */
static inline __attribute__((always_inline)) bool z_zero(double complex x) {
    return creal(x) > -1e-6 && creal(x) < 1e-6 && cimag(x) > -1e-6 && cimag(x) < 1e-6;
}

#define MAKE_API(V0, V1, V2, V3, V4, V5)                                         \
    switch (type) {                                                              \
        case 1:                                                                  \
            MAKE_PROG(float, s##V0, s##V1, s##V2, s##V3, s##V4, s##V5);          \
            break;                                                               \
        case 2:                                                                  \
            MAKE_PROG(double, d##V0, d##V1, d##V2, d##V3, d##V4, d##V5);         \
        case 3:                                                                  \
            MAKE_PROG(float complex, c##V0, c##V1, c##V2, c##V3, c##V4, c##V5);  \
            break;                                                               \
        case 4:                                                                  \
            MAKE_PROG(double complex, z##V0, z##V1, z##V2, z##V3, z##V4, z##V5); \
        default:                                                                 \
            break;                                                               \
    }

#define MAKE_API_REAL(V0, V1, V2, V3, V4, V5)                            \
    switch (type) {                                                      \
        case 1:                                                          \
            MAKE_PROG(float, s##V0, s##V1, s##V2, s##V3, s##V4, s##V5);  \
            break;                                                       \
        case 2:                                                          \
            MAKE_PROG(double, d##V0, d##V1, d##V2, d##V3, d##V4, d##V5); \
        default:                                                         \
            break;                                                       \
    }

#define MAKE_API_COMPLEX(V0, V1, V2, V3, V4, V5)                                 \
    switch (type) {                                                              \
        case 3:                                                                  \
            MAKE_PROG(float complex, c##V0, c##V1, c##V2, c##V3, c##V4, c##V5);  \
            break;                                                               \
        case 4:                                                                  \
            MAKE_PROG(double complex, z##V0, z##V1, z##V2, z##V3, z##V4, z##V5); \
        default:                                                                 \
            break;                                                               \
    }

int identity(int type, void *r, int m, int n) {
    int i, loop = m > n ? n : m;
    debug("identity is called.");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5) \
    for (i = 0; i < loop; ++i) {             \
        ((T *)r)[i * n + i] = 1.;            \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int random_(int type, void *r, int m, int n) {
    int i, loop = m * n;
    srand(time(NULL));
    debug("random is called.");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5)    \
    for (i = 0; i < loop; ++i) {                \
        ((T *)r)[i] = (float)rand() / RAND_MAX; \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int diag(int type, void *r, int row, int column, const void *src) {
    int i;
    debug("diag is called.");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5)      \
    memset(r, 0x00, row *column * sizeof(T));     \
    for (i = 0; i < row; ++i) {                   \
        ((T *)r)[i * column + i] = ((T *)src)[i]; \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int diagonal(int type, void *r, const void *src, int row, int column) {
    int i, loop = row > column ? column : row;
    debug("diagonal is called.");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5)      \
    for (i = 0; i < loop; ++i) {                  \
        ((T *)r)[i] = ((T *)src)[i * column + i]; \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int sum(int type, void *r, const void *src, int row, int column) {
    int i;
    float srv = 0.;
    double drv = 0.;
    float complex crv = 0. + 0. * I;
    double complex zrv = 0. + 0. * I;
    debug("sum is called.");
#define MAKE_PROG(T, ACC, V1, V2, V3, V4, V5) \
    for (i = 0; i < row * column; ++i) {      \
        ACC += ((T *)src)[i];                 \
    }                                         \
    *((T *)r) = ACC;
    MAKE_API(rv, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int product(int type, void *r, const void *src, int row, int column) {
    int i;
    float srv = 1.;
    double drv = 1.;
    float complex crv = 1. + 0. * I;
    double complex zrv = 1. + 0. * I;
    debug("product is called.");
#define MAKE_PROG(T, ACC, V1, V2, V3, V4, V5) \
    for (i = 0; i < row * column; ++i) {      \
        ACC *= ((T *)src)[i];                 \
    }                                         \
    *((T *)r) = ACC;
    MAKE_API(rv, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int mean(int type, void *r, const void *src, int row, int column) {
    int i;
    float srv = 0.;
    double drv = 0.;
    float complex crv = 0. + 0. * I;
    double complex zrv = 0. + 0. * I;
    debug("mean is called.");
#define MAKE_PROG(T, ACC, V1, V2, V3, V4, V5) \
    for (i = 0; i < row * column; ++i) {      \
        ACC += ((T *)src)[i];                 \
    }                                         \
    *((T *)r) = ACC / (row * column);
    MAKE_API(rv, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int transpose(int type, void *r, const void *src, int row, int column) {
    int i, j;
    debug("transpose is called.");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5)                    \
    for (i = 0; i < row; ++i) {                                 \
        for (j = 0; j < column; ++j) {                          \
            ((T *)r)[j * row + i] = ((T *)src)[i * column + j]; \
        }                                                       \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int lower(int type, void *r, const void *src, int row, int column) {
    int i, j, loop = row > column ? column : row;
    debug("lower triangularize is called");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5) \
    memcpy(r, src, row *column * sizeof(T)); \
    for (i = 0; i < loop; ++i) {             \
        for (j = i; j < column; ++j) {       \
            ((T *)r)[i * column + j] = 0.;   \
        }                                    \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int upper(int type, void *r, const void *src, int row, int column) {
    int i, j, loop = row > column ? column : row;
    debug("upper triangularize is called");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5)                       \
    for (i = 0; i < loop; ++i) {                                   \
        for (j = i; j < column; ++j) {                             \
            ((T *)r)[i * column + j] = ((T *)src)[i * column + j]; \
        }                                                          \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int shift(int type, void *r, void *x, const void *src, int row, int column) {
    int i;
    debug("shift is called.");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5)   \
    for (i = 0; i < row * column; ++i) {       \
        ((T *)r)[i] = ((T *)src)[i] + *(T *)x; \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int times(int type, void *r, void *x, const void *src, int row, int column) {
    int i;
    debug("times is called.");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5)   \
    for (i = 0; i < row * column; ++i) {       \
        ((T *)r)[i] = ((T *)src)[i] * *(T *)x; \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int add(int type, void *r, int m, int n, int k, const void *A, const void *B) {
    int i, loop = m * n;
    assert(m * k == k * n);
    debug("add is called.");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5)     \
    for (i = 0; i < loop; ++i) {                 \
        ((T *)r)[i] = ((T *)A)[i] + ((T *)B)[i]; \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int minus(int type, void *r, int m, int n, int k, const void *A, const void *B) {
    int i, loop = m * n;
    assert(m * k == k * n);
    debug("minus is called.");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5)     \
    for (i = 0; i < loop; ++i) {                 \
        ((T *)r)[i] = ((T *)A)[i] - ((T *)B)[i]; \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int mult(int type, void *r, int m, int n, int k, const void *A, const void *B) {
    int i, loop = m * n;
    assert(m * k == k * n);
    debug("mult is called.");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5)     \
    for (i = 0; i < loop; ++i) {                 \
        ((T *)r)[i] = ((T *)A)[i] * ((T *)B)[i]; \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int division(int type, void *r, int m, int n, int k, const void *A, const void *B) {
    int i, loop = m * n;
    assert(m * k == k * n);
    debug("division is called.");
#define MAKE_PROG(T, V0, V1, V2, V3, V4, V5)     \
    for (i = 0; i < loop; ++i) {                 \
        ((T *)r)[i] = ((T *)A)[i] / ((T *)B)[i]; \
    }
    MAKE_API(NULL, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int dot(int type, void *r, int m, int n, int k, const void *A, const void *B) {
    debug("dot is called: (%d, %d) x (%d, %d) -> (%d, %d).", m, k, k, n, m, n);
    float salpha = 1., sbeta = 0.;
    double dalpha = 1., dbeta = 0.;
    float complex _calpha = 1. + 0. * I, _cbeta = 0. + 0. * I;
    double complex _zalpha = 1. + 0. * I, _zbeta = 0. + 0. * I;
    float complex *calpha = &_calpha, *cbeta = &_cbeta;
    double complex *zalpha = &_zalpha, *zbeta = &_zbeta;
#define MAKE_PROG(T, GEMM, ALPHA, BETA, V3, V4, V5)                                          \
    cblas_##GEMM(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, ALPHA, (const T *)A, k, \
                 (const T *)B, n, BETA, (T *)r, n);
    MAKE_API(gemm, alpha, beta, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int inner(int type, void *r, int n, const void *A, int stepa, int offa, const void *B, int stepb,
          int offb) {
    debug("inner is called: length of the vectors is %d.", n);
#define MAKE_PROG(T, DOT, V1, V2, V3, V4, V5) \
    *((T *)r) = cblas_##DOT(n, ((const T *)A) + offa, stepa, ((const T *)B) + offb, stepb);
    MAKE_API_REAL(dot, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
#define MAKE_PROG(T, DOTC_SUB, V1, V2, V3, V4, V5) \
    cblas_##DOTC_SUB(n, ((const T *)A) + offa, stepa, ((const T *)B) + offb, stepb, r);
    MAKE_API_COMPLEX(dotc_sub, NULL, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int det(int type, void *r, const void *src, int row, int column) {
    float srv = 1.;
    double drv = 1.;
    float complex crv = 1. + 0. * I;
    double complex zrv = 1. + 0. * I;
    int i, singular = 0;
    int *ipiv = (int *)malloc(column * sizeof(int));
    memset(ipiv, 0x00, column * sizeof(int));
    float *sp = NULL;
    double *dp = NULL;
    float complex *cp = NULL;
    double complex *zp = NULL;
    debug("det is called.");

// step 1: copy
//
// step 2: LU factorization by getrf, if singular > 0, then the input matrix is singular.
//
// step 3: multiply all pivots
//
// step 4: clear intermediate memory

#define MAKE_PROG(T, COPY, GETRF, MA, ACC, V4, V5)                            \
    MA = (T *)malloc(row * column * sizeof(T));                               \
    cblas_##COPY(row *column, (const T *)src, 1, MA, 1);                      \
    singular = LAPACKE_##GETRF(CblasRowMajor, row, column, MA, column, ipiv); \
    if (singular == 0) {                                                      \
        for (i = 0; i < column; ++i) {                                        \
            if (ipiv[i] != i) {                                               \
                ACC *= -1.0;                                                  \
            }                                                                 \
            ACC *= MA[i * column + i];                                        \
        }                                                                     \
    } else {                                                                  \
        ACC = 0.;                                                             \
    }                                                                         \
    *((T *)r) = ACC;                                                          \
    if (MA != NULL) {                                                         \
        free(MA);                                                             \
    }
    MAKE_API(copy, getrf, p, rv, NULL, NULL);
#undef MAKE_PROG

    free(ipiv);

    return 0;
}

int trace(int type, void *r, const void *src, int row, int column) {
    int i, loop = row > column ? column : row;
    float srv = 0.;
    double drv = 0.;
    float complex crv = 0. + 0. * I;
    double complex zrv = 0. + 0. * I;
    debug("trace is called.");
#define MAKE_PROG(T, F, ACC, V2, V3, V4, V5) \
    for (i = 0; i < loop; ++i) {             \
        ACC += ((T *)src)[i * column + i];   \
    }                                        \
    *((T *)r) = ACC;
    MAKE_API(NULL, rv, NULL, NULL, NULL, NULL);
#undef MAKE_PROG
    return 0;
}

int rank(int type, int *r, const void *src, int row, int column) {
    int i;
    int *ipiv = (int *)malloc(column * sizeof(int));
    memset(ipiv, 0x00, column * sizeof(int));

    float *sp = NULL;
    double *dp = NULL;
    float complex *cp = NULL;
    double complex *zp = NULL;

    float *stau = NULL;
    double *dtau = NULL;
    float complex *ctau = NULL;
    double complex *ztau = NULL;

    debug("rank is called, original matrix: %d x %d.", row, column);
    *r = 0;

// ref: http://people.sc.fsu.edu/~jburkardt/f_src/geqp3/geqp3.html
//
// step 1: copy
//
// step 2: LU factorization by getrf, if singular > 0, then the input matrix is singular.

#define MAKE_PROG(T, COPY, GEQP3, _ZERO, MA, TAU, V5)                   \
    MA = (T *)malloc(row * column * sizeof(T));                         \
    cblas_##COPY(row *column, (T *)src, 1, MA, 1);                      \
    TAU = (T *)malloc(min(row, column) * sizeof(T));                    \
    LAPACKE_##GEQP3(CblasRowMajor, row, column, MA, column, ipiv, TAU); \
    for (i = 0; i < row && i < column; ++i) {                           \
        if (!_ZERO(MA[i * column + i])) {                               \
            *r += 1;                                                    \
        }                                                               \
    }                                                                   \
    if (MA != NULL) {                                                   \
        free(MA);                                                       \
    }                                                                   \
    if (TAU != NULL) {                                                  \
        free(TAU);                                                      \
    }
    MAKE_API(copy, geqp3, _zero, p, tau, NULL)
#undef MAKE_PROG

    free(ipiv);

    return 0;
}

// TODO
int norm(int type, void *r, const void *src, int row, int column) {
    debug("norm is called.");
    return 0;
}

int inverse(int type, void *r, const void *src, int row, int column) {
    int singular = 0;
    int *ipiv = (int *)malloc(column * sizeof(int));
    debug("inverse is called.");

// step 1: do LU decomposition via getrf
//
// step 2: do inversion via getri

#define MAKE_PROG(T, COPY, GETRF, GETRI, V3, V4, V5)                              \
    cblas_##COPY(row *column, (T *)src, 1, (T *)r, 1);                            \
    singular = LAPACKE_##GETRF(CblasRowMajor, row, column, (T *)r, column, ipiv); \
    if (singular == 0) {                                                          \
        LAPACKE_##GETRI(CblasRowMajor, row, (T *)r, column, ipiv);                \
    } else {                                                                      \
        memset(r, 0xff, row *column * sizeof(T));                                 \
    }
    MAKE_API(copy, getrf, getri, NULL, NULL, NULL);
#undef MAKE_PROG

    free(ipiv);
    return 0;
}

// Eigen system of ordinary matrix, AX = \Lambda X.
int eigen(int type, void *rm1, int r1, int c1, void *rm2, int r2, int c2, void *rm3, int r3, int c3,
          const void *m) {
    return 0;
}

// AP = LU, where rm1 is L, rm2 is U, rm3 is P and m is A.
int lu(int type, void *rm1, int r1, int c1, void *rm2, int r2, int c2, void *rm3, int r3, int c3,
       const void *m) {
    int i, j, t;
    float *sp;
    double *dp;
    //     int *ipiv = (int *)malloc(c2 * sizeof(int));
    //     int *permutation = (int *)malloc(c2 * sizeof(int));
    //     assert(c1 == r2);
    //     assert(c2 == c3);
    //     assert(r3 == c3);
    //     memset(ipiv, 0xff, sizeof(int) * c2);
    //     memset(permutation, 0x00, sizeof(int) * c2);
    //     debug("lu decompose is called.");
    // #define MAKE_PROG(T, V, F) \
//     if (r1 > c1) {         \
//         V = (T *)rm1;      \
//     } else {               \
//         V = (T *)rm2;      \
//     }                      \
//     cblas_##F(r1 *c2, (T *)m, 1, V, 1);
    //     MAKE_API(p, copy);
    // #undef MAKE_PROG

    // #define MAKE_PROG(T, V, F)                         \
// /* do LU decomposition, compute matrix L and U. */ \
// // if (r1 > c1) {                                       \
//     //     clapack_##F(CblasRowMajor, r1, c1, V, c1, ipiv); \
//     //     for (i = 0; i < r1 && i < c2; ++i) {             \
//     //         j = i;                                       \
//     //         while (++j < c1) {                           \
//     //             ((T *)rm2)[i * c2 + j] = V[i * c1 + j];  \
//     //             V[i * c1 + j] = 0.0;                     \
//     //         }                                            \
//     //         ((T *)rm2)[i * c2 + i] = 1.0;                \
//     //     }                                                \
//     // } else {                                             \
//     //     clapack_##F(CblasRowMajor, r2, c2, V, c2, ipiv); \
//     //     for (i = 0; i < r1 && i < c2; ++i) {             \
//     //         for (j = 0; j <= i; ++j) {                   \
//     //             ((T *)rm1)[i * c1 + j] = V[i * c2 + j];  \
//     //             V[i * c2 + j] = 0.0;                     \
//     //         }                                            \
//     //         ((T *)rm2)[i * c2 + i] = 1.0;                \
//     //     }                                                \
//     // }
    // // MAKE_API(p, getrf);
    // #undef MAKE_PROG

    //     /* Compute permutation matrix P according to ipiv. */
    //     for (i = 0; i < c2; ++i) {
    //         permutation[i] = i;
    //     }
    //     for (i = 0; i < c2; ++i) {
    //         if (ipiv[i] >= 0 && ipiv[i] != i) {
    //             t = permutation[i];
    //             permutation[i] = permutation[ipiv[i]];
    //             permutation[ipiv[i]] = t;
    //         }
    //     }

    // #define MAKE_PROG(T, V, F)                         \
//     for (i = 0; i < r3; ++i) {                     \
//         ((T *)rm3)[permutation[i] * c3 + i] = 1.0; \
//     }
    //     MAKE_API(NULL, NULL);
    // #undef MAKE_PROG
    //     free(ipiv);
    //     free(permutation);
    return 0;
}

int qr(int type, void *rm1, int r1, int c1, void *rm2, int r2, int c2, void *rm3, int r3, int c3,
       const void *m) {
    // #define MAKE_PROG(T, V, F) LAPACK_##F(CblasRowMajor, r1, c2, (const T
    // *)m, c2, (T *)rm1);
    //     MAKE_API(NULL, geqrf)
    // #undef MAKE_PROG
    return 0;
}
int svd(int type, void *rm1, int r1, int c1, void *rm2, int r2, int c2, void *rm3, int r3, int c3,
        const void *m) {
    return 0;
}
int jordan(int type, void *rm1, int r1, int c1, void *rm2, int r2, int c2, void *rm3, int r3,
           int c3, const void *m) {
    return 0;
}
int cholesky(int type, void *rm1, int r1, int c1, void *rm2, int r2, int c2, void *rm3, int r3,
             int c3, const void *m) {
    return 0;
}
int schur(int type, void *rm1, int r1, int c1, void *rm2, int r2, int c2, void *rm3, int r3, int c3,
          const void *m) {
    return 0;
}

// r = Av
int transform(int type, void *r, int row, int column, const void *A, const void *v) {
    //     debug("transform is called: matrix has size (%d, %d).", row, column);
    // #define MAKE_PROG(T, V, F)                                                 \
//     cblas_##F(CblasRowMajor, CblasNoTrans, row, column, 1.0, (const T *)A, \
//               column, (const T *)v, 1, 0.0, (T *)r, 1);
    //     MAKE_API(NULL, gemv)
    // #undef MAKE_PROG
    //     return 0;
}

#if defined(__cplusplus)
}
#endif  // __cplusplus
