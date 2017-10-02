#ifndef __LINEAR_H__
#define __LINEAR_H__

#if defined(__cplusplus)
extern "C" {
#endif  // __cplusplus

int identity(int type, void *r, int m, int n);
int random_(int type, void *r, int m, int n);
int diag(int type, void *r, int row, int column, const void *src);
int diagonal(int type, void *r, const void *src, int row, int column);

int sum(int type, void *r, const void *src, int row, int column);
int product(int type, void *r, const void *src, int row, int column);
int mean(int type, void *r, const void *src, int row, int column);

int transpose(int type, void *r, const void *src, int row, int column);
int lower(int type, void *r, const void *src, int row, int column);
int upper(int type, void *r, const void *src, int row, int column);

int shift(int type, void *r, void *x, const void *src, int row, int column);
int times(int type, void *r, void *x, const void *src, int row, int column);
int add(int type, void *r, int m, int n, int k, const void *A, const void *B);
int minus(int type, void *r, int m, int n, int k, const void *A, const void *B);
int mult(int type, void *r, int m, int n, int k, const void *A, const void *B);
int division(int type, void *r, int m, int n, int k, const void *A, const void *B);
int dot(int type, void *r, int m, int n, int k, const void *A, const void *B);

int inner(int type, void *r, int n, const void *A, int stepa, int offa, const void *B, int stepb,
          int offb);

int det(int type, void *r, const void *src, int row, int column);
int trace(int type, void *r, const void *src, int row, int column);
int rank(int type, int *r, const void *src, int row, int column);
int norm(int type, void *r, const void *src, int row, int column);
int inverse(int type, void *r, const void *src, int row, int column);

int eigen(int type, void *A, int r0, int c0, void *Lambda, int r1, int c1, void *UT, int r2, int c2,
          void *V, int r3, int c3);
int eigenh(int type, void *A, int r0, int c0, void *Lambda, int r1, int c1, void *Z, int r2,
           int c2);

int lu(int type, void *rm1, int r1, int c1, void *rm2, int r2, int c2, void *rm3, int r3, int c3,
       const void *m);
int qr(int type, void *A, int r0, int c0, void *Q, int r1, int c1, void *R, int r2, int c2);
int svd(int type, void *A, int r0, int c0, void *U, int r1, int c1, void *Sigma, int r2, int c2,
        void *VT, int r3, int c3);
int jordan(int type, void *rm1, int r1, int c1, void *rm2, int r2, int c2, void *rm3, int r3,
           int c3, const void *m);
int cholesky(int type, char uplo, void *A, int r0, int c0, void *U, int r1, int c1);
int schur(int type, void *rm1, int r1, int c1, void *rm2, int r2, int c2, void *rm3, int r3, int c3,
          const void *m);

int transform(int type, void *r, int row, int column, const void *A, const void *v);

#if defined(__cplusplus)
}
#endif  // __cplusplus

#endif  // __LINEAR_H__
