#ifndef SUB_PERM_H
#define SUB_PERM_H

int cmp(uint64_t a, uint64_t b);
double complex prod(uint64_t n, double complex *vec);
int sum0axis(uint64_t n, double complex **mat, double complex *vec);
double complex perm(uint64_t n, double complex **M);

#endif
