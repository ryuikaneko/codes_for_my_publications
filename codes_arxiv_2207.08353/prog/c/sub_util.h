#ifndef SUB_UTIL_H
#define SUB_UTIL_H

double gettimeofday_sec();
void *malloc2d(size_t size, int n1, int n2);
int print_mat(double complex **matA, int L);
int print_vec(double complex *vecA, int L);

#endif
