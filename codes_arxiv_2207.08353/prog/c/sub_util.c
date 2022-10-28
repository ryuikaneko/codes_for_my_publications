#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <sys/time.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

// get time
double gettimeofday_sec()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

// memory allocation
void *malloc2d(size_t size, uint64_t n1, uint64_t n2){
  uint64_t i;
  uint64_t t=size*n2;
  char **a1, *a2;
  a1 = (char**)malloc((sizeof(*a1) + t) * n1);
  if(a1){
    a2 = (char*)(a1 + n1);
    for(i=0; i<n1; i++){
      a1[i] = a2;
      a2 += t;
    }
    return a1;
  }
  return NULL;
}

// print matrix
int print_mat(double complex **matA, uint64_t L){
  uint64_t i,j;
  for(i=0; i<L; i++){
    for(j=0; j<L; j++){
      printf("%+.6f%+.6fI ",creal(matA[i][j]),cimag(matA[i][j]));
    }
    printf("\n");
  }
  printf("\n");
  return 0;
}

// print vector
int print_vec(double complex *vecA, uint64_t L){
  uint64_t i;
  for(i=0; i<L; i++){
    printf("%+.6f%+.6fI ",creal(vecA[i]),cimag(vecA[i]));
  }
  printf("\n\n");
  return 0;
}
