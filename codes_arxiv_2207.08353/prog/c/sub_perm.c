#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <sys/time.h>
#include <math.h>
#include <complex.h>
#include "sub_util.h"

// https://stackoverflow.com/questions/10996418/efficient-integer-compare-function
//  -1 if a < b
//   0 if a = b
//  +1 if a > b
int cmp(uint64_t a, uint64_t b){
  return (a>b)-(a<b);
}

// calculate product
double complex prod(uint64_t n, double complex *vec){
  uint64_t i;
  double complex tmp;
  tmp = vec[0];
  for(i=1; i<n; i++){
    tmp *= vec[i];
  }
  return tmp;
}

// calculate sum along 0 axis
int sum0axis(uint64_t n, double complex **mat, double complex *vec){
  uint64_t i;
  uint64_t j;
  for(i=0; i<n; i++){
    vec[i] = 0.0;
  }
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      vec[j] += mat[i][j];
    }
  }
  return 0;
}

// permnanent
double complex perm(uint64_t n, double complex **M){
  uint64_t i;
  uint64_t old_gray;
  uint64_t new_gray;
  uint64_t gray_diff;
  uint64_t gray_diff_index;
  uint64_t num_loops;
  uint64_t bin_index;
  int sign;// takes -1,1
  int direction;// takes -2,0,+2
  double complex total;
  double complex *row_comb = (double complex*)malloc(n*sizeof(double complex));
  double complex reduced;

  if(n==0) return 1.0;
  sum0axis(n,M,row_comb);
  total = 0.0;
  old_gray = 0;
  sign = 1;
  num_loops = 1LLU << (n-1);
  for(bin_index=1; bin_index<=num_loops; bin_index++){
    reduced = prod(n,row_comb);
    total += sign * reduced;
    new_gray = bin_index ^ (bin_index / 2);
    gray_diff = old_gray ^ new_gray;
    gray_diff_index = __builtin_ctzll(gray_diff);
    direction = 2 * cmp(old_gray,new_gray);
    for(i=0; i<n; i++){
      row_comb[i] += M[gray_diff_index][i] * direction;
    }
    sign = -sign;
    old_gray = new_gray;
  }
  free(row_comb);
  return total / num_loops;
}
