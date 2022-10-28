#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <sys/time.h>
#include <math.h>
#include <complex.h>
#include "sub_util.h"

// calculate matrix
int calc_matA(double complex **matA, double complex **matX,
  double complex **matY, double complex **matZ, double complex *vecEPS,
  int L, int L_A, double tstep){
  int i,j,k;
  for(i=0; i<L; i++){
    for(j=0; j<L; j++){
      matX[i][j] = sin((i+1.0)*(j+1.0)*M_PI/(L+1.0)) * sqrt(2.0/(L+1.0));
    }
  }
  for(i=0; i<L; i++){
    vecEPS[i] = -2.0*cos((i+1.0)*M_PI/(L+1.0));
  }
  for(i=0; i<L; i++){
    for(j=0; j<L; j++){
      matY[i][j] = 0.0;
    }
  }
  for(i=0; i<L; i++){
    for(j=0; j<L; j++){
      for(k=0; k<L; k++){
        matY[i][j] += cexp(-I*vecEPS[k]*tstep) * conj(matX[k][i]) * matX[k][j];
      }
    }
  }
//// MI
//  for(i=0; i<L; i++){
//    for(j=0; j<L; j++){
//// CDW
  for(i=0; i<L/2; i++){
    for(j=0; j<L/2; j++){
      matZ[i][j] = 0.0;
    }
  }
//// MI
//  for(i=0; i<L; i++){
//    for(j=0; j<L; j++){
//// CDW
  for(i=0; i<L/2; i++){
    for(j=0; j<L/2; j++){
      for(k=0; k<L_A; k++){
//        matZ[i][j] += conj(matY[i][k]) * matY[j][k];// MI
        matZ[i][j] += conj(matY[2*i][k]) * matY[2*j][k];// CDW
      }
    }
  }
//// MI
//  for(i=0; i<L; i++){
//    for(j=0; j<L; j++){
//// CDW
  for(i=0; i<L/2; i++){
    for(j=0; j<L/2; j++){
//      matA[i][j] = matZ[i][j];
//      matA[i][j+L] = -matZ[i][j];
//      matA[i+L][j] = -matZ[i][j];
//      matA[i+L][j+L] = matZ[i][j];
      matA[i][j] = matZ[i][j];
      matA[i][j+L/2] = -matZ[i][j];
      matA[i+L/2][j] = -matZ[i][j];
      matA[i+L/2][j+L/2] = matZ[i][j];
    }
//    matA[i][i+L] += 1.0;
//    matA[i+L][i] += 1.0;
    matA[i][i+L/2] += 1.0;
    matA[i+L/2][i] += 1.0;
  }
  return 0;
}
