#ifndef SUB_CALMAT_H
#define SUB_CALMAT_H

int calc_matA(double complex **matA, double complex **matX,
  double complex **matY, double complex **matZ, double complex *vecEPS,
  int L, int L_A, double tstep);

#endif
