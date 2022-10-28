#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <sys/time.h>
#include <math.h>
#include <complex.h>
#include "sub_util.h"
#include "sub_calmat.h"
#include "sub_perm.h"

int main(int argc, char *argv[]){
  int i;
  char opt;
  int L,L_A;
  double tstep;
  double complex val;
  double time_count[2];

//  L = 40;
  L = 24;
  tstep = 0.0;

//// prameter settings
  for(i = 0; i < argc; ++i){
    if(*argv[i] == '-'){
      opt = *(argv[i]+1);
      if(opt == 'L'){
        L = atoi(argv[i+1]);
      }else if(opt == 't'){
        tstep = atof(argv[i+1]);
      }
    }
  }

  L_A = L/2;

//// prepare matrix
//  double complex **matA = (double complex**)malloc2d(sizeof(double complex),2*L,2*L);// MI
  double complex **matA = (double complex**)malloc2d(sizeof(double complex),L,L);// CDW
  double complex **matX = (double complex**)malloc2d(sizeof(double complex),L,L);
  double complex **matY = (double complex**)malloc2d(sizeof(double complex),L,L);
//  double complex **matZ = (double complex**)malloc2d(sizeof(double complex),L,L);// MI
  double complex **matZ = (double complex**)malloc2d(sizeof(double complex),L/2,L/2);// CDW
  double complex *vecEPS = (double complex*)malloc(sizeof(double complex)*L);

//// open files
  FILE *fp_setup;
  FILE *fp_time;
  FILE *fp_renyi;
  char buf_setup[256];
  char buf_time[256];
  char buf_renyi[256];
  sprintf(buf_setup,"dat_setup_L%d_t%.10f",L,tstep);
  sprintf(buf_time,"dat_time_L%d_t%.10f",L,tstep);
  sprintf(buf_renyi,"dat_renyi_L%d_t%.10f",L,tstep);
  fp_setup = fopen(buf_setup,"w");
  fp_time = fopen(buf_time,"w");
  fp_renyi = fopen(buf_renyi,"w");

//// output parameters
  fprintf(fp_setup,"# L,tstep\n");
  fprintf(fp_setup,"%d %.10f\n",L,tstep);

  time_count[0] = gettimeofday_sec();

//// calculate entanglement
  calc_matA(matA,matX,matY,matZ,vecEPS,L,L_A,tstep);
//  val = perm(2*L,matA);// MI
  val = perm(L,matA);// CDW
  printf("%g %g %g %g\n",tstep,-log(fabs(creal(val))),creal(val),cimag(val));

//// output data
//// print tstep, -ln(perm(A)), perm(A).re, perm(A).im
  fprintf(fp_renyi,"%.10g %.10g %.10g %.10g\n",
    tstep,-log(fabs(creal(val))),creal(val),cimag(val)
    );

  time_count[1] = gettimeofday_sec();

  fprintf(fp_time,"%13.6f\n",time_count[1]-time_count[0]);

//// close files
  fclose(fp_setup);
  fclose(fp_time);
  fclose(fp_renyi);

//// free matrix
  free(matA);
  free(matX);
  free(matY);
  free(matZ);
  free(vecEPS);

  return 0;
}
