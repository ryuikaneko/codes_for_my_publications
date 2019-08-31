/*-------------------------------------------------------------------*
 * exchange monte carlo
 *-------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "../dSFMT/dSFMT.h"

//#define M_PI 3.14159265358979323846264338327

#define min(a,b) (((a)<(b)) ? (a):(b))
dsfmt_t dsfmt;

/*-------------------------------------------------------------------*
 * get time
 *-------------------------------------------------------------------*/
double gettimeofday_sec()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

/*-------------------------------------------------------------------*
 * energy
 *-------------------------------------------------------------------*/
int calc_energy
(
  int size_L,
  int m,
  double parJ1,
  double parD,
  int ***spins,
  double *ene
)
{
  int i,ip,im,ipp,j,jp,jpp;
  double tmp_ene;
  tmp_ene = 0.0;
  for(i=0; i<size_L; i++){
    ip  = (i+1)%size_L;
    im  = (i-1+size_L)%size_L;
    ipp = (i+2)%size_L;
    for(j=0; j<size_L; j++){
      jp  = (j+1)%size_L;
      jpp = (j+2)%size_L;
      tmp_ene +=
        - parJ1 * spins[m][i][j]
          * (spins[m][ip][j] + spins[m][i][jp] + spins[m][ip][jp])
        - parD * spins[m][i][j] * spins[m][i][j];
    }
  }
  ene[m] = tmp_ene;
  return 0;
}

/*-------------------------------------------------------------------*
 * magnetization
 *-------------------------------------------------------------------*/
double calc_mag
(
  int size_L,
  int m,
  int ***spins,
  double **mag
)
{
  int i,j;
  double tmp_mag[4];
  tmp_mag[0] = 0.0;
  tmp_mag[1] = 0.0;
  tmp_mag[2] = 0.0;
  tmp_mag[3] = 0.0;
  for(i=0; i<size_L; i++){
    for(j=0; j<size_L; j++){
      tmp_mag[0] += spins[m][i][j];
      tmp_mag[1+(i+j)%3] += spins[m][i][j];
    }
  }
  mag[m][0] = tmp_mag[0];
  mag[m][1] = tmp_mag[1];
  mag[m][2] = tmp_mag[2];
  mag[m][3] = tmp_mag[3];
  return 0;
}

/*-------------------------------------------------------------------*
 * one monte carlo step
 *-------------------------------------------------------------------*/
int one_flip
(
  int size_L,
  int m,
  double ***BoltzFac,
  int *list_temp,
  int ***spins
)
{
  int egap1,egap2;
  int i,j,ip,im,ipp,imm,jp,jm,jpp,jmm;
  int sold,snew;
  /* spin flip candidate */
  i = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  j = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  /* calc egap */
  ip  = (i+1)%size_L;
  im  = (i-1+size_L)%size_L;
  ipp = (i+2)%size_L;
  imm = (i-2+size_L)%size_L;
  jp  = (j+1)%size_L;
  jm  = (j-1+size_L)%size_L;
  jpp = (j+2)%size_L;
  jmm = (j-2+size_L)%size_L;
  sold = spins[m][i][j];
  snew = ((int)(sold + 2*dsfmt_genrand_close_open(&dsfmt) + 2))%3 - 1;
  egap1 = (snew - sold)
    *(spins[m][ip][j] + spins[m][im][j]
    + spins[m][i][jp] + spins[m][i][jm]
    + spins[m][ip][jp] + spins[m][im][jm]);
  egap2 = snew*snew - sold*sold;
  /* metropolis method */
  if(dsfmt_genrand_close_open(&dsfmt) <= BoltzFac[list_temp[m]][egap1+12][egap2+1]){
    spins[m][i][j] = snew;
  }
  return 0;
}

/*-------------------------------------------------------------------*
 * one monte carlo step
 * (calculate energy and magnetization differences)
 *-------------------------------------------------------------------*/
int one_flip_ene_mag
(
  int size_L,
  int m,
  double parJ1,
  double parD,
  double ***BoltzFac,
  int *list_temp,
  int ***spins,
  double *ene,
  double **mag
)
{
  int egap1,egap2;
  int i,j,ip,im,ipp,imm,jp,jm,jpp,jmm;
  int sold,snew;
  /* spin flip candidate */
  i = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  j = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  /* calc egap */
  ip  = (i+1)%size_L;
  im  = (i-1+size_L)%size_L;
  ipp = (i+2)%size_L;
  imm = (i-2+size_L)%size_L;
  jp  = (j+1)%size_L;
  jm  = (j-1+size_L)%size_L;
  jpp = (j+2)%size_L;
  jmm = (j-2+size_L)%size_L;
  sold = spins[m][i][j];
  snew = ((int)(sold + 2*dsfmt_genrand_close_open(&dsfmt) + 2))%3 - 1;
  egap1 = (snew - sold)
    *(spins[m][ip][j] + spins[m][im][j]
    + spins[m][i][jp] + spins[m][i][jm]
    + spins[m][ip][jp] + spins[m][im][jm]);
  egap2 = snew*snew - sold*sold;
  /* metropolis method */
  if(dsfmt_genrand_close_open(&dsfmt) <= BoltzFac[list_temp[m]][egap1+12][egap2+1]){
    spins[m][i][j] = snew;
    ene[m] += (-parJ1) * egap1 + (-parD) * egap2;
    mag[m][0] += (snew-sold);
    mag[m][1+(i+j)%3] += (snew-sold);
  }
  return 0;
}

/*-------------------------------------------------------------------*
 * temperature exchange
 *-------------------------------------------------------------------*/
int exchange_states
(
  int size_L,
  int m,
  double parJ1,
  double parD,
  double *T,
  int *list_temp,
  int *list_conf,
  int ***spins,
  double *ene
)
{
  int conf0,conf1;
  double egap,w;
  /* calculate gap */
  conf0 = list_conf[m];
  conf1 = list_conf[m+1];
  egap = (1.0/T[m+1] - 1.0/T[m]) * (ene[conf0] - ene[conf1]);
  /* metropolis method */
  w = exp(-egap);
  if(dsfmt_genrand_close_open(&dsfmt) <= w){
    list_temp[conf0] = m+1;
    list_temp[conf1] = m;
    list_conf[m]     = conf1;
    list_conf[m+1]   = conf0;
  }
  return 0;
}

/*-------------------------------------------------------------------*
 * calculate physical quantity during sampling
 *-------------------------------------------------------------------*/
int calc_phys_quant
(
  int m,
  int index,
  int *list_temp,
  double **PhysQuant,
  double PhysAdd
)
{
  double d_temp;
  d_temp = PhysAdd;
  PhysQuant[list_temp[m]][index+0] += d_temp;
  d_temp = d_temp*d_temp;
  PhysQuant[list_temp[m]][index+1] += d_temp;
  d_temp = d_temp*d_temp;
  PhysQuant[list_temp[m]][index+2] += d_temp;
  return 0;
}

/*-------------------------------------------------------------------*
 * cos(6*atan2(y,x))
 *-------------------------------------------------------------------*/
double cos6atan
(
  double y,
  double x
)
{
  int sign;
  double t2;
  double d_temp;
  if(fabs(x)+fabs(y)<1.0e-8){
    d_temp = 0.0;
  }else{
    if(fabs(y)>fabs(x)){
      // t = x/y
      // ret = +(t^6-15t^4+15t^2-1)/((t^2+1)^3)
      sign = +1;
      t2 = x/y*x/y;
    }else{
      // t = y/x
      // ret = -(t^6-15t^4+15t^2-1)/((t^2+1)^3)
      sign = -1;
      t2 = y/x*y/x;
    }
    d_temp = t2-15;
    d_temp = d_temp*t2+15;
    d_temp = d_temp*t2-1;
    d_temp /= (t2+1);
    d_temp /= (t2+1);
    d_temp /= (t2+1);
    d_temp *= sign;
  }
//printf("### y=%f x=%f ret=%f ###",y,x,d_temp);
  return(d_temp);
}

/*-------------------------------------------------------------------*
 * malloc2d
 *-------------------------------------------------------------------*/
void *malloc2d(size_t size, int n1, int n2)
{
  int i;
  int t=size*n2;
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

/*-------------------------------------------------------------------*
 * malloc3d
 *-------------------------------------------------------------------*/
void *malloc3d(size_t size, int n1, int n2, int n3)
{
  int i,j;
  int t=size*n3;
  char ***a1, **a2, *a3;
  a1 = (char***)malloc((sizeof(*a1) + sizeof(**a1) * n2 + t * n2) * n1);
  if(a1){
    a2 = (char**)(a1 + n1);
    a3 = (char *)(a2 + n1 * n2);
    for(i=0; i<n1; i++){
      a1[i] = a2;
      for(j=0; j<n2; j++){
        a2[j] = a3;
        a3 += t;
      }
      a2 += n2;
    }
    return a1;
  }
  return NULL;
}

/*-------------------------------------------------------------------*
 * main
 *-------------------------------------------------------------------*/
int main
(
  int argc,
  char *argv[]
)
{
  int i,j,m;
  char opt;
  int size_L = 12;
  int seed;
  int size_Ns = size_L*size_L;
  int N_t = 1;/* # of temperature */
  int N_smp = 10000;
  int N_equ = 10000;
  int N_phys = 60;
  double T_i = 0.1;
  double T_f = 1.0;
  double T_width = 0.02;
  double parJ1 = -1.0;/* antiferro */
  double parD  = -1.0;/* antiferro */
  double time_count[10];

  int N_hist_mag = 120;
  int N_hist_theta = 120;
  int N_hist_cos6theta = 120;
//  int N_hist_mag = size_L*2;
//  int N_hist_theta = size_L*2;
//  int N_hist_cos6theta = size_L*2;
  double R_hist_mag = 2.0;
  double R_hist_theta = M_PI+1.0e-6;
  double R_hist_cos6theta = 1.0;

  int i_temp1,i_temp2;
  double d_temp[N_phys/5];

  time_count[0] = gettimeofday_sec();

  seed = (int)(time_count[0]
    + (time_count[0]-(int)time_count[0])*100000000);
//  printf("%f %d\n",time_count[0],seed);
//  seed=time(NULL);
//  seed=1234567;

  /* prameter settings */
  for(i = 0; i < argc; ++i){
    if(*argv[i] == '-'){
      opt = *(argv[i]+1);
      if(opt == 'l'){
        size_L = atoi(argv[i+1]);
        size_Ns = size_L*size_L;
      }else if(opt == 's'){
        N_smp = atoi(argv[i+1]);
      }else if(opt == 'e'){
        N_equ = atoi(argv[i+1]);
      }else if(opt == 'i'){
        T_i = atof(argv[i+1]);
      }else if(opt == 'f'){
        T_f = atof(argv[i+1]);
      }else if(opt == 'w'){
        T_width = atof(argv[i+1]);
      }else if(opt == 'd'){
        parD = atof(argv[i+1]);
      }else if(opt == 'r'){
        seed = atoi(argv[i+1]);
      }
    }
  }

  N_t = (int)((T_f-T_i)/T_width + 1+1.0e-10);

  int list_temp[N_t];
  int list_conf[N_t];
  double T[N_t];
  double Eadd[N_t];
  double **PhysQuant = (double**)malloc2d(sizeof(double),N_t,N_phys);
  double **Madd = (double**)malloc2d(sizeof(double),N_t,4);
  double **HistTheta = (double**)malloc2d(sizeof(double),N_t,(N_hist_theta+1));
  double **HistCos6Theta = (double**)malloc2d(sizeof(double),N_t,(N_hist_cos6theta+1));
  int ***spins = (int***)malloc3d(sizeof(int),N_t,size_L,size_L);
  double ***BoltzFac = (double***)malloc3d(sizeof(double),N_t,25,25);
  double ***HistMag = (double***)malloc3d(sizeof(double),N_t,(N_hist_mag+1),(N_hist_mag+1));

  FILE *fp_time;
  FILE *fp_para;
  FILE *fp_phys;
  FILE *fp_phys_chosen;
  FILE *fp_hist_mag;
  FILE *fp_hist_theta;
  FILE *fp_hist_cos6theta;
  char buf_time[256];
  char buf_para[256];
  char buf_phys[256];
  char buf_phys_chosen[256];
  char **buf_hist_mag = (char**)malloc2d(sizeof(char),N_t,256);
  char **buf_hist_theta = (char**)malloc2d(sizeof(char),N_t,256);
  char **buf_hist_cos6theta = (char**)malloc2d(sizeof(char),N_t,256);
  sprintf(buf_time,"ZZZ_time_seed_%010d.dat",seed);
  sprintf(buf_para,"ZZZ_para_seed_%010d.dat",seed);
  sprintf(buf_phys,"ZZZ_phys_seed_%010d.dat",seed);
  sprintf(buf_phys_chosen,"ZZZ_phys_chosen_seed_%010d.dat",seed);
  fp_time=fopen(buf_time,"w");
  fp_para=fopen(buf_para,"w");
  fp_phys=fopen(buf_phys,"w");
  fp_phys_chosen=fopen(buf_phys_chosen,"w");
/*
  for(m=0; m<N_t; m++){
    sprintf(buf_hist_mag[m],"ZZZ_hist_mag_seed_%010d_temp_%06d.dat",seed,m);
    sprintf(buf_hist_theta[m],"ZZZ_hist_theta_seed_%010d_temp_%06d.dat",seed,m);
    sprintf(buf_hist_cos6theta[m],"ZZZ_hist_cos6theta_seed_%010d_temp_%06d.dat",seed,m);
  }
*/

  /* set temperature */
  for(m=0; m<N_t; m++){
    list_temp[m] = m;
    list_conf[m] = m;
    T[m] = T_i + T_width * m;
  }
  /* end of set temperature */

  /* table of Boltzmann factors */
  for(m=0; m<N_t; m++){
    for(i=-12; i<=12; i++){
      for(j=-1; j<=1; j++){
        BoltzFac[m][i+12][j+1] = exp((parJ1*i+parD*j)/T[m]);
      }
    }
  }
  /* end of table of Boltzmann factors */

  for(m=0; m<N_t; m++){
    for(i=0; i<N_hist_mag+1; i++){
      for(j=0; j<N_hist_mag+1; j++){
        HistMag[m][i][j] = 0.0;
      }
    }
    for(i=0; i<N_hist_theta+1; i++){
      HistTheta[m][i] = 0.0;
    }
    for(i=0; i<N_hist_cos6theta+1; i++){
      HistCos6Theta[m][i] = 0.0;
    }
  }

  fprintf(fp_para,"# seed=%d\n",seed);
  fprintf(fp_para,"J1=%13.6e\n",parJ1);
  fprintf(fp_para,"D=%13.6e\n",parD);
  fprintf(fp_para,"size_L=%d\n",size_L);
  fprintf(fp_para,"size_Ns=%d\n",size_Ns);
  fprintf(fp_para,"N_smp=%d\n",N_smp);
  fprintf(fp_para,"N_equ=%d\n",N_equ);
  fprintf(fp_para,"T_i=%13.6e\n",T_i);
  fprintf(fp_para,"T_f=%13.6e\n",T_f);
  fprintf(fp_para,"T_width=%13.6e\n",T_width);
  fprintf(fp_para,"N_t=%d\n",N_t);
  fprintf(fp_para,"#\n");
  fprintf(fp_para,"# 1:     T\n");
  fprintf(fp_para,"# 2-6:   E\n");
  fprintf(fp_para,"# 7-11:  M\n");
  fprintf(fp_para,"# 12-16: Mx\n");
  fprintf(fp_para,"# 17-21: My\n");
  fprintf(fp_para,"# 22-26: Mz\n");
  fprintf(fp_para,"# 27-31: Mxy=sqrt(Mx^2+My^2)\n");
  fprintf(fp_para,"# 32-36: |Mz|\n");
  fprintf(fp_para,"# 37-41: Msub\n");
  fprintf(fp_para,"# 42-46: theta\n");
  fprintf(fp_para,"# 47-51: cos6atan\n");
  fprintf(fp_para,"# 52-56: Mxy^6\n");
  fprintf(fp_para,"# 57-61: Mxy^6 cos6atan\n");

  time_count[1] = gettimeofday_sec();

  dsfmt_init_gen_rand(&dsfmt,seed);

  /* initialization */
  for(m=0; m<N_t; m++){// loop for each temperature
    for(i=0; i<size_L; i++){
      for(j=0; j<size_L; j++){
        spins[m][i][j] = ((int)(dsfmt_genrand_close_open(&dsfmt) * 3)) - 1;
      }
    }
  }

  time_count[2] = gettimeofday_sec();

  /* neglecting and exchanging */
  for(i=0; i<N_equ; i++){
    for(m=0; m<N_t-1; m++){/* loop for each temperature */
      /* m: parameter of environment */
      calc_energy(size_L,list_conf[m],parJ1,parD,spins,Eadd);
      calc_energy(size_L,list_conf[m+1],parJ1,parD,spins,Eadd);
      exchange_states(size_L,m,parJ1,parD,T,list_temp,list_conf,spins,Eadd);
    }
    for(m=0; m<N_t; m++){/* loop for each temperature */
      /* m: parameter of environment */
      for(j=0; j<size_Ns; j++){
        one_flip(size_L,m,BoltzFac,list_temp,spins);
      }
    }
  }
  /* end of neglecting and exchanging */

  time_count[3] = gettimeofday_sec();

  /* parameter initialization */
  for(m=0; m<N_t; m++){/* loop for each temperature */
    /* m: parameter of environment */
    calc_energy(size_L,m,parJ1,parD,spins,Eadd);
    calc_mag(size_L,m,spins,Madd);
    for(i=0; i<N_phys; i++){
      PhysQuant[m][i] = 0.0;
    }
  }
  /* end of parameter initialization */

  time_count[4] = gettimeofday_sec();

  /* measurement loop */
  for(i=0; i<N_smp; i++){
    for(m=0; m<N_t-1; m++){/* loop for each temperature */
      /* m: parameter of environment */
      exchange_states(size_L,m,parJ1,parD,T,list_temp,list_conf,spins,Eadd);
    }
    for(m=0; m<N_t; m++){/* loop for each temperature */
      /* m: parameter of environment */
      for(j=0; j<size_Ns; j++){
        one_flip_ene_mag(size_L,m,parJ1,parD,BoltzFac,list_temp,spins,Eadd,Madd);
      }
    }

    for(m=0; m<N_t; m++){/* loop for each temperature */
      /* m:            parameter of environment */
      /* list_temp[m]: parameter of temperature */
      d_temp[0] = Eadd[m];
      d_temp[1] = Madd[m][0];
      d_temp[2] = 3.0*(2.0*Madd[m][1]-Madd[m][2]-Madd[m][3])/sqrt(6.0);// Mx
      d_temp[3] = 3.0*(Madd[m][2]-Madd[m][3])/sqrt(2.0);// My
      d_temp[4] = 3.0*(Madd[m][1]+Madd[m][2]+Madd[m][3])/sqrt(3.0);// Mz
      d_temp[5] = sqrt(d_temp[2]*d_temp[2] + d_temp[3]*d_temp[3]);// Mxy=sqrt(Mx^2+My^2)
      d_temp[6] = fabs(d_temp[4]);// |Mz|
      d_temp[7] = 3.0*sqrt(Madd[m][1]*Madd[m][1]+Madd[m][2]*Madd[m][2]
        +Madd[m][3]*Madd[m][3]);// Msub
      d_temp[8] = atan2(d_temp[3],d_temp[2]);// theta
      d_temp[9] = cos(6.0*d_temp[8]);// cos6atan
//      d_temp8 = cos6atan(d_temp2,d_temp1);// cos6atan
      d_temp[10] = pow(d_temp[5],6.0);// Mxy^6
      d_temp[11] = d_temp[9]*d_temp[10];// Mxy^6 cos6atan

      for(j=0; j<N_phys/5; j++){
        calc_phys_quant(m,5*j,list_temp,PhysQuant,d_temp[j]);
      }

/*
//// index = (int) 0.5*(Nx/R+N+1)
//// x in [-R,R]
//// index in {0,1,2,...,N} (N+1 components)

      i_temp1 = (int)(0.5*(N_hist_mag*(d_temp[2]/size_Ns)/R_hist_mag+N_hist_mag+1));
      i_temp2 = (int)(0.5*(N_hist_mag*(d_temp[3]/size_Ns)/R_hist_mag+N_hist_mag+1));
      HistMag[list_temp[m]][i_temp1][i_temp2] += 1.0;

      i_temp1 = (int)(0.5*(N_hist_theta*d_temp[8]/R_hist_theta+N_hist_theta+1));
      HistTheta[list_temp[m]][i_temp1] += 1.0;

      i_temp1 = (int)(0.5*(N_hist_cos6theta*d_temp[9]/R_hist_cos6theta+N_hist_cos6theta+1));
      HistCos6Theta[list_temp[m]][i_temp1] += 1.0;
*/

    }
  }
  /* end of measurement loop */

  time_count[5] = gettimeofday_sec();

  for(m=0; m<N_t; m++){/* loop for each temperature */
    /* m: parameter of temperature */
    /*                 ^^^^^^^^^^^ */
    for(i=0; i<N_phys/5; i++){// Ene, Mag, |Mag|, ...
      for(j=0; j<3; j++){// Phys, Phys2, Phys4
        PhysQuant[m][5*i+j] /= (double)N_smp;
      }
      // susceptibility
      PhysQuant[m][5*i+3]
        = PhysQuant[m][5*i+1]-PhysQuant[m][5*i+0]*PhysQuant[m][5*i+0];
      // Binder parameter
      PhysQuant[m][5*i+4]
        = 1.0-PhysQuant[m][5*i+2]/PhysQuant[m][5*i+1]/PhysQuant[m][5*i+1]/3.0;
    }
    for(i=0; i<N_phys; i++){
      PhysQuant[m][i] /= (double)size_Ns;
    }
  }

  time_count[6] = gettimeofday_sec();

  for(m=0; m<N_t; m++){/* loop for each temperature */
    /* m: parameter of temperature */
    /*                 ^^^^^^^^^^^ */
    fprintf(fp_phys,"%13.6e ",T[m]);
    for(i=0; i<N_phys; i++){
      fprintf(fp_phys,"%13.6e ",PhysQuant[m][i]);
    }
    fprintf(fp_phys,"\n");

    fprintf(fp_phys_chosen,"%13.6e ",T[m]);
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][3]/T[m]/T[m]);// C
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][8]/T[m]);// X of m
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][28]/T[m]);// X of sqrt(Mx^2+My^2)
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][33]/T[m]);// X of |Mz|
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][38]/T[m]);// X of Msub
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][26]);// X' of sqrt(Mx^2+My^2)
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][31]);// X' of |Mz|
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][36]);// X' of Msub
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][45]*size_Ns);// cos6atan
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][55]/PhysQuant[m][50]);// <Mxy^6 cos6atan>/<Mxy^6>
    fprintf(fp_phys_chosen,"\n");
  }

/*
  for(m=0; m<N_t; m++){// loop for each temperature
    fp_hist_mag=fopen(buf_hist_mag[m],"w");
    fp_hist_theta=fopen(buf_hist_theta[m],"w");
    fp_hist_cos6theta=fopen(buf_hist_cos6theta[m],"w");

    fprintf(fp_hist_mag,"# seed=%d\n",seed);
    fprintf(fp_hist_mag,"# %13.6e\n",T[m]);
    for(i=0; i<N_hist_mag+1; i++){
      for(j=0; j<N_hist_mag+1; j++){
        fprintf(fp_hist_mag,"%3d %3d %13.6e %13.6e %13.6e\n",
          i,j,
          (i*2.0/N_hist_mag-1.0)*R_hist_mag,
          (j*2.0/N_hist_mag-1.0)*R_hist_mag,
          HistMag[m][i][j]/(double)N_smp);
      }
      fprintf(fp_hist_mag,"\n");
    }

    fprintf(fp_hist_theta,"# seed=%d\n",seed);
    fprintf(fp_hist_theta,"# %13.6e\n",T[m]);
    for(i=0; i<N_hist_theta+1; i++){
      fprintf(fp_hist_theta,"%3d %13.6e %13.6e\n",
        i,(i*2.0/N_hist_theta-1.0)*R_hist_theta,
        HistTheta[m][i]/(double)N_smp);
    }

    fprintf(fp_hist_cos6theta,"# seed=%d\n",seed);
    fprintf(fp_hist_cos6theta,"# %13.6e\n",T[m]);
    for(i=0; i<N_hist_cos6theta+1; i++){
      fprintf(fp_hist_cos6theta,"%3d %13.6e %13.6e\n",
        i,(i*2.0/N_hist_cos6theta-1.0)*R_hist_cos6theta,
        HistCos6Theta[m][i]/(double)N_smp);
    }

    fclose(fp_hist_mag);
    fclose(fp_hist_theta);
    fclose(fp_hist_cos6theta);
  }
*/

  time_count[7] = gettimeofday_sec();

  fprintf(fp_time,"# seed=%d\n",seed);
  fprintf(fp_time,"total time: %13.6f\n",time_count[7]-time_count[0]);
  fprintf(fp_time,"  parameter setting:           %13.6f\n",time_count[1]-time_count[0]);
  fprintf(fp_time,"  initialization (spin):       %13.6f\n",time_count[2]-time_count[1]);
  fprintf(fp_time,"  run until equilibrium:       %13.6f\n",time_count[3]-time_count[2]);
  fprintf(fp_time,"  initialization (parameters): %13.6f\n",time_count[4]-time_count[3]);
  fprintf(fp_time,"  measurement:                 %13.6f\n",time_count[5]-time_count[4]);
  fprintf(fp_time,"  calculate suscep. etc:       %13.6f\n",time_count[6]-time_count[5]);
  fprintf(fp_time,"  print results:               %13.6f\n",time_count[7]-time_count[6]);

  fclose(fp_time);
  fclose(fp_para);
  fclose(fp_phys);
  fclose(fp_phys_chosen);

  free(PhysQuant);
  free(Madd);
  free(spins);
  free(BoltzFac);
  free(HistMag);
  free(HistTheta);
  free(HistCos6Theta);
  free(buf_hist_mag);
  free(buf_hist_theta);
  free(buf_hist_cos6theta);

  return 0;
}
