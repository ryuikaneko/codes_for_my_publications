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
#define max(a,b) (((a)<(b)) ? (b):(a))
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
int calc_energy(
  int size_L,
  int m,
  double parU,
  double parV,
  double parmu,
  int ***n_up,
  int ***n_dn,
  double *ene
)
{
  int i,ip,j,jp;
  double tmp_ene;
  tmp_ene = 0.0;
  for(i=0; i<size_L; i++){
    ip  = (i+1)%size_L;
    for(j=0; j<size_L; j++){
      jp  = (j+1)%size_L;
      tmp_ene +=
        + parU * n_up[m][i][j] * n_dn[m][i][j]
        + parV * (n_up[m][i][j] + n_dn[m][i][j])
          *(n_up[m][ip][j] + n_dn[m][ip][j]
          + n_up[m][i][jp] + n_dn[m][i][jp]
          + n_up[m][ip][jp] + n_dn[m][ip][jp])
        - parmu * (n_up[m][i][j] + n_dn[m][i][j]);
    }
  }
  ene[m] = tmp_ene;
  return 0;
}

/*-------------------------------------------------------------------*
 * the number of electrons
 *-------------------------------------------------------------------*/
int calc_num(
  int size_L,
  int m,
  int ***n_up,
  int ***n_dn,
  double **num
)
{
  int i,j;
  double tmp_num[4];
  tmp_num[0] = 0.0;
  tmp_num[1] = 0.0;
  tmp_num[2] = 0.0;
  tmp_num[3] = 0.0;
  for(i=0; i<size_L; i++){
    for(j=0; j<size_L; j++){
      tmp_num[0] += n_up[m][i][j] + n_dn[m][i][j];
      tmp_num[1+(i+j)%3] += n_up[m][i][j] + n_dn[m][i][j];// 3 sublattice
    }
  }
  num[m][0] = tmp_num[0];
  num[m][1] = tmp_num[1];
  num[m][2] = tmp_num[2];
  num[m][3] = tmp_num[3];
  return 0;
}

/*-------------------------------------------------------------------*
 * charge and spin corelations
 *-------------------------------------------------------------------*/
int calc_corr(
  int size_L,
  int m,
  int ***n_up,
  int ***n_dn,
  double ***corr_c,
  double ***corr_s
)
{
  int i,j;
  for(i=0; i<size_L; i++){
    for(j=0; j<size_L; j++){
      // N(ri,rj) = N[(0,0),(xj-xi,yj-yi)]
      // S(ri,rj) = S[(0,0),(xj-xi,yj-yi)]
      corr_c[m][i][j] = (n_up[m][0][0] + n_dn[m][0][0])
        * (n_up[m][i][j] + n_dn[m][i][j]);
      corr_s[m][i][j] = (n_up[m][0][0] - n_dn[m][0][0])
        * (n_up[m][i][j] - n_dn[m][i][j]) * 0.25;
    }
  }
  return 0;
}

/*-------------------------------------------------------------------*
 * one monte carlo step (up)
 *-------------------------------------------------------------------*/
int one_flip_up(
  int size_L,
  int m,
  double ****BoltzFac,
  int *list_temp,
  int ***n_up,
  int ***n_dn
)
{
  int i,ip,im,j,jp,jm;
  int n_s,n_ms,n_nn;
  /* spin flip candidate */
  i = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  j = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  /* calc egap */
  ip  = (i+1)%size_L;
  im  = (i-1+size_L)%size_L;
  jp  = (j+1)%size_L;
  jm  = (j-1+size_L)%size_L;
  n_s  = n_up[m][i][j];
  n_ms = n_dn[m][i][j];
  n_nn =
    + n_up[m][ip][j] + n_dn[m][ip][j]
    + n_up[m][i][jp] + n_dn[m][i][jp]
    + n_up[m][im][j] + n_dn[m][im][j]
    + n_up[m][i][jm] + n_dn[m][i][jm]
    + n_up[m][ip][jp] + n_dn[m][ip][jp]
    + n_up[m][im][jm] + n_dn[m][im][jm];
  /* metropolis method */
  if(dsfmt_genrand_close_open(&dsfmt) <= BoltzFac[list_temp[m]][n_s][n_ms][n_nn]){
    n_up[m][i][j] = 1 - n_up[m][i][j];
  }
  return 0;
}

/*-------------------------------------------------------------------*
 * one monte carlo step (down)
 *-------------------------------------------------------------------*/
int one_flip_dn(
  int size_L,
  int m,
  double ****BoltzFac,
  int *list_temp,
  int ***n_up,
  int ***n_dn
)
{
  int i,ip,im,j,jp,jm;
  int n_s,n_ms,n_nn;
  /* spin flip candidate */
  i = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  j = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  /* calc egap */
  ip  = (i+1)%size_L;
  im  = (i-1+size_L)%size_L;
  jp  = (j+1)%size_L;
  jm  = (j-1+size_L)%size_L;
  n_s  = n_dn[m][i][j];
  n_ms = n_up[m][i][j];
  n_nn =
    + n_up[m][ip][j] + n_dn[m][ip][j]
    + n_up[m][i][jp] + n_dn[m][i][jp]
    + n_up[m][im][j] + n_dn[m][im][j]
    + n_up[m][i][jm] + n_dn[m][i][jm]
    + n_up[m][ip][jp] + n_dn[m][ip][jp]
    + n_up[m][im][jm] + n_dn[m][im][jm];
  /* metropolis method */
  if(dsfmt_genrand_close_open(&dsfmt) <= BoltzFac[list_temp[m]][n_s][n_ms][n_nn]){
    n_dn[m][i][j] = 1 - n_dn[m][i][j];
  }
  return 0;
}

/*-------------------------------------------------------------------*
 * one monte carlo step (up)
 * (calculate energy and charge differences)
 *-------------------------------------------------------------------*/
int one_flip_phys_up(
  int size_L,
  int m,
  double ****BoltzFac,
  double ****Egap,
  int *list_temp,
  int ***n_up,
  int ***n_dn,
  double *ene,
  double **num
)
{
  int i,ip,im,j,jp,jm;
  int n_s,n_ms,n_nn;
  /* spin flip candidate */
  i = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  j = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  /* calc egap */
  ip  = (i+1)%size_L;
  im  = (i-1+size_L)%size_L;
  jp  = (j+1)%size_L;
  jm  = (j-1+size_L)%size_L;
  n_s  = n_up[m][i][j];
  n_ms = n_dn[m][i][j];
  n_nn =
    + n_up[m][ip][j] + n_dn[m][ip][j]
    + n_up[m][i][jp] + n_dn[m][i][jp]
    + n_up[m][im][j] + n_dn[m][im][j]
    + n_up[m][i][jm] + n_dn[m][i][jm]
    + n_up[m][ip][jp] + n_dn[m][ip][jp]
    + n_up[m][im][jm] + n_dn[m][im][jm];
  /* metropolis method */
  if(dsfmt_genrand_close_open(&dsfmt) <= BoltzFac[list_temp[m]][n_s][n_ms][n_nn]){
    n_up[m][i][j] = 1 - n_up[m][i][j];
    ene[m] += Egap[m][n_s][n_ms][n_nn];
    num[m][0] += 2*n_up[m][i][j] - 1;
    num[m][1+(i+j)%3] += 2*n_up[m][i][j] - 1;// 3 sublattice
  }
  return 0;
}


/*-------------------------------------------------------------------*
 * one monte carlo step (down)
 * (calculate energy and charge differences)
 *-------------------------------------------------------------------*/
int one_flip_phys_dn(
  int size_L,
  int m,
  double ****BoltzFac,
  double ****Egap,
  int *list_temp,
  int ***n_up,
  int ***n_dn,
  double *ene,
  double **num
)
{
  int i,ip,im,j,jp,jm;
  int n_s,n_ms,n_nn;
  /* spin flip candidate */
  i = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  j = (int) (size_L * dsfmt_genrand_close_open(&dsfmt));
  /* calc egap */
  ip  = (i+1)%size_L;
  im  = (i-1+size_L)%size_L;
  jp  = (j+1)%size_L;
  jm  = (j-1+size_L)%size_L;
  n_s  = n_dn[m][i][j];
  n_ms = n_up[m][i][j];
  n_nn =
    + n_up[m][ip][j] + n_dn[m][ip][j]
    + n_up[m][i][jp] + n_dn[m][i][jp]
    + n_up[m][im][j] + n_dn[m][im][j]
    + n_up[m][i][jm] + n_dn[m][i][jm]
    + n_up[m][ip][jp] + n_dn[m][ip][jp]
    + n_up[m][im][jm] + n_dn[m][im][jm];
  /* metropolis method */
  if(dsfmt_genrand_close_open(&dsfmt) <= BoltzFac[list_temp[m]][n_s][n_ms][n_nn]){
    n_dn[m][i][j] = 1 - n_dn[m][i][j];
    ene[m] += Egap[m][n_s][n_ms][n_nn];
    num[m][0] += 2*n_dn[m][i][j] - 1;
    num[m][1+(i+j)%3] += 2*n_dn[m][i][j] - 1;// 3 sublattice
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
  double parU,
  double parV,
  double parmu,
  double *parT,
  int *list_temp,
  int *list_conf,
  int ***n_up,
  int ***n_dn,
  double *ene
)
{
  int conf0,conf1;
  double egap,w;
  /* calculate gap */
  conf0 = list_conf[m];
  conf1 = list_conf[m+1];
  egap = (1.0/parT[m+1] - 1.0/parT[m]) * (ene[conf0] - ene[conf1]);
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
 * calculate correlations during sampling
 *-------------------------------------------------------------------*/
int calc_phys_quant_corr
(
  int size_L,
  int m,
  int *list_temp,
  double ***PhysQuant,
  double ***corr
)
{
  int i,j;
  for(i=0; i<size_L; i++){
    for(j=0; j<size_L; j++){
      PhysQuant[list_temp[m]][i][j] += corr[m][i][j];
    }
  }
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
 * malloc4d
 *-------------------------------------------------------------------*/
void *malloc4d(size_t size, int n1, int n2, int n3, int n4)
{
  int i,j,k;
  int  t=size*n4;
  char ****a1,***a2,**a3,*a4;
  a1 = (char****)malloc((sizeof(*a1)+sizeof(**a1)*n2+sizeof(***a1)*n2*n3+t*n2*n3)*n1);
  if(a1){
    a2 = (char***)(a1 + n1);
    a3 = (char **)(a2 + n1 * n2);
    a4 = (char  *)(a3 + n1 * n2 * n3);
    for(i=0; i<n1; i++){
      a1[i] = a2;
      for(j=0; j<n2; j++){
        a2[j] = a3;
        for(k=0; k<n3; k++){
          a3[k] = a4;
          a4 += t;
        }
        a3 += n3;
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
  int i,j,k,m;
  char opt;
  int size_L = 12;
  int seed;
  int Ns = size_L*size_L;
  int N_t = 1;/* # of temperature */
  int N_smp = 10000;
  int N_equ = 10000;
  int N_phys = 60;
  double parTi = 0.1;
  double parTf = 2.0;
  double parTw = 0.1;
  double parU = 1.0;
  double parV = 1.0;
  int parz = 6;/* 6 n.n. for a triangular lattice */
  double parmu0;
  double parmu = 0.0;
  double time_count[10];

  int N_hist_num = 120;
  int N_hist_theta = 120;
  int N_hist_cos6theta = 120;
  double R_hist_num = 2.0;
  double R_hist_theta = M_PI+1.0e-6;
  double R_hist_cos6theta = 1.0;

  int i_temp1,i_temp2;
  double d_temp[N_phys/5];

  int ix,iy;
  int kx,ky;
  double qx,qy;
  double inv_size_L;
  double inv_Ns;
  double ave_n2;

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
        Ns = size_L*size_L;
      }else if(opt == 's'){
        N_smp = atoi(argv[i+1]);
      }else if(opt == 'e'){
        N_equ = atoi(argv[i+1]);
      }else if(opt == 'i'){
        parTi = atof(argv[i+1]);
      }else if(opt == 'f'){
        parTf = atof(argv[i+1]);
      }else if(opt == 'w'){
        parTw = atof(argv[i+1]);
      }else if(opt == 'u'){
        parU = atof(argv[i+1]);
      }else if(opt == 'v'){
        parV = atof(argv[i+1]);
      }else if(opt == 'm'){
        parmu = atof(argv[i+1]);
      }else if(opt == 'r'){
        seed = atoi(argv[i+1]);
      }
    }
  }

  parmu0 = 0.5*parU + parz*parV;/* chemical potential @ n=1 */
  parmu += parmu0;

  N_t = (int)((parTf-parTi)/parTw + 1+1.0e-10);

  int list_temp[N_t];
  int list_conf[N_t];
  double parT[N_t];
  double Eadd[N_t];
  double **PhysQuant = (double**)malloc2d(sizeof(double),N_t,N_phys);
  double **Nadd = (double**)malloc2d(sizeof(double),N_t,4);
  double **HistTheta = (double**)malloc2d(sizeof(double),N_t,(N_hist_theta+1));
  double **HistCos6Theta = (double**)malloc2d(sizeof(double),N_t,(N_hist_cos6theta+1));
  int ***n_up = (int***)malloc3d(sizeof(int),N_t,size_L,size_L);
  int ***n_dn = (int***)malloc3d(sizeof(int),N_t,size_L,size_L);
  double ***HistNum = (double***)malloc3d(sizeof(double),N_t,(N_hist_num+1),(N_hist_num+1));
  double ****BoltzFac = (double****)malloc4d(sizeof(double),N_t,2,2,13);
  double ****Egap = (double****)malloc4d(sizeof(double),N_t,2,2,13);

  double ****cosqr = (double****)malloc4d(sizeof(double),size_L,size_L,size_L,size_L);
  double ****sinqr = (double****)malloc4d(sizeof(double),size_L,size_L,size_L,size_L);

  double ***corr_c = (double***)malloc3d(sizeof(double),N_t,size_L,size_L);
  double ***corr_s = (double***)malloc3d(sizeof(double),N_t,size_L,size_L);
  double ***PhysQuantNr = (double***)malloc3d(sizeof(double),N_t,size_L,size_L);
  double ***PhysQuantSr = (double***)malloc3d(sizeof(double),N_t,size_L,size_L);
  double ***PhysQuantShiftNr = (double***)malloc3d(sizeof(double),N_t,size_L,size_L);
  double ***PhysQuantNqRe = (double***)malloc3d(sizeof(double),N_t,size_L,size_L);
  double ***PhysQuantNqIm = (double***)malloc3d(sizeof(double),N_t,size_L,size_L);
  double ***PhysQuantSqRe = (double***)malloc3d(sizeof(double),N_t,size_L,size_L);
  double ***PhysQuantSqIm = (double***)malloc3d(sizeof(double),N_t,size_L,size_L);
  double ***PhysQuantShiftNqRe = (double***)malloc3d(sizeof(double),N_t,size_L,size_L);
  double ***PhysQuantShiftNqIm = (double***)malloc3d(sizeof(double),N_t,size_L,size_L);

  double enegap;
  int intgap_U;
  int intgap_V;
  int intgap_mu;

  FILE *fp_time;
  FILE *fp_para;
  FILE *fp_phys;
  FILE *fp_phys_chosen;
  FILE *fp_phys_nr;
  FILE *fp_phys_sr;
  FILE *fp_phys_shiftnr;
  FILE *fp_phys_nq;
  FILE *fp_phys_sq;
  FILE *fp_phys_shiftnq;
  FILE *fp_hist_num;
  FILE *fp_hist_theta;
  FILE *fp_hist_cos6theta;
  char buf_time[256];
  char buf_para[256];
  char buf_phys[256];
  char buf_phys_chosen[256];
  char **buf_hist_num = (char**)malloc2d(sizeof(char),N_t,256);
  char **buf_hist_theta = (char**)malloc2d(sizeof(char),N_t,256);
  char **buf_hist_cos6theta = (char**)malloc2d(sizeof(char),N_t,256);
  char **buf_phys_nr = (char**)malloc2d(sizeof(char),N_t,256);
  char **buf_phys_sr = (char**)malloc2d(sizeof(char),N_t,256);
  char **buf_phys_shiftnr = (char**)malloc2d(sizeof(char),N_t,256);
  char **buf_phys_nq = (char**)malloc2d(sizeof(char),N_t,256);
  char **buf_phys_sq = (char**)malloc2d(sizeof(char),N_t,256);
  char **buf_phys_shiftnq = (char**)malloc2d(sizeof(char),N_t,256);
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
    sprintf(buf_hist_num[m],"ZZZ_hist_num_seed_%010d_temp_%06d.dat",seed,m);
    sprintf(buf_hist_theta[m],"ZZZ_hist_theta_seed_%010d_temp_%06d.dat",seed,m);
    sprintf(buf_hist_cos6theta[m],"ZZZ_hist_cos6theta_seed_%010d_temp_%06d.dat",seed,m);
  }
*/
  for(m=0; m<N_t; m++){
    sprintf(buf_phys_nr[m],"ZZZ_phys_nr_seed_%010d_temp_%06d.dat",seed,m);
    sprintf(buf_phys_sr[m],"ZZZ_phys_sr_seed_%010d_temp_%06d.dat",seed,m);
    sprintf(buf_phys_shiftnr[m],"ZZZ_phys_shiftnr_seed_%010d_temp_%06d.dat",seed,m);
    sprintf(buf_phys_nq[m],"ZZZ_phys_nq_seed_%010d_temp_%06d.dat",seed,m);
    sprintf(buf_phys_sq[m],"ZZZ_phys_sq_seed_%010d_temp_%06d.dat",seed,m);
    sprintf(buf_phys_shiftnq[m],"ZZZ_phys_shiftnq_seed_%010d_temp_%06d.dat",seed,m);
  }

  /* set temperature */
  for(m=0; m<N_t; m++){
    list_temp[m] = m;
    list_conf[m] = m;
    parT[m] = parTi + parTw * m;
  }
  /* end of set temperature */

  /* table of Boltzmann factors */
  for(m=0; m<N_t; m++){
    for(i=0; i<=1; i++){
      /* 1-2*n_target_spin */
      intgap_mu = 1-2*i;
      for(j=0; j<=1; j++){
        /* (1-2*n_target_spin)*n_opposite_spin */
        intgap_U = intgap_mu*j;
        for(k=0; k<=12; k++){
          intgap_V = intgap_mu*k;
          enegap = parU*intgap_U + parV*intgap_V - parmu*intgap_mu;
          Egap[m][i][j][k] = enegap;
          BoltzFac[m][i][j][k] = exp(-enegap/parT[m]);
        }
      }
    }
  }
  /* end of table of Boltzmann factors */

  for(m=0; m<N_t; m++){
    for(i=0; i<N_hist_num+1; i++){
      for(j=0; j<N_hist_num+1; j++){
        HistNum[m][i][j] = 0.0;
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
  fprintf(fp_para,"parU=%13.6e\n",parU);
  fprintf(fp_para,"parV=%13.6e\n",parV);
  fprintf(fp_para,"parmu=%13.6e\n",parmu);
  fprintf(fp_para,"size_L=%d\n",size_L);
  fprintf(fp_para,"Ns=%d\n",Ns);
  fprintf(fp_para,"N_smp=%d\n",N_smp);
  fprintf(fp_para,"N_equ=%d\n",N_equ);
  fprintf(fp_para,"parTi=%13.6e\n",parTi);
  fprintf(fp_para,"parTf=%13.6e\n",parTf);
  fprintf(fp_para,"parTw=%13.6e\n",parTw);
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
        n_up[m][i][j] = (dsfmt_genrand_close_open(&dsfmt) > 0.5) ? 0 : 1;
        n_dn[m][i][j] = (dsfmt_genrand_close_open(&dsfmt) > 0.5) ? 0 : 1;
      }
    }
  }

  time_count[2] = gettimeofday_sec();

  /* neglecting and exchanging */
  for(i=0; i<N_equ; i++){
    for(m=0; m<N_t-1; m++){/* loop for each temperature */
      /* m: parameter of environment */
      calc_energy(size_L,list_conf[m],parU,parV,parmu,n_up,n_dn,Eadd);
      calc_energy(size_L,list_conf[m+1],parU,parV,parmu,n_up,n_dn,Eadd);
      exchange_states(size_L,m,parU,parV,parmu,parT,list_temp,list_conf,n_up,n_dn,Eadd);
    }
    for(m=0; m<N_t; m++){/* loop for each temperature */
      /* m: parameter of environment */
      for(j=0; j<Ns; j++){
        one_flip_up(size_L,m,BoltzFac,list_temp,n_up,n_dn);
        one_flip_dn(size_L,m,BoltzFac,list_temp,n_up,n_dn);
      }
    }
  }
  /* end of neglecting and exchanging */

  time_count[3] = gettimeofday_sec();

  /* parameter initialization */
  for(m=0; m<N_t; m++){/* loop for each temperature */
    /* m: parameter of environment */
    calc_energy(size_L,m,parU,parV,parmu,n_up,n_dn,Eadd);
    calc_num(size_L,m,n_up,n_dn,Nadd);
    for(i=0; i<N_phys; i++){
      PhysQuant[m][i] = 0.0;
    }
    for(i=0; i<size_L; i++){
      for(j=0; j<size_L; j++){
        PhysQuantNr[m][i][j] = 0.0;
        PhysQuantSr[m][i][j] = 0.0;
        PhysQuantShiftNr[m][i][j] = 0.0;
      }
    }
  }
  /* end of parameter initialization */

  time_count[4] = gettimeofday_sec();

  /* measurement loop */
  for(i=0; i<N_smp; i++){
    for(m=0; m<N_t-1; m++){/* loop for each temperature */
      /* m: parameter of environment */
      exchange_states(size_L,m,parU,parV,parmu,parT,list_temp,list_conf,n_up,n_dn,Eadd);
    }
    for(m=0; m<N_t; m++){/* loop for each temperature */
      /* m: parameter of environment */
      for(j=0; j<Ns; j++){
        one_flip_phys_up(size_L,m,BoltzFac,Egap,list_temp,n_up,n_dn,Eadd,Nadd);
        one_flip_phys_dn(size_L,m,BoltzFac,Egap,list_temp,n_up,n_dn,Eadd,Nadd);
      }
    }
    for(m=0; m<N_t; m++){/* loop for each temperature */
      calc_corr(size_L,m,n_up,n_dn,corr_c,corr_s);
    }

    for(m=0; m<N_t; m++){/* loop for each temperature */
      /* m:            parameter of environment */
      /* list_temp[m]: parameter of temperature */
      d_temp[0] = Eadd[m];
      d_temp[1] = Nadd[m][0];
      d_temp[2] = 3.0*(2.0*Nadd[m][1]-Nadd[m][2]-Nadd[m][3])/sqrt(6.0);// Mx
      d_temp[3] = 3.0*(Nadd[m][2]-Nadd[m][3])/sqrt(2.0);// My
      d_temp[4] = 3.0*(Nadd[m][1]+Nadd[m][2]+Nadd[m][3])/sqrt(3.0);// Mz
      d_temp[5] = sqrt(d_temp[2]*d_temp[2] + d_temp[3]*d_temp[3]);// Mxy=sqrt(Mx^2+My^2)
      d_temp[6] = fabs(d_temp[4]);// |Mz|
      d_temp[7] = 3.0*sqrt(Nadd[m][1]*Nadd[m][1]+Nadd[m][2]*Nadd[m][2]
        +Nadd[m][3]*Nadd[m][3]);// Msub
      d_temp[8] = atan2(d_temp[3],d_temp[2]);// theta
      d_temp[9] = cos(6.0*d_temp[8]);// cos6atan
//      d_temp8 = cos6atan(d_temp2,d_temp1);// cos6atan
      d_temp[10] = pow(d_temp[5],6.0);// Mxy^6
      d_temp[11] = d_temp[9]*d_temp[10];// Mxy^6 cos6atan

      for(j=0; j<N_phys/5; j++){
        calc_phys_quant(m,5*j,list_temp,PhysQuant,d_temp[j]);
      }

      calc_phys_quant_corr(size_L,m,list_temp,PhysQuantNr,corr_c);
      calc_phys_quant_corr(size_L,m,list_temp,PhysQuantSr,corr_s);

/*

//// index = (int) 0.5*(Nx/R+N+1)
//// x in [-R,R]
//// index in {0,1,2,...,N} (N+1 components)

      i_temp1 = (int)(0.5*(N_hist_num*(d_temp[2]/Ns)/R_hist_num+N_hist_num+1));
      i_temp2 = (int)(0.5*(N_hist_num*(d_temp[3]/Ns)/R_hist_num+N_hist_num+1));
      HistNum[list_temp[m]][i_temp1][i_temp2] += 1.0;

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
    for(i=0; i<N_phys/5; i++){// Ene, Num, |Num|, ...
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
      PhysQuant[m][i] /= (double)Ns;
    }
    for(i=0; i<size_L; i++){
      for(j=0; j<size_L; j++){
        PhysQuantNr[m][i][j] /= (double)N_smp;
        PhysQuantSr[m][i][j] /= (double)N_smp;
      }
    }
    ave_n2 = PhysQuant[m][5]*PhysQuant[m][5];
    for(i=0; i<size_L; i++){
      for(j=0; j<size_L; j++){
        PhysQuantShiftNr[m][i][j] = PhysQuantNr[m][i][j] - ave_n2;
      }
    }
  }

  inv_size_L = 1.0/size_L;
  inv_Ns = 1.0/Ns;
  for(ix=0; ix<size_L; ix++){
    for(iy=0; iy<size_L; iy++){
      for(kx=0; kx<size_L; kx++){
        qx = 2.0*M_PI*inv_size_L*kx;
        for(ky=0; ky<size_L; ky++){
          qy = 2.0*M_PI*inv_size_L*ky;
          cosqr[ix][iy][kx][ky] = cos(ix*qx+iy*qy);
          sinqr[ix][iy][kx][ky] = sin(ix*qx+iy*qy);
        }
      }
    }
  }

  for(m=0; m<N_t; m++){/* loop for each temperature */
    for(kx=0; kx<size_L; kx++){
      for(ky=0; ky<size_L; ky++){
        PhysQuantNqRe[m][kx][ky] = 0.0;
        PhysQuantNqIm[m][kx][ky] = 0.0;
        PhysQuantSqRe[m][kx][ky] = 0.0;
        PhysQuantSqIm[m][kx][ky] = 0.0;
      }
    }
    for(ix=0; ix<size_L; ix++){
      for(iy=0; iy<size_L; iy++){
        for(kx=0; kx<size_L; kx++){
          for(ky=0; ky<size_L; ky++){
            PhysQuantNqRe[m][kx][ky] += PhysQuantNr[m][ix][iy]
              * cosqr[ix][iy][kx][ky] * inv_Ns;
            PhysQuantNqIm[m][kx][ky] += PhysQuantNr[m][ix][iy]
              * sinqr[ix][iy][kx][ky] * inv_Ns;
            PhysQuantSqRe[m][kx][ky] += PhysQuantSr[m][ix][iy]
              * cosqr[ix][iy][kx][ky] * inv_Ns;
            PhysQuantSqIm[m][kx][ky] += PhysQuantSr[m][ix][iy]
              * sinqr[ix][iy][kx][ky] * inv_Ns;
          }
        }
      }
    }
    for(kx=0; kx<size_L; kx++){
      for(ky=0; ky<size_L; ky++){
        PhysQuantShiftNqRe[m][kx][ky] = PhysQuantNqRe[m][kx][ky];
        PhysQuantShiftNqIm[m][kx][ky] = PhysQuantNqIm[m][kx][ky];
      }
    }
    PhysQuantShiftNqRe[m][0][0] = 0.0;
    PhysQuantShiftNqIm[m][0][0] = 0.0;
  }

  time_count[6] = gettimeofday_sec();

  for(m=0; m<N_t; m++){/* loop for each temperature */
    /* m: parameter of temperature */
    /*                 ^^^^^^^^^^^ */
    fprintf(fp_phys,"%13.6e ",parT[m]);
    for(i=0; i<N_phys; i++){
      fprintf(fp_phys,"%13.6e ",PhysQuant[m][i]);
    }
    fprintf(fp_phys,"\n");

    fprintf(fp_phys_chosen,"%13.6e ",parT[m]);
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][3]/parT[m]/parT[m]);// C
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][8]/parT[m]);// X of m
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][28]/parT[m]);// X of sqrt(Mx^2+My^2)
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][33]/parT[m]);// X of |Mz|
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][38]/parT[m]);// X of Msub
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][26]);// X' of sqrt(Mx^2+My^2)
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][31]);// X' of |Mz|
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][36]);// X' of Msub
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][45]*Ns);// cos6atan
    fprintf(fp_phys_chosen,"%13.6e ",PhysQuant[m][55]/PhysQuant[m][50]);// <Mxy^6 cos6atan>/<Mxy^6>
    fprintf(fp_phys_chosen,"\n");
  }

/*
  for(m=0; m<N_t; m++){// loop for each temperature
    fp_hist_num=fopen(buf_hist_num[m],"w");
    fp_hist_theta=fopen(buf_hist_theta[m],"w");
    fp_hist_cos6theta=fopen(buf_hist_cos6theta[m],"w");

    fprintf(fp_hist_num,"# seed=%d\n",seed);
    fprintf(fp_hist_num,"# %13.6e\n",parT[m]);
    for(i=0; i<N_hist_num+1; i++){
      for(j=0; j<N_hist_num+1; j++){
        fprintf(fp_hist_num,"%3d %3d %13.6e %13.6e %13.6e\n",
          i,j,
          (i*2.0/N_hist_num-1.0)*R_hist_num,
          (j*2.0/N_hist_num-1.0)*R_hist_num,
          HistNum[m][i][j]/(double)N_smp);
      }
      fprintf(fp_hist_num,"\n");
    }

    fprintf(fp_hist_theta,"# seed=%d\n",seed);
    fprintf(fp_hist_theta,"# %13.6e\n",parT[m]);
    for(i=0; i<N_hist_theta+1; i++){
      fprintf(fp_hist_theta,"%3d %13.6e %13.6e\n",
        i,(i*2.0/N_hist_theta-1.0)*R_hist_theta,
        HistTheta[m][i]/(double)N_smp);
    }

    fprintf(fp_hist_cos6theta,"# seed=%d\n",seed);
    fprintf(fp_hist_cos6theta,"# %13.6e\n",parT[m]);
    for(i=0; i<N_hist_cos6theta+1; i++){
      fprintf(fp_hist_cos6theta,"%3d %13.6e %13.6e\n",
        i,(i*2.0/N_hist_cos6theta-1.0)*R_hist_cos6theta,
        HistCos6Theta[m][i]/(double)N_smp);
    }

    fclose(fp_hist_num);
    fclose(fp_hist_theta);
    fclose(fp_hist_cos6theta);
  }
*/

  for(m=0; m<N_t; m++){// loop for each temperature
    fp_phys_nr=fopen(buf_phys_nr[m],"w");
    fp_phys_sr=fopen(buf_phys_sr[m],"w");
    fp_phys_shiftnr=fopen(buf_phys_shiftnr[m],"w");
    fp_phys_nq=fopen(buf_phys_nq[m],"w");
    fp_phys_sq=fopen(buf_phys_sq[m],"w");
    fp_phys_shiftnq=fopen(buf_phys_shiftnq[m],"w");

    fprintf(fp_phys_nr,"# seed=%d\n",seed);
    fprintf(fp_phys_nr,"# %13.6e\n",parT[m]);
    for(i=0; i<size_L+1; i++){
      for(j=0; j<size_L+1; j++){
        fprintf(fp_phys_nr,"%3d %3d %13.6e\n",
          i,j,PhysQuantNr[m][i%size_L][j%size_L]);
      }
      fprintf(fp_phys_nr,"\n");
    }

    fprintf(fp_phys_sr,"# seed=%d\n",seed);
    fprintf(fp_phys_sr,"# %13.6e\n",parT[m]);
    for(i=0; i<size_L+1; i++){
      for(j=0; j<size_L+1; j++){
        fprintf(fp_phys_sr,"%3d %3d %13.6e\n",
          i,j,PhysQuantSr[m][i%size_L][j%size_L]);
      }
      fprintf(fp_phys_sr,"\n");
    }

    fprintf(fp_phys_shiftnr,"# seed=%d\n",seed);
    fprintf(fp_phys_shiftnr,"# %13.6e\n",parT[m]);
    for(i=0; i<size_L+1; i++){
      for(j=0; j<size_L+1; j++){
        fprintf(fp_phys_shiftnr,"%3d %3d %13.6e\n",
          i,j,PhysQuantShiftNr[m][i%size_L][j%size_L]);
      }
      fprintf(fp_phys_shiftnr,"\n");
    }

    fprintf(fp_phys_nq,"# seed=%d\n",seed);
    fprintf(fp_phys_nq,"# %13.6e\n",parT[m]);
    for(i=0; i<size_L+1; i++){
      for(j=0; j<size_L+1; j++){
        fprintf(fp_phys_nq,"%3d %3d %13.6e %13.6e\n",
          i,j,PhysQuantNqRe[m][i%size_L][j%size_L],
          PhysQuantNqIm[m][i%size_L][j%size_L]);
      }
      fprintf(fp_phys_nq,"\n");
    }

    fprintf(fp_phys_sq,"# seed=%d\n",seed);
    fprintf(fp_phys_sq,"# %13.6e\n",parT[m]);
    for(i=0; i<size_L+1; i++){
      for(j=0; j<size_L+1; j++){
        fprintf(fp_phys_sq,"%3d %3d %13.6e %13.6e\n",
          i,j,PhysQuantSqRe[m][i%size_L][j%size_L],
          PhysQuantSqIm[m][i%size_L][j%size_L]);
      }
      fprintf(fp_phys_sq,"\n");
    }

    fprintf(fp_phys_shiftnq,"# seed=%d\n",seed);
    fprintf(fp_phys_shiftnq,"# %13.6e\n",parT[m]);
    for(i=0; i<size_L+1; i++){
      for(j=0; j<size_L+1; j++){
        fprintf(fp_phys_shiftnq,"%3d %3d %13.6e %13.6e\n",
          i,j,PhysQuantShiftNqRe[m][i%size_L][j%size_L],
          PhysQuantShiftNqIm[m][i%size_L][j%size_L]);
      }
      fprintf(fp_phys_shiftnq,"\n");
    }

    fclose(fp_phys_nr);
    fclose(fp_phys_sr);
    fclose(fp_phys_shiftnr);
    fclose(fp_phys_nq);
    fclose(fp_phys_sq);
    fclose(fp_phys_shiftnq);
  }

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
  free(Nadd);
  free(n_up);
  free(n_dn);
  free(BoltzFac);
  free(Egap);
  free(cosqr);
  free(sinqr);
  free(corr_c);
  free(corr_s);
  free(PhysQuantNr);
  free(PhysQuantSr);
  free(PhysQuantShiftNr);
  free(PhysQuantNqRe);
  free(PhysQuantNqIm);
  free(PhysQuantSqRe);
  free(PhysQuantSqIm);
  free(PhysQuantShiftNqRe);
  free(PhysQuantShiftNqIm);
  free(HistNum);
  free(HistTheta);
  free(HistCos6Theta);
  free(buf_hist_num);
  free(buf_hist_theta);
  free(buf_hist_cos6theta);
  free(buf_phys_nr);
  free(buf_phys_sr);
  free(buf_phys_shiftnr);
  free(buf_phys_nq);
  free(buf_phys_sq);
  free(buf_phys_shiftnq);

  return 0;
}
