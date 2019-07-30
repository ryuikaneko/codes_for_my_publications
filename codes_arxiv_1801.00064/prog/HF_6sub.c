#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#define Nspin 2 // spin: up and down
#define N 6 // orbital: a,b,c,d,e,f
#define Ndirec 3 // direction: 1,2,3

double complex A[N*N];
double complex copyA[N*N];

double complex Guij[N*N];
double complex Gdij[N*N];
double complex oldGuij[N*N];
double complex oldGdij[N*N];

double complex parN[Nspin][N];// parN[spin][orbital]
double complex oldparN[Nspin][N];
double complex parNtotal;
double complex parChi[Nspin][N][Ndirec];// parChi[spin][orbital(a,b,c,d,e,f)][direction(1,2,3)]
double complex oldparChi[Nspin][N][Ndirec];

double w[N];
double complex work[3*N];
double rwork[3*N];
double doublekmax=2.0*M_PI;
int Niter=500;
double mixing=0.9;
double orderpara_eps=1.0e-12;
int init_state=0;

//int Nspin=2;// spin up and down
//int Ndirec=3;
int intkmax=24;// L=intkmax, Ns=6*L^2
double filling=0.5;// 1/2 filling

double part=1.0;
double parU=0.0;
double parV2=0.0;

double complex tmp;
double tmp_err;

int zheev_(char *jobz, char *uplo, int *n, double complex *a, int *lda, double *w, double complex *work, int *lwork, double *rwork, int *info);

struct myData {
  double data;
  double complex vec[N];
  double px;
  double py;
  int ix;
  int iy;
  int sigma;
  int orig_pos; // will hold the original position of data on your array
};

int myData_compare (const void* a, const void* b) {
  double xx = ((struct myData*)a)->data;
  double yy = ((struct myData*)b)->data;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}

int main(int argc, char *argv[])
{
  int c;
  int show_band = 0;

  if (argv[optind] == NULL || argv[optind + 1] == NULL) {
    fprintf(stderr,"arguments missing\n");
    fprintf(stderr,"usage: %s -U [U] -V [V2] -k [kmax] -i [Niter] -m [mix] -s [initial state] -b (show DOS)\n",argv[0]);
    fprintf(stderr,"  initial state = 0: uniform metal \n");
    fprintf(stderr,"  initial state = 1: AF \n");
    fprintf(stderr,"  initial state = 2: 001122 CDW (nonmag) \n");
    fprintf(stderr,"  initial state = 3: 001122 CDW (magnetic) \n");
    fprintf(stderr,"  initial state = 4: 220200 CDW \n");
    fprintf(stderr,"  initial state = 5: spin orbit \n");
    exit(1);
  }
  while((c = getopt(argc, argv, "U:P:V:k:i:m:s:b")) != -1){
    switch (c){
      case 'U':
        parU = strtod(optarg,NULL);
        break;
      case 'V':
        parV2 = strtod(optarg,NULL);
        break;
      case 'k':
        intkmax = strtol(optarg,NULL,10);
        break;
      case 'i':
        Niter = strtol(optarg,NULL,10);
        break;
      case 'm':
        mixing = strtod(optarg,NULL);
        break;
      case 's':
        init_state = strtol(optarg,NULL,10);
        break;
      case 'b':
        show_band = 1;
        break;
      default:
        fprintf(stderr,"usage: %s -U [U] -V [V2] -k [kmax] -i [Niter] -m [mix] -s [initial state] -b (show DOS)\n",argv[0]);
        exit(1);
    }
  }

  int i,j,k;
  char jobz = 'V';
  char uplo ='U';
  int n=N,lda=N,lwork=3*N,info;

  int intkx,intky;
  int numk;
  int spin;
  int Norbital=Nspin*N*intkmax*intkmax;// 2*6*L^2 = 12*L^2 = (# eigenvalue)
  int Ne=(int)Norbital*filling;
  int Nexpik=3;
  int Nepsilon=3;
  int Nfourier=intkmax*intkmax;
  double invintkmax=1.0/intkmax;
  double kx,ky;
  double Ene=0.0;
  double Ene_0=0.0;
  double Ene_U=0.0;

//  int Ndiv=intkmax;
  int Ndiv=(int)sqrt(Nspin*N)*intkmax;
  int intdos[Ndiv];
  int intdosu[Ndiv];
  int intdosd[Ndiv];
  double Efermi=0.0;
  double Emin;
  double Emax;
  double Edelta;
  double dos[Ndiv];
  double dosu[Ndiv];
  double dosd[Ndiv];

  struct myData* array = (struct myData*) malloc(Norbital*sizeof(struct myData));

  // Expik
  double complex ****Expik;
  Expik = (double complex****)malloc(N*sizeof(double complex***));
  Expik[0] = (double complex***)malloc(N*Nexpik*sizeof(double complex**));
  Expik[0][0] = (double complex**)malloc(N*Nexpik*intkmax*sizeof(double complex*));
  Expik[0][0][0] = (double complex*)malloc(N*Nexpik*intkmax*intkmax*sizeof(double complex));
  for(i=0; i<N; i++){
    Expik[i] = Expik[0] + i*Nexpik;
    for(j=0; j<Nexpik; j++){
      Expik[i][j] = Expik[0][0] + (i*Nexpik + j)*intkmax;
      for(k=0; k<intkmax; k++){
        Expik[i][j][k] = Expik[0][0][0] + ((i*Nexpik + j)*intkmax + k)*intkmax;
      }
    }
  }
  // epsilon: EB, DA, FC
  double complex ***Epsilon;
  Epsilon = (double complex***)malloc(Nepsilon*sizeof(double complex**));
  Epsilon[0] = (double complex**)malloc(Nepsilon*intkmax*sizeof(double complex*));
  Epsilon[0][0] = (double complex*)malloc(Nepsilon*intkmax*intkmax*sizeof(double complex));
  for(i=0; i<Nepsilon; i++){
    Epsilon[i] = Epsilon[0] + i*intkmax;
    for(j=0; j<intkmax; j++){
      Epsilon[i][j] = Epsilon[0][0] + (i*intkmax + j)*intkmax;
    }
  }

  int step;

  for(i=0; i<N*N; i++){
    oldGuij[i]=0.0;
    oldGdij[i]=0.0;
  }

  // initial condition: normal metal
  for(spin=0; spin<2; spin++){
    for(i=0; i<N; i++){
      parN[spin][i] = filling;
      oldparN[spin][i] = parN[spin][i];
      for(j=0; j<Ndirec; j++){
        parChi[spin][i][j] = 0.0;
        oldparChi[spin][i][j] = parChi[spin][i][j];
      }
    }
  }

  if(init_state==1){
    // AF
    parN[0][0] = 1.0;
    parN[1][0] = 0.0;
    parN[0][1] = 0.0;
    parN[1][1] = 1.0;
    //
    parN[0][2] = 1.0;
    parN[1][2] = 0.0;
    parN[0][3] = 0.0;
    parN[1][3] = 1.0;
    //
    parN[0][4] = 1.0;
    parN[1][4] = 0.0;
    parN[0][5] = 0.0;
    parN[1][5] = 1.0;
    //
    // add real cc term for stable optimization
    parChi[0][0][0] = -0.05;
    parChi[0][0][1] = -0.05;
    parChi[0][0][2] = -0.05;
    //
    parChi[0][1][0] = +0.05;
    parChi[0][1][1] = +0.05;
    parChi[0][1][2] = +0.05;
    //
    parChi[0][2][0] = -0.05;
    parChi[0][2][1] = -0.05;
    parChi[0][2][2] = -0.05;
    //
    parChi[0][3][0] = +0.05;
    parChi[0][3][1] = +0.05;
    parChi[0][3][2] = +0.05;
    //
    parChi[0][4][0] = -0.05;
    parChi[0][4][1] = -0.05;
    parChi[0][4][2] = -0.05;
    //
    parChi[0][5][0] = +0.05;
    parChi[0][5][1] = +0.05;
    parChi[0][5][2] = +0.05;
    //
    parChi[1][0][0] = +0.05;
    parChi[1][0][1] = +0.05;
    parChi[1][0][2] = +0.05;
    //
    parChi[1][1][0] = -0.05;
    parChi[1][1][1] = -0.05;
    parChi[1][1][2] = -0.05;
    //
    parChi[1][2][0] = +0.05;
    parChi[1][2][1] = +0.05;
    parChi[1][2][2] = +0.05;
    //
    parChi[1][3][0] = -0.05;
    parChi[1][3][1] = -0.05;
    parChi[1][3][2] = -0.05;
    //
    parChi[1][4][0] = +0.05;
    parChi[1][4][1] = +0.05;
    parChi[1][4][2] = +0.05;
    //
    parChi[1][5][0] = -0.05;
    parChi[1][5][1] = -0.05;
    parChi[1][5][2] = -0.05;
  }else if(init_state==2){
    // initial condition: 012 charge order (nonmag)
    //   |  up |  dn
    // --------------
    // A | 0.0 | 0.0
    // B | 0.0 | 0.0
    // C | 0.5 | 0.5
    // D | 0.5 | 0.5
    // E | 1.0 | 1.0
    // F | 1.0 | 1.0
    parN[0][0] = 0.0;
    parN[1][0] = 0.0;
    parN[0][1] = 0.0;
    parN[1][1] = 0.0;
    //
    parN[0][2] = 0.5;
    parN[1][2] = 0.5;
    parN[0][3] = 0.5;
    parN[1][3] = 0.5;
    //
    parN[0][4] = 1.0;
    parN[1][4] = 1.0;
    parN[0][5] = 1.0;
    parN[1][5] = 1.0;
  }else if(init_state==3){
    // initial condition: 012 charge order (AF)
    //   |  up |  dn
    // --------------
    // A | 0.0 | 0.0
    // B | 0.0 | 0.0
    // C | 1.0 | 0.0
    // D | 0.0 | 1.0
    // E | 1.0 | 1.0
    // F | 1.0 | 1.0
    parN[0][0] = 0.0;
    parN[1][0] = 0.0;
    parN[0][1] = 0.0;
    parN[1][1] = 0.0;
    //
    parN[0][2] = 1.0;
    parN[1][2] = 0.0;
    parN[0][3] = 0.0;
    parN[1][3] = 1.0;
    //
    parN[0][4] = 1.0;
    parN[1][4] = 1.0;
    parN[0][5] = 1.0;
    parN[1][5] = 1.0;
  }else if(init_state==4){
    // initial condition: 220200 charge order (Kurita state)
    parN[0][0] = 1.0;
    parN[1][0] = 1.0;
    parN[0][1] = 0.75;
    parN[1][1] = 0.75;
    //
    parN[0][2] = 0.25;
    parN[1][2] = 0.25;
    parN[0][3] = 0.75;
    parN[1][3] = 0.75;
    //
    parN[0][4] = 0.25;
    parN[1][4] = 0.25;
    parN[0][5] = 0.0;
    parN[1][5] = 0.0;
  }else if(init_state==5){
    // SO init
    parChi[0][0][0] = I * 1.0;
    parChi[0][0][1] = I * 1.0;
    parChi[0][0][2] = I * 1.0;
    //
    parChi[0][1][0] = I * 1.0;
    parChi[0][1][1] = I * 1.0;
    parChi[0][1][2] = I * 1.0;
    //
    parChi[0][2][0] = I * 1.0;
    parChi[0][2][1] = I * 1.0;
    parChi[0][2][2] = I * 1.0;
    //
    parChi[0][3][0] = I * 1.0;
    parChi[0][3][1] = I * 1.0;
    parChi[0][3][2] = I * 1.0;
    //
    parChi[0][4][0] = I * 1.0;
    parChi[0][4][1] = I * 1.0;
    parChi[0][4][2] = I * 1.0;
    //
    parChi[0][5][0] = I * 1.0;
    parChi[0][5][1] = I * 1.0;
    parChi[0][5][2] = I * 1.0;
    //
    parChi[1][0][0] = -I * 1.0;
    parChi[1][0][1] = -I * 1.0;
    parChi[1][0][2] = -I * 1.0;
    //
    parChi[1][1][0] = -I * 1.0;
    parChi[1][1][1] = -I * 1.0;
    parChi[1][1][2] = -I * 1.0;
    //
    parChi[1][2][0] = -I * 1.0;
    parChi[1][2][1] = -I * 1.0;
    parChi[1][2][2] = -I * 1.0;
    //
    parChi[1][3][0] = -I * 1.0;
    parChi[1][3][1] = -I * 1.0;
    parChi[1][3][2] = -I * 1.0;
    //
    parChi[1][4][0] = -I * 1.0;
    parChi[1][4][1] = -I * 1.0;
    parChi[1][4][2] = -I * 1.0;
    //
    parChi[1][5][0] = -I * 1.0;
    parChi[1][5][1] = -I * 1.0;
    parChi[1][5][2] = -I * 1.0;
  }

  // initialize oldparN[spin][i]
  // initialize oldparChi[spin][i]
  for(spin=0; spin<2; spin++){
    for(i=0; i<N; i++){
      oldparN[spin][i] = parN[spin][i];
      for(j=0; j<Ndirec; j++){
        oldparChi[spin][i][j] = parChi[spin][i][j];
      }
    }
  }

  printf("# filling: %f\n",filling);
  printf("# t: %f\n",part);
  printf("# U: %f\n",parU);
  printf("# V2: %f\n",parV2);
  printf("# intkmax: %d\n",intkmax);
  printf("# mixing: %f\n",mixing);
  printf("# orderpara_eps: %e\n",orderpara_eps);
  printf("# initial state: %d\n",init_state);
  printf("\n\n");
  printf("# ene ene_0 ene_U");
  printf("\n");

  for(intkx=0; intkx<intkmax; intkx++){
    for(intky=0; intky<intkmax; intky++){
      kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
      ky = doublekmax*((intky+0.0)*invintkmax-0.5);
      Epsilon[0][intkx][intky] = -part * cexp(-I*kx);
      Epsilon[1][intkx][intky] = -part * cexp(-I*(kx+ky));
      Epsilon[2][intkx][intky] = -part * cexp(-I*ky);
      Expik[0][0][intkx][intky] = 1.0;
      Expik[0][1][intkx][intky] = cexp(-I*ky);
      Expik[0][2][intkx][intky] = cexp(+I*kx);
      Expik[1][0][intkx][intky] = cexp(+I*kx);
      Expik[1][1][intkx][intky] = cexp(+I*(kx+ky));
      Expik[1][2][intkx][intky] = 1.0;
      Expik[2][0][intkx][intky] = cexp(+I*ky);
      Expik[2][1][intkx][intky] = 1.0;
      Expik[2][2][intkx][intky] = cexp(+I*(kx+ky));
      Expik[3][0][intkx][intky] = conj(Expik[0][0][intkx][intky]);
      Expik[3][1][intkx][intky] = conj(Expik[0][1][intkx][intky]);
      Expik[3][2][intkx][intky] = conj(Expik[0][2][intkx][intky]);
      Expik[4][0][intkx][intky] = conj(Expik[1][0][intkx][intky]);
      Expik[4][1][intkx][intky] = conj(Expik[1][1][intkx][intky]);
      Expik[4][2][intkx][intky] = conj(Expik[1][2][intkx][intky]);
      Expik[5][0][intkx][intky] = conj(Expik[2][0][intkx][intky]);
      Expik[5][1][intkx][intky] = conj(Expik[2][1][intkx][intky]);
      Expik[5][2][intkx][intky] = conj(Expik[2][2][intkx][intky]);
    }
  }

  for(step=0; step<Niter; step++){// start step loop
    numk=0;
    for(spin=0; spin<2; spin++){
      for(intkx=0; intkx<intkmax; intkx++){
        for(intky=0; intky<intkmax; intky++){
          kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
          ky = doublekmax*((intky+0.0)*invintkmax-0.5);
// 6x6 matrix
//----
//  0x  1o  2o  3o  4o  5o
//  6x  7x  8o  9o 10o 11o
// 12x 13x 14x 15o 16o 17o
// 18x 19x 20x 21x 22o 23o
// 24x 25x 26x 27x 28x 29o
// 30x 31x 32x 33x 34x 35x
            for(i=0; i<N*N; i++){
              A[i]=0.0;
              copyA[i]=0.0;
            }

            A[0*N+0] = parU * parN[1-spin][0] + 3.0*parV2 * (parN[0][2]+parN[1][2]+parN[0][4]+parN[1][4]);
            A[1*N+1] = parU * parN[1-spin][1] + 3.0*parV2 * (parN[0][3]+parN[1][3]+parN[0][5]+parN[1][5]);
            A[2*N+2] = parU * parN[1-spin][2] + 3.0*parV2 * (parN[0][4]+parN[1][4]+parN[0][0]+parN[1][0]);
            A[3*N+3] = parU * parN[1-spin][3] + 3.0*parV2 * (parN[0][5]+parN[1][5]+parN[0][1]+parN[1][1]);
            A[4*N+4] = parU * parN[1-spin][4] + 3.0*parV2 * (parN[0][0]+parN[1][0]+parN[0][2]+parN[1][2]);
            A[5*N+5] = parU * parN[1-spin][5] + 3.0*parV2 * (parN[0][1]+parN[1][1]+parN[0][3]+parN[1][3]);

            A[1*N+0] = -part;
            A[2*N+1] = -part;
            A[3*N+2] = -part;
            A[4*N+3] = -part;
            A[5*N+4] = -part;
            A[5*N+0] = -part;

            A[3*N+0] = Epsilon[0][intkx][intky];
            A[4*N+1] = Epsilon[1][intkx][intky];
            A[5*N+2] = Epsilon[2][intkx][intky];

            A[2*N+0] = -parV2 * (
              + parChi[spin][0][0] * conj(Expik[0][0][intkx][intky])
              + parChi[spin][0][1] * conj(Expik[0][1][intkx][intky])
              + parChi[spin][0][2] * conj(Expik[0][2][intkx][intky])
              );
            A[3*N+1] = -parV2 * (
              + parChi[spin][1][0] * conj(Expik[1][0][intkx][intky])
              + parChi[spin][1][1] * conj(Expik[1][1][intkx][intky])
              + parChi[spin][1][2] * conj(Expik[1][2][intkx][intky])
              );
            A[4*N+2] = -parV2 * (
              + parChi[spin][2][0] * conj(Expik[2][0][intkx][intky])
              + parChi[spin][2][1] * conj(Expik[2][1][intkx][intky])
              + parChi[spin][2][2] * conj(Expik[2][2][intkx][intky])
              );
            A[5*N+3] = -parV2 * (
              + parChi[spin][3][0] * conj(Expik[3][0][intkx][intky])
              + parChi[spin][3][1] * conj(Expik[3][1][intkx][intky])
              + parChi[spin][3][2] * conj(Expik[3][2][intkx][intky])
              );
            A[4*N+0] = -parV2 * conj(
              + parChi[spin][4][0] * conj(Expik[4][0][intkx][intky])
              + parChi[spin][4][1] * conj(Expik[4][1][intkx][intky])
              + parChi[spin][4][2] * conj(Expik[4][2][intkx][intky])
              );
            A[5*N+1] = -parV2 * conj(
              + parChi[spin][5][0] * conj(Expik[5][0][intkx][intky])
              + parChi[spin][5][1] * conj(Expik[5][1][intkx][intky])
              + parChi[spin][5][2] * conj(Expik[5][2][intkx][intky])
              );

            for(i=0; i<N; i++){
              for(j=0; j<i; j++){
                A[i+j*N] = conj(A[i*N+j]);
              }
            }

            for(i=0; i<N; i++){
              for(j=0; j<N; j++){
                copyA[i+j*N]=A[i*N+j];
              }
            }

            zheev_(&jobz, &uplo, &n, copyA, &lda, w, work, &lwork, rwork, &info);

            for(i=0; i<N; i++){
              array[i + numk*N].data = w[i];
              array[i + numk*N].orig_pos = i + numk*N;
              array[i + numk*N].sigma = spin;
              array[i + numk*N].px = kx;
              array[i + numk*N].py = ky;
              array[i + numk*N].ix = intkx;
              array[i + numk*N].iy = intky;
              for(j=0; j<N; j++){
                array[i + numk*N].vec[j] = copyA[i*N+j];
              }
            }
            numk++;
        }
      }
    }

    qsort(array, Norbital, sizeof(struct myData), &myData_compare);

    Ene_0=0.0;
    for(i=0; i<Ne; i++){
      Ene_0 += array[i].data;
    }

    Ene_U=0.0;
    Ene_U +=
      -parU*Nfourier*(
      + conj(parN[0][0]) * parN[1][0]
      + conj(parN[0][1]) * parN[1][1]
      + conj(parN[0][2]) * parN[1][2]
      + conj(parN[0][3]) * parN[1][3]
      + conj(parN[0][4]) * parN[1][4]
      + conj(parN[0][5]) * parN[1][5]
      )
      -3.0*parV2*Nfourier*(
      + conj(parN[0][0]+parN[1][0]) * (parN[0][2]+parN[1][2])
      + conj(parN[0][1]+parN[1][1]) * (parN[0][3]+parN[1][3])
      + conj(parN[0][2]+parN[1][2]) * (parN[0][4]+parN[1][4])
      + conj(parN[0][3]+parN[1][3]) * (parN[0][5]+parN[1][5])
      + conj(parN[0][4]+parN[1][4]) * (parN[0][0]+parN[1][0])
      + conj(parN[0][5]+parN[1][5]) * (parN[0][1]+parN[1][1])
      )
      +parV2*Nfourier*(
      + parChi[0][0][0] * conj(parChi[0][0][0])
      + parChi[0][0][1] * conj(parChi[0][0][1])
      + parChi[0][0][2] * conj(parChi[0][0][2])
      + parChi[0][1][0] * conj(parChi[0][1][0])
      + parChi[0][1][1] * conj(parChi[0][1][1])
      + parChi[0][1][2] * conj(parChi[0][1][2])
      + parChi[0][2][0] * conj(parChi[0][2][0])
      + parChi[0][2][1] * conj(parChi[0][2][1])
      + parChi[0][2][2] * conj(parChi[0][2][2])
      + parChi[0][3][0] * conj(parChi[0][3][0])
      + parChi[0][3][1] * conj(parChi[0][3][1])
      + parChi[0][3][2] * conj(parChi[0][3][2])
      + parChi[0][4][0] * conj(parChi[0][4][0])
      + parChi[0][4][1] * conj(parChi[0][4][1])
      + parChi[0][4][2] * conj(parChi[0][4][2])
      + parChi[0][5][0] * conj(parChi[0][5][0])
      + parChi[0][5][1] * conj(parChi[0][5][1])
      + parChi[0][5][2] * conj(parChi[0][5][2])
      + parChi[1][0][0] * conj(parChi[1][0][0])
      + parChi[1][0][1] * conj(parChi[1][0][1])
      + parChi[1][0][2] * conj(parChi[1][0][2])
      + parChi[1][1][0] * conj(parChi[1][1][0])
      + parChi[1][1][1] * conj(parChi[1][1][1])
      + parChi[1][1][2] * conj(parChi[1][1][2])
      + parChi[1][2][0] * conj(parChi[1][2][0])
      + parChi[1][2][1] * conj(parChi[1][2][1])
      + parChi[1][2][2] * conj(parChi[1][2][2])
      + parChi[1][3][0] * conj(parChi[1][3][0])
      + parChi[1][3][1] * conj(parChi[1][3][1])
      + parChi[1][3][2] * conj(parChi[1][3][2])
      + parChi[1][4][0] * conj(parChi[1][4][0])
      + parChi[1][4][1] * conj(parChi[1][4][1])
      + parChi[1][4][2] * conj(parChi[1][4][2])
      + parChi[1][5][0] * conj(parChi[1][5][0])
      + parChi[1][5][1] * conj(parChi[1][5][1])
      + parChi[1][5][2] * conj(parChi[1][5][2])
      );
    Ene = Ene_0 + Ene_U;

    if(step%10==0){
      parNtotal = 0.0;
      for(spin=0; spin<2; spin++){
        for(i=0; i<N; i++){
          parNtotal += parN[spin][i];
        }
      }
      printf("%6d %13.10f ",step,Ene/Nfourier/N);
      printf("%13.10f %13.10f ",Ene_0/Nfourier/N,Ene_U/Nfourier/N);
      printf("%13.10f ",creal(parNtotal));
      for(spin=0; spin<2; spin++){
        for(i=0; i<N; i++){
          printf("%13.10f ",creal(parN[spin][i]));
        }
      }
      for(spin=0; spin<2; spin++){
        for(i=0; i<N; i++){
          for(j=0; j<Ndirec; j++){
            printf("%13.10f %13.10f ",creal(parChi[spin][i][j]),cimag(parChi[spin][i][j]));
          }
        }
      }
      printf("\n");
    }

    for(i=0; i<N*N; i++){
      Guij[i]=0.0;
      Gdij[i]=0.0;
    }
    for(k=0; k<Ne; k++){
      if(array[k].sigma==0){
        for(i=0; i<N; i++){
          tmp = conj(array[k].vec[i]);
          for(j=0; j<N; j++){
            Guij[i*N+j] += tmp * array[k].vec[j];
          }
        }
      }else{
        for(i=0; i<N; i++){
          tmp = conj(array[k].vec[i]);
          for(j=0; j<N; j++){
            Gdij[i*N+j] += tmp * array[k].vec[j];
          }
        }
      }
    }
    for(i=0; i<N; i++){
      parN[0][i] = Guij[i*N+i]/(1.0*Nfourier);
      parN[1][i] = Gdij[i*N+i]/(1.0*Nfourier);
    }

    for(spin=0; spin<2; spin++){
      for(i=0; i<N; i++){// orbital
        for(j=0; j<Ndirec; j++){// direction
          parChi[spin][i][j] = 0.0;
        }
      }
    }
    for(k=0; k<Ne; k++){
      spin = array[k].sigma;
      intkx = array[k].ix;
      intky = array[k].iy;
      for(i=0; i<N; i++){// orbital
        for(j=0; j<Ndirec; j++){// direction
          parChi[spin][i][j] +=
            conj(array[k].vec[i]) * array[k].vec[(i+2)%N] // < a^dag c >: orbital i -> orbital (i+2)%N
            * Expik[i][j][intkx][intky]; // e^ik.index
        }
      }
    }
    for(spin=0; spin<2; spin++){
      for(i=0; i<N; i++){// orbital
        for(j=0; j<Ndirec; j++){// direction
          parChi[spin][i][j] /= (1.0*Nfourier);
        }
      }
    }

    tmp_err = 0.0;
    for(i=0; i<N*N; i++){
      tmp_err += cabs(oldGuij[i] - Guij[i])*cabs(oldGuij[i] - Guij[i]);
      tmp_err += cabs(oldGdij[i] - Gdij[i])*cabs(oldGdij[i] - Gdij[i]);
      for(spin=0; spin<2; spin++){
        for(j=0; j<N; j++){// orbital
          for(k=0; k<Ndirec; k++){// direction
            tmp_err += cabs(oldparChi[spin][j][k] - parChi[spin][j][k])*cabs(oldparChi[spin][j][k] - parChi[spin][j][k]);
          }
        }
      }
    }

    if(sqrt(tmp_err)/Nfourier < orderpara_eps){
      break;
    }

    for(spin=0; spin<2; spin++){
      for(i=0; i<N; i++){
        parN[spin][i] = parN[spin][i]*mixing + oldparN[spin][i]*(1.0-mixing);
        for(j=0; j<Ndirec; j++){
          parChi[spin][i][j] = parChi[spin][i][j]*mixing + oldparChi[spin][i][j]*(1.0-mixing);
        }
      }
    }

    for(i=0; i<N*N; i++){
      oldGuij[i] = Guij[i];
      oldGdij[i] = Gdij[i];
    }
    for(spin=0; spin<2; spin++){
      for(i=0; i<N; i++){
        oldparN[spin][i] = parN[spin][i];
        for(j=0; j<Ndirec; j++){
          oldparChi[spin][i][j] = parChi[spin][i][j];
        }
      }
    }
  }// end step loop

  if(step==Niter){
    printf("\n\n");
    printf("# not converged within %d steps !!!\n",Niter);
    printf("# ");
  }else{
    printf("\n\n");
    printf("# converged in %d steps\n",step);
  }

  parNtotal = 0.0;
  for(spin=0; spin<2; spin++){
    for(i=0; i<N; i++){
      parNtotal += parN[spin][i];
    }
  }
  printf("%6d %13.10f ",step,Ene/Nfourier/N);
  printf("%13.10f %13.10f ",Ene_0/Nfourier/N,Ene_U/Nfourier/N);
  printf("%13.10f ",creal(parNtotal));
  for(spin=0; spin<2; spin++){
    for(i=0; i<N; i++){
      printf("%13.10f ",creal(parN[spin][i]));
    }
  }
  for(spin=0; spin<2; spin++){
    for(i=0; i<N; i++){
      for(j=0; j<Ndirec; j++){
        printf("%13.10f %13.10f ",creal(parChi[spin][i][j]),cimag(parChi[spin][i][j]));
      }
    }
  }
  printf("\n");

  if(fabs(array[Ne-1].data-array[Ne].data) < 1.0e-16){
    printf("\n\n");
    printf("# !!! OPEN SHELL !!!\n");
  }

  if(show_band==1){
    Efermi = array[Ne-1].data;
    printf("\n\n");
    printf("# 3: dos\n");
    printf("# Efermi=%f\n",Efermi);
    printf("# Egap=%20.17f\n",array[Ne].data-Efermi);
    Emin=array[0].data-(array[Norbital-1].data-array[0].data)/Ndiv*2.0;
    Emax=array[Norbital-1].data+(array[Norbital-1].data-array[0].data)/Ndiv*2.0;
    Edelta=(Emax-Emin)/Ndiv;
    for(i=0; i<Ndiv; i++){
      intdos[i]=0;
      intdosu[i]=0;
      intdosd[i]=0;
    }
    for(i=0; i<Norbital; i++){
      j=(int)((array[i].data-Emin)/Edelta);
      intdos[j]++;
      if(array[i].sigma==0){
        j=(int)((array[i].data-Emin)/Edelta);
        intdosu[j]++;
      }else{
        j=(int)((array[i].data-Emin)/Edelta);
        intdosd[j]++;
      }
    }
    printf("\n");
    for(i=0; i<Ndiv; i++){
      dos[i]=intdos[i]/(double)(Norbital*Edelta);
      dosu[i]=intdosu[i]/(double)(Norbital*Edelta);
      dosd[i]=intdosd[i]/(double)(Norbital*Edelta);
      printf("%d %f %d %f %d %f %d %f\n",i,Emin+i*Edelta-Efermi,intdos[i],dos[i],intdosu[i],dosu[i],intdosd[i],dosd[i]);
    }
  }

  free(array);
  free(Expik[0][0][0]);
  free(Expik[0][0]);
  free(Expik[0]);
  free(Expik);
  free(Epsilon[0][0]);
  free(Epsilon[0]);
  free(Epsilon);

  return 0;
}
