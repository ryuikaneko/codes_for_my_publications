#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#define N 2

double complex A[N*N];
double complex copyA[N*N];

double complex Guij[N*N];
double complex Gdij[N*N];
double complex oldGuij[N*N];
double complex oldGdij[N*N];

double complex GuijExpPIKx[N*N];
double complex GuijExpPIKy[N*N];
double complex GuijExpPIKxMIKy[N*N];
double complex GuijExpMIKx[N*N];
double complex GuijExpMIKy[N*N];
double complex GuijExpMIKxPIKy[N*N];
double complex GdijExpPIKx[N*N];
double complex GdijExpPIKy[N*N];
double complex GdijExpPIKxMIKy[N*N];
double complex GdijExpMIKx[N*N];
double complex GdijExpMIKy[N*N];
double complex GdijExpMIKxPIKy[N*N];
double complex oldGuijExpPIKx[N*N];
double complex oldGuijExpPIKy[N*N];
double complex oldGuijExpPIKxMIKy[N*N];
double complex oldGuijExpMIKx[N*N];
double complex oldGuijExpMIKy[N*N];
double complex oldGuijExpMIKxPIKy[N*N];
double complex oldGdijExpPIKx[N*N];
double complex oldGdijExpPIKy[N*N];
double complex oldGdijExpPIKxMIKy[N*N];
double complex oldGdijExpMIKx[N*N];
double complex oldGdijExpMIKy[N*N];
double complex oldGdijExpMIKxPIKy[N*N];

double complex parN[2][N];
double complex oldparN[2][N];
double complex parNtotal;
double complex parChi[2][N][3];// parChi[spin][orbital(c,d)][direction(i->i+x,i+x->i+y,i+y->i)]
double complex oldparChi[2][N][3];

double w[N];
double complex work[3*N];
double rwork[3*N];
double doublekmax=2.0*M_PI;
int Niter=500;
double mixing=0.9;
double orderpara_eps=1.0e-12;
int init_state=0;

int Nspin=2;// spin up and down
int Ndirec=3;// i->i+x, i+x->i+y, i+y->i
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
    fprintf(stderr,"  initial state = 2: QAH \n");
    fprintf(stderr,"  initial state = 3: QSH \n");
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
  int Norbital=Nspin*N*intkmax*intkmax;// 2*2*L^2 = 4*L^2 = (# eigenvalue)
  int Ne=(int)Norbital*filling;
  int Nexpik=6;
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

  // Expik: ix, -ix, iy, -iy, -ix+iy, ix-iy
  double complex ***Expik;
  Expik = (double complex***)malloc(Nexpik*sizeof(double complex**));
  Expik[0] = (double complex**)malloc(Nexpik*intkmax*sizeof(double complex*));
  Expik[0][0] = (double complex*)malloc(Nexpik*intkmax*intkmax*sizeof(double complex));
  for(i=0; i<Nexpik; i++){
    Expik[i] = Expik[0] + i*intkmax;
    for(j=0; j<intkmax; j++){
      Expik[i][j] = Expik[0][0] + (i*intkmax + j)*intkmax;
    }
  }
  // E: up, down
  double complex ***Eks;
  Eks = (double complex***)malloc(Nspin*sizeof(double complex**));
  Eks[0] = (double complex**)malloc(Nspin*intkmax*sizeof(double complex*));
  Eks[0][0] = (double complex*)malloc(Nspin*intkmax*intkmax*sizeof(double complex));
  for(i=0; i<Nspin; i++){
    Eks[i] = Eks[0] + i*intkmax;
    for(j=0; j<intkmax; j++){
      Eks[i][j] = Eks[0][0] + (i*intkmax + j)*intkmax;
    }
  }

  int step;

  for(i=0; i<N*N; i++){
    oldGuij[i]=0.0;
    oldGdij[i]=0.0;
    oldGuijExpPIKx[i]=0.0;
    oldGuijExpPIKy[i]=0.0;
    oldGuijExpPIKxMIKy[i]=0.0;
    oldGuijExpMIKx[i]=0.0;
    oldGuijExpMIKy[i]=0.0;
    oldGuijExpMIKxPIKy[i]=0.0;
    oldGdijExpPIKx[i]=0.0;
    oldGdijExpPIKy[i]=0.0;
    oldGdijExpPIKxMIKy[i]=0.0;
    oldGdijExpMIKx[i]=0.0;
    oldGdijExpMIKy[i]=0.0;
    oldGdijExpMIKxPIKy[i]=0.0;
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
    parN[0][0] = 0.5;
    parN[1][0] = 0.0;
    parN[0][1] = 0.0;
    parN[1][1] = 0.5;
  }else if(init_state==2){
    // QAH
    parChi[0][0][0] = I * 1.0;
    parChi[0][0][1] = I * (-1.0);
    parChi[0][0][2] = I * 1.0;
    //
    parChi[0][1][0] = I * 1.0;
    parChi[0][1][1] = I * (-1.0);
    parChi[0][1][2] = I * 1.0;
    //
    parChi[1][0][0] = I * 1.0;
    parChi[1][0][1] = I * (-1.0);
    parChi[1][0][2] = I * 1.0;
    //
    parChi[1][1][0] = I * 1.0;
    parChi[1][1][1] = I * (-1.0);
    parChi[1][1][2] = I * 1.0;
  }else if(init_state==3){
    // QSH
    parChi[0][0][0] = I * 1.0;
    parChi[0][0][1] = I * (-1.0);
    parChi[0][0][2] = I * 1.0;
    //
    parChi[0][1][0] = I * 1.0;
    parChi[0][1][1] = I * (-1.0);
    parChi[0][1][2] = I * 1.0;
    //
    parChi[1][0][0] = -I * 1.0;
    parChi[1][0][1] = -I * (-1.0);
    parChi[1][0][2] = -I * 1.0;
    //
    parChi[1][1][0] = -I * 1.0;
    parChi[1][1][1] = -I * (-1.0);
    parChi[1][1][2] = -I * 1.0;
  }

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
      Eks[0][intkx][intky] = -part * (1.0 + cexp(-I*kx) + cexp(-I*ky));
      Eks[1][intkx][intky] = -part * (1.0 + cexp(-I*kx) + cexp(-I*ky));
      Expik[0][intkx][intky] = cexp(+I*kx);
      Expik[1][intkx][intky] = cexp(-I*kx);
      Expik[2][intkx][intky] = cexp(+I*ky);
      Expik[3][intkx][intky] = cexp(-I*ky);
      Expik[4][intkx][intky] = cexp(-I*(kx-ky));
      Expik[5][intkx][intky] = cexp(+I*(kx-ky));
    }
  }

  for(step=0; step<Niter; step++){// start step loop
    numk=0;
    for(spin=0; spin<2; spin++){
      for(intkx=0; intkx<intkmax; intkx++){
        for(intky=0; intky<intkmax; intky++){
          kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
          ky = doublekmax*((intky+0.0)*invintkmax-0.5);
// 2x2 matrix
//----
//  0x  1o
//  2x  3x
            for(i=0; i<N*N; i++){
              A[i]=0.0;
              copyA[i]=0.0;
            }

            A[0] =
              + parU * parN[1-spin][0]
              - parV2 * (
              - 3.0
              +      parChi[spin][0][0]  * Expik[0][intkx][intky]
              + conj(parChi[spin][0][0]) * Expik[1][intkx][intky]
              +      parChi[spin][0][1]  * Expik[2][intkx][intky]
              + conj(parChi[spin][0][1]) * Expik[3][intkx][intky]
              +      parChi[spin][0][2]  * Expik[4][intkx][intky]
              + conj(parChi[spin][0][2]) * Expik[5][intkx][intky]
              );
            A[1] = Eks[spin][intkx][intky];
            A[2] = conj(Eks[spin][intkx][intky]);
            A[3] =
              + parU * parN[1-spin][1]
              - parV2 * (
              - 3.0
              +      parChi[spin][1][0]  * Expik[1][intkx][intky]
              + conj(parChi[spin][1][0]) * Expik[0][intkx][intky]
              +      parChi[spin][1][1]  * Expik[3][intkx][intky]
              + conj(parChi[spin][1][1]) * Expik[2][intkx][intky]
              +      parChi[spin][1][2]  * Expik[5][intkx][intky]
              + conj(parChi[spin][1][2]) * Expik[4][intkx][intky]
              );

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
      - parU*Nfourier*(
      + parN[0][0] * conj(parN[1][0])
      + parN[0][1] * conj(parN[1][1])
      )
      + parV2*Nfourier*(
      + parChi[0][0][0] * conj(parChi[0][0][0])
      + parChi[0][0][1] * conj(parChi[0][0][1])
      + parChi[0][0][2] * conj(parChi[0][0][2])
      + parChi[0][1][0] * conj(parChi[0][1][0])
      + parChi[0][1][1] * conj(parChi[0][1][1])
      + parChi[0][1][2] * conj(parChi[0][1][2])
      + parChi[1][0][0] * conj(parChi[1][0][0])
      + parChi[1][0][1] * conj(parChi[1][0][1])
      + parChi[1][0][2] * conj(parChi[1][0][2])
      + parChi[1][1][0] * conj(parChi[1][1][0])
      + parChi[1][1][1] * conj(parChi[1][1][1])
      + parChi[1][1][2] * conj(parChi[1][1][2])
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

    for(i=0; i<N*N; i++){
      GuijExpPIKx[i]=0.0;
      GuijExpPIKy[i]=0.0;
      GuijExpPIKxMIKy[i]=0.0;
      GuijExpMIKx[i]=0.0;
      GuijExpMIKy[i]=0.0;
      GuijExpMIKxPIKy[i]=0.0;
      GdijExpPIKx[i]=0.0;
      GdijExpPIKy[i]=0.0;
      GdijExpPIKxMIKy[i]=0.0;
      GdijExpMIKx[i]=0.0;
      GdijExpMIKy[i]=0.0;
      GdijExpMIKxPIKy[i]=0.0;
    }
    for(k=0; k<Ne; k++){
      if(array[k].sigma==0){
        for(i=0; i<N; i++){
          tmp = conj(array[k].vec[i]) * cexp(+I*array[k].px);
          for(j=0; j<N; j++){
            GuijExpPIKx[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(+I*array[k].py);
          for(j=0; j<N; j++){
            GuijExpPIKy[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(+I*array[k].px-I*array[k].py);
          for(j=0; j<N; j++){
            GuijExpPIKxMIKy[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(-I*array[k].px);
          for(j=0; j<N; j++){
            GuijExpMIKx[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(-I*array[k].py);
          for(j=0; j<N; j++){
            GuijExpMIKy[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(-I*array[k].px+I*array[k].py);
          for(j=0; j<N; j++){
            GuijExpMIKxPIKy[i*N+j] += tmp * array[k].vec[j];
          }
        }
      }else{
        for(i=0; i<N; i++){
          tmp = conj(array[k].vec[i]) * cexp(+I*array[k].px);
          for(j=0; j<N; j++){
            GdijExpPIKx[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(+I*array[k].py);
          for(j=0; j<N; j++){
            GdijExpPIKy[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(+I*array[k].px-I*array[k].py);
          for(j=0; j<N; j++){
            GdijExpPIKxMIKy[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(-I*array[k].px);
          for(j=0; j<N; j++){
            GdijExpMIKx[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(-I*array[k].py);
          for(j=0; j<N; j++){
            GdijExpMIKy[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(-I*array[k].px+I*array[k].py);
          for(j=0; j<N; j++){
            GdijExpMIKxPIKy[i*N+j] += tmp * array[k].vec[j];
          }
        }
      }
    }

    for(i=0; i<N; i++){
      parN[0][i] = Guij[i*N+i]/(1.0*Nfourier);
      parN[1][i] = Gdij[i*N+i]/(1.0*Nfourier);
    }

    // spin, orbital, direction
    parChi[0][0][0] = GuijExpMIKx[0*N+0]/(1.0*Nfourier);
    parChi[0][0][1] = GuijExpMIKy[0*N+0]/(1.0*Nfourier);
    parChi[0][0][2] = GuijExpPIKxMIKy[0*N+0]/(1.0*Nfourier);
    parChi[0][1][0] = GuijExpPIKx[1*N+1]/(1.0*Nfourier);
    parChi[0][1][1] = GuijExpPIKy[1*N+1]/(1.0*Nfourier);
    parChi[0][1][2] = GuijExpMIKxPIKy[1*N+1]/(1.0*Nfourier);
    parChi[1][0][0] = GdijExpMIKx[0*N+0]/(1.0*Nfourier);
    parChi[1][0][1] = GdijExpMIKy[0*N+0]/(1.0*Nfourier);
    parChi[1][0][2] = GdijExpPIKxMIKy[0*N+0]/(1.0*Nfourier);
    parChi[1][1][0] = GdijExpPIKx[1*N+1]/(1.0*Nfourier);
    parChi[1][1][1] = GdijExpPIKy[1*N+1]/(1.0*Nfourier);
    parChi[1][1][2] = GdijExpMIKxPIKy[1*N+1]/(1.0*Nfourier);

    tmp_err = 0.0;
    for(i=0; i<N*N; i++){
      tmp_err += cabs(oldGuij[i] - Guij[i])*cabs(oldGuij[i] - Guij[i]);
      tmp_err += cabs(oldGdij[i] - Gdij[i])*cabs(oldGdij[i] - Gdij[i]);
      tmp_err += cabs(oldGuijExpPIKx[i] - GuijExpPIKx[i]) * cabs(oldGuijExpPIKx[i] - GuijExpPIKx[i]);
      tmp_err += cabs(oldGuijExpPIKy[i] - GuijExpPIKy[i]) * cabs(oldGuijExpPIKy[i] - GuijExpPIKy[i]);
      tmp_err += cabs(oldGuijExpPIKxMIKy[i] - GuijExpPIKxMIKy[i]) * cabs(oldGuijExpPIKxMIKy[i] - GuijExpPIKxMIKy[i]);
      tmp_err += cabs(oldGuijExpMIKx[i] - GuijExpMIKx[i]) * cabs(oldGuijExpMIKx[i] - GuijExpMIKx[i]);
      tmp_err += cabs(oldGuijExpMIKy[i] - GuijExpMIKy[i]) * cabs(oldGuijExpMIKy[i] - GuijExpMIKy[i]);
      tmp_err += cabs(oldGuijExpMIKxPIKy[i] - GuijExpMIKxPIKy[i]) * cabs(oldGuijExpMIKxPIKy[i] - GuijExpMIKxPIKy[i]);
      tmp_err += cabs(oldGdijExpPIKx[i] - GdijExpPIKx[i]) * cabs(oldGdijExpPIKx[i] - GdijExpPIKx[i]);
      tmp_err += cabs(oldGdijExpPIKy[i] - GdijExpPIKy[i]) * cabs(oldGdijExpPIKy[i] - GdijExpPIKy[i]);
      tmp_err += cabs(oldGdijExpPIKxMIKy[i] - GdijExpPIKxMIKy[i]) * cabs(oldGdijExpPIKxMIKy[i] - GdijExpPIKxMIKy[i]);
      tmp_err += cabs(oldGdijExpMIKx[i] - GdijExpMIKx[i]) * cabs(oldGdijExpMIKx[i] - GdijExpMIKx[i]);
      tmp_err += cabs(oldGdijExpMIKy[i] - GdijExpMIKy[i]) * cabs(oldGdijExpMIKy[i] - GdijExpMIKy[i]);
      tmp_err += cabs(oldGdijExpMIKxPIKy[i] - GdijExpMIKxPIKy[i]) * cabs(oldGdijExpMIKxPIKy[i] - GdijExpMIKxPIKy[i]);
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
      oldGuijExpPIKx[i] = GuijExpPIKx[i];
      oldGuijExpPIKy[i] = GuijExpPIKy[i];
      oldGuijExpPIKxMIKy[i] = GuijExpPIKxMIKy[i];
      oldGuijExpMIKx[i] = GuijExpMIKx[i];
      oldGuijExpMIKy[i] = GuijExpMIKy[i];
      oldGuijExpMIKxPIKy[i] = GuijExpMIKxPIKy[i];
      oldGdijExpPIKx[i] = GdijExpPIKx[i];
      oldGdijExpPIKy[i] = GdijExpPIKy[i];
      oldGdijExpPIKxMIKy[i] = GdijExpPIKxMIKy[i];
      oldGdijExpMIKx[i] = GdijExpMIKx[i];
      oldGdijExpMIKy[i] = GdijExpMIKy[i];
      oldGdijExpMIKxPIKy[i] = GdijExpMIKxPIKy[i];
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
    printf("# dos\n");
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
  free(Expik[0][0]);
  free(Expik[0]);
  free(Expik);
  free(Eks[0][0]);
  free(Eks[0]);
  free(Eks);

  return 0;
}
