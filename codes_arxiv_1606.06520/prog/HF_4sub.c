#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#define N 4

double complex A[N*N];
double complex copyA[N*N];

double complex Guij[N*N];
double complex Gdij[N*N];

double complex oldGuij[N*N];
double complex oldGdij[N*N];

double complex GuijExpPIKx[N*N];
double complex GuijExpPIKy[N*N];
double complex GuijExpMIKx[N*N];
double complex GuijExpMIKy[N*N];
double complex GdijExpPIKx[N*N];
double complex GdijExpPIKy[N*N];
double complex GdijExpMIKx[N*N];
double complex GdijExpMIKy[N*N];

double complex Oc[2];
double complex Of[2];
double complex Deltac[2];
double complex Deltaf[2];
double complex EksAdd[2];
double complex YksAdd[2];

double complex Nc[2];
double complex Nf[2];
double complex DeltaNc[2];
double complex DeltaNf[2];
double complex Chi[2];
double complex DeltaChi[2];
double complex Eta[2];
double complex DeltaEta[2];

double complex oldNc[2];
double complex oldNf[2];
double complex oldDeltaNc[2];
double complex oldDeltaNf[2];
double complex oldChi[2];
double complex oldDeltaChi[2];
double complex oldEta[2];
double complex oldDeltaEta[2];

double w[N];
double complex work[3*N];
double rwork[3*N];
double doublekmax=2.0*M_PI;
int Niter=500;
double mixing=0.9;
double orderpara_eps=1.0e-12;
int init_state=0;

double Qx=M_PI;
double Qy=M_PI;
int RBZ=2;// reduced BZ
int Nspin=2;// spin up and down
int intkmax=120;// L=intkmax, Ns=2*L^2
double filling=0.75;// 3/4 filling

double tb1=-1.0;// t1
double tb2=0.0;
double tp=0.0;
double tq=-1.0;// t2

/*
// organic charge transfer salts
double tb1=-1.0;
double tb2=-0.5;
double tp=0.0;
double tq=-0.4;
*/

double U=0.0;
double Up=0.0;
double V=0.0;

double complex tmp;
double tmp_err;

int zheev_(char *jobz, char *uplo, int *n, double complex *a, int *lda, double *w, double complex *work, int *lwork, double *rwork, int *info);

struct myData {
  double data;
  double complex vec[N];
  double px;
  double py;
  int sigma;
  int orig_pos;
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
    fprintf(stderr,"usage: %s -U [U] -P [Up] -V [V] -k [kmax] -i [Niter] -m [mix] -s [initial state] -b (show band)\n",argv[0]);
    fprintf(stderr,"  initial state = 0: normal metal \n");
    fprintf(stderr,"  initial state = 1: CO metal \n");
    fprintf(stderr,"  initial state = 2: AF insulator (up-up-down-down magnetic order) \n");
    fprintf(stderr,"  initial state = 3: FM metal (up-up-up-up magnetic order) \n");
    fprintf(stderr,"  initial state = 4: CO+AF insulator (UP-down-DOWN-up magnetic order) \n");
    fprintf(stderr,"  initial state = 5: CO+FM insulator (UP-up-UP-up magnetic order) \n");
    exit(1);
  }
  while((c = getopt(argc, argv, "U:P:V:k:i:m:s:b")) != -1){
    switch (c){
      case 'U':
        U = strtod(optarg,NULL);
        break;
      case 'P':
        Up = strtod(optarg,NULL);
        break;
      case 'V':
        V = strtod(optarg,NULL);
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
        fprintf(stderr,"usage: %s -U [U] -P [Up] -V [V] -k [kmax] -i [Niter] -m [mix] -s [initial state] -b (show band)\n",argv[0]);
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
  int Norbital=Nspin*N*intkmax*intkmax/RBZ;// 2*4*L^2/2 = (2*L)^2 = (# eigenvalue)
  int Ne=(int)Norbital*filling;
  int Nfourier=intkmax*intkmax;
  double invintkmax=1.0/intkmax;
  double kx,ky;
  double Ene=0.0;
  double Ene_0=0.0;
  double Ene_U=0.0;

  int Ndiv=intkmax;
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

  double complex **Epsilon;
  Epsilon = (double complex**)malloc((intkmax)*sizeof(double complex*));
  Epsilon[0] = (double complex*)malloc((intkmax)*(intkmax)*sizeof(double complex));
  for(i=0; i<(intkmax); i++){
    Epsilon[i] = Epsilon[0] + i*(intkmax);
  }
  double complex **EpsilonQ;
  EpsilonQ = (double complex**)malloc((intkmax)*sizeof(double complex*));
  EpsilonQ[0] = (double complex*)malloc((intkmax)*(intkmax)*sizeof(double complex));
  for(i=0; i<(intkmax); i++){
    EpsilonQ[i] = EpsilonQ[0] + i*(intkmax);
  }
  double complex ***Eks;
  Eks = (double complex***)malloc((Nspin)*sizeof(double complex**));
  Eks[0] = (double complex**)malloc((Nspin)*(intkmax)*sizeof(double complex*));
  Eks[0][0] = (double complex*)malloc((Nspin)*(intkmax)*(intkmax)*sizeof(double complex));
  for(i=0; i<(Nspin); i++){
    Eks[i] = Eks[0] + i*(intkmax);
    for(j=0; j<(intkmax); j++){
      Eks[i][j] = Eks[0][0] + (i*(intkmax) + j)*(intkmax);
    }
  }
  double complex ***EQks;
  EQks = (double complex***)malloc((Nspin)*sizeof(double complex**));
  EQks[0] = (double complex**)malloc((Nspin)*(intkmax)*sizeof(double complex*));
  EQks[0][0] = (double complex*)malloc((Nspin)*(intkmax)*(intkmax)*sizeof(double complex));
  for(i=0; i<(Nspin); i++){
    EQks[i] = EQks[0] + i*(intkmax);
    for(j=0; j<(intkmax); j++){
      EQks[i][j] = EQks[0][0] + (i*(intkmax) + j)*(intkmax);
    }
  }
  double complex ***Yks;
  Yks = (double complex***)malloc((Nspin)*sizeof(double complex**));
  Yks[0] = (double complex**)malloc((Nspin)*(intkmax)*sizeof(double complex*));
  Yks[0][0] = (double complex*)malloc((Nspin)*(intkmax)*(intkmax)*sizeof(double complex));
  for(i=0; i<(Nspin); i++){
    Yks[i] = Yks[0] + i*(intkmax);
    for(j=0; j<(intkmax); j++){
      Yks[i][j] = Yks[0][0] + (i*(intkmax) + j)*(intkmax);
    }
  }
  double complex ***YQks;
  YQks = (double complex***)malloc((Nspin)*sizeof(double complex**));
  YQks[0] = (double complex**)malloc((Nspin)*(intkmax)*sizeof(double complex*));
  YQks[0][0] = (double complex*)malloc((Nspin)*(intkmax)*(intkmax)*sizeof(double complex));
  for(i=0; i<(Nspin); i++){
    YQks[i] = YQks[0] + i*(intkmax);
    for(j=0; j<(intkmax); j++){
      YQks[i][j] = YQks[0][0] + (i*(intkmax) + j)*(intkmax);
    }
  }

  int step;

  for(i=0; i<N*N; i++){
    oldGuij[i]=0.0;
    oldGdij[i]=0.0;
  }

  if(init_state==0){
    // normal metal
    Nc[0] = 0.75;
    Nc[1] = 0.75;
    Nf[0] = 0.75;
    Nf[1] = 0.75;
    DeltaNc[0] = 0.0;
    DeltaNc[1] = 0.0;
    DeltaNf[0] = 0.0;
    DeltaNf[1] = 0.0;
    Chi[0] = 0.0;
    Chi[1] = 0.0;
    DeltaChi[0] = 0.0;
    DeltaChi[1] = 0.0;
    Eta[0] = 0.0;
    Eta[1] = 0.0;
    DeltaEta[0] = 0.0;
    DeltaEta[1] = 0.0;
  }else if(init_state==1){
    // CO metal
    Nc[0] = 1.0;
    Nc[1] = 1.0;
    Nf[0] = 0.5;
    Nf[1] = 0.5;
    DeltaNc[0] = 0.0;
    DeltaNc[1] = 0.0;
    DeltaNf[0] = 0.0;
    DeltaNf[1] = 0.0;
    Chi[0] = 0.0;
    Chi[1] = 0.0;
    DeltaChi[0] = 0.0;
    DeltaChi[1] = 0.0;
    Eta[0] = 0.0;
    Eta[1] = 0.0;
    DeltaEta[0] = 0.0;
    DeltaEta[1] = 0.0;
  }else if(init_state==2){
    // AF insulator (up-up-down-down magnetic order)
    Nc[0] = 0.75;
    Nc[1] = 0.75;
    Nf[0] = 0.75;
    Nf[1] = 0.75;
    DeltaNc[0] = 0.25;
    DeltaNc[1] = -0.25;
    DeltaNf[0] = 0.25;
    DeltaNf[1] = -0.25;
    Chi[0] = 0.0;
    Chi[1] = 0.0;
    DeltaChi[0] = 0.0;
    DeltaChi[1] = 0.0;
    Eta[0] = 0.0;
    Eta[1] = 0.0;
    DeltaEta[0] = 0.0;
    DeltaEta[1] = 0.0;
  }else if(init_state==3){
    // FM metal (up-up-up-up magnetic order)
    Nc[0] = 1.0;
    Nc[1] = 0.5;
    Nf[0] = 1.0;
    Nf[1] = 0.5;
    DeltaNc[0] = 0.0;
    DeltaNc[1] = 0.0;
    DeltaNf[0] = 0.0;
    DeltaNf[1] = 0.0;
    Chi[0] = 0.0;
    Chi[1] = 0.0;
    DeltaChi[0] = 0.0;
    DeltaChi[1] = 0.0;
    Eta[0] = 0.0;
    Eta[1] = 0.0;
    DeltaEta[0] = 0.0;
    DeltaEta[1] = 0.0;
  }else if(init_state==4){
    // CO+AF insulator (UP-down-DOWN-up magnetic order)
    Nc[0] = 0.5;
    Nc[1] = 0.5;
    Nf[0] = 1.0;
    Nf[1] = 1.0;
    DeltaNc[0] = 0.48;
    DeltaNc[1] = -0.48;
    DeltaNf[0] = -0.02;
    DeltaNf[1] = 0.02;
    Chi[0] = 0.0;
    Chi[1] = 0.0;
    DeltaChi[0] = 0.0;
    DeltaChi[1] = 0.0;
    Eta[0] = 0.0;
    Eta[1] = 0.0;
    DeltaEta[0] = 0.0;
    DeltaEta[1] = 0.0;
  }else if(init_state==5){
    // CO+FM insulator (UP-up-UP-up magnetic order)
    Nc[0] = 0.98;
    Nc[1] = 0.02;
    Nf[0] = 1.02;
    Nf[1] = 0.98;
    DeltaNc[0] = 0.0;
    DeltaNc[1] = 0.0;
    DeltaNf[0] = 0.0;
    DeltaNf[1] = 0.0;
    Chi[0] = 0.0;
    Chi[1] = 0.0;
    DeltaChi[0] = 0.0;
    DeltaChi[1] = 0.0;
    Eta[0] = 0.0;
    Eta[1] = 0.0;
    DeltaEta[0] = 0.0;
    DeltaEta[1] = 0.0;
  }

  oldNc[0] = 0.0;
  oldNc[1] = 0.0;
  oldNf[0] = 0.0;
  oldNf[1] = 0.0;
  oldDeltaNc[0] = 0.0;
  oldDeltaNc[1] = 0.0;
  oldDeltaNf[0] = 0.0;
  oldDeltaNf[1] = 0.0;
  oldChi[0] = 0.0;
  oldChi[1] = 0.0;
  oldDeltaChi[0] = 0.0;
  oldDeltaChi[1] = 0.0;
  oldEta[0] = 0.0;
  oldEta[1] = 0.0;
  oldDeltaEta[0] = 0.0;
  oldDeltaEta[1] = 0.0;

  printf("# filling: %f\n",filling);
  printf("# tb1: %f\n",tb1);
  printf("# tb2: %f\n",tb2);
  printf("# tp: %f\n",tp);
  printf("# tq: %f\n",tq);
  printf("# U: %f\n",U);
  printf("# Up: %f\n",Up);
  printf("# V: %f\n",V);
  printf("# intkmax: %d\n",intkmax);
  printf("# mixing: %f\n",mixing);
  printf("# orderpara_eps: %e\n",orderpara_eps);
  printf("# initial state: %d\n",init_state);
  printf("\n\n");
  printf("# ene ene_0 ene_U Ave(N) Nc[0] Nc[1] Nf[0] Nf[1] ");
  printf("Ave(DeltaN) DeltaNc[0] DeltaNc[1] DeltaNf[0] DeltaNf[1] ");
  printf("Chi[0] Chi[1] DeltaChi[0] DeltaChi[1]\n");

  for(intkx=0; intkx<intkmax; intkx++){
    for(intky=0; intky<intkmax; intky++){
      kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
      ky = doublekmax*((intky+0.0)*invintkmax-0.5);
      Epsilon[intkx][intky] = tb1+tq*cexp(-I*kx)+tq*cexp(-I*ky)+tb2*cexp(-I*(kx+ky));
      EpsilonQ[intkx][intky] = tb1+tq*cexp(-I*(kx+Qx))+tq*cexp(-I*(ky+Qy))+tb2*cexp(-I*((kx+Qx)+(ky+Qy)));
    }
  }

  for(step=0; step<Niter; step++){// start step loop
    Oc[0] = U*Nc[1] + Up*(Nf[0]+Nf[1]) + 2.0*V*(Nf[0]+Nf[1]);
    Oc[1] = U*Nc[0] + Up*(Nf[0]+Nf[1]) + 2.0*V*(Nf[0]+Nf[1]);
    Of[0] = U*Nf[1] + Up*(Nc[0]+Nc[1]) + 2.0*V*(Nc[0]+Nc[1]);
    Of[1] = U*Nf[0] + Up*(Nc[0]+Nc[1]) + 2.0*V*(Nc[0]+Nc[1]);
    Deltac[0] = U*DeltaNc[1] + Up*(DeltaNf[0]+DeltaNf[1]) + (cexp(-I*Qx)+cexp(-I*Qy))*V*(DeltaNf[0]+DeltaNf[1]);
    Deltac[1] = U*DeltaNc[0] + Up*(DeltaNf[0]+DeltaNf[1]) + (cexp(-I*Qx)+cexp(-I*Qy))*V*(DeltaNf[0]+DeltaNf[1]);
    Deltaf[0] = U*DeltaNf[1] + Up*(DeltaNc[0]+DeltaNc[1]) - 2.0*V*(DeltaNc[0]+DeltaNc[1]);
    Deltaf[1] = U*DeltaNf[0] + Up*(DeltaNc[0]+DeltaNc[1]) - 2.0*V*(DeltaNc[0]+DeltaNc[1]);
    EksAdd[0] = -Up*conj(Chi[0]);
    EksAdd[1] = -Up*conj(Chi[1]);
    YksAdd[0] = -Up*conj(DeltaChi[0]);
    YksAdd[1] = -Up*conj(DeltaChi[1]);

    for(spin=0; spin<2; spin++){
      for(intkx=0; intkx<intkmax; intkx++){
        for(intky=0; intky<intkmax; intky++){
          kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
          ky = doublekmax*((intky+0.0)*invintkmax-0.5);
          Eks[spin][intkx][intky] = Epsilon[intkx][intky] + EksAdd[spin]
            + (-V*conj(Eta[spin])*(cexp(-I*kx)+cexp(-I*ky)) );
          EQks[spin][intkx][intky] = EpsilonQ[intkx][intky] + EksAdd[spin]
            + (-V*conj(Eta[spin])*(cexp(-I*(kx+Qx))+cexp(-I*(ky+Qy))) );
          Yks[spin][intkx][intky] = YksAdd[spin]
            + (-V*conj(DeltaEta[spin])*(cexp(-I*kx)+cexp(-I*ky)) );
          YQks[spin][intkx][intky] = YksAdd[spin]
            + (-V*conj(DeltaEta[spin])*(cexp(-I*(kx+Qx))+cexp(-I*(ky+Qy))) );
        }
      }
    }

    numk=0;
    for(spin=0; spin<2; spin++){
      for(intkx=0; intkx<intkmax; intkx++){
        for(intky=0; intky<intkmax; intky++){
          kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
          ky = doublekmax*((intky+0.0)*invintkmax-0.5);
          if(fabs(kx)+fabs(ky) < doublekmax*0.5+1.0e-10){// reduced BZ for Q=(pi,pi)

            //  0x  1o  2o  3o
            //  4x  5x  6o  7o
            //  8x  9x 10x 11o
            // 12x 13x 14x 15x
            for(i=0; i<N*N; i++){
              A[i]=0.0;
              copyA[i]=0.0;
            }

            A[0]=Oc[spin];
            A[4]=Deltac[spin];
            A[8]=Eks[spin][intkx][intky];
            A[12]=Yks[spin][intkx][intky];
            A[5]=Oc[spin];
            A[9]=YQks[spin][intkx][intky];
            A[13]=EQks[spin][intkx][intky];
            A[10]=Of[spin];
            A[14]=Deltaf[spin];
            A[15]=Of[spin];

            A[1]=conj(A[4]);
            A[2]=conj(A[8]);
            A[3]=conj(A[12]);
            A[6]=conj(A[9]);
            A[7]=conj(A[13]);
            A[11]=conj(A[14]);

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
    }

    qsort(array, Norbital, sizeof(struct myData), &myData_compare);

    Ene_0=0.0;
    for(i=0; i<Ne; i++){
      Ene_0 += array[i].data;
    }

    Ene_U=0.0;
    Ene_U += -U*Nfourier*(
      + conj(Nc[0])*Nc[1] + conj(Nf[0])*Nf[1]
      + conj(DeltaNc[0])*DeltaNc[1] + conj(DeltaNf[0])*DeltaNf[1])
      +Up*Nfourier*(
      - conj(Nc[0])*Nf[0] - conj(Nc[0])*Nf[1]
      - conj(Nc[1])*Nf[0] - conj(Nc[1])*Nf[1]
      - conj(DeltaNc[0])*DeltaNf[0] - conj(DeltaNc[0])*DeltaNf[1]
      - conj(DeltaNc[1])*DeltaNf[0] - conj(DeltaNc[1])*DeltaNf[1]
      + conj(Chi[0])*Chi[0] + conj(Chi[1])*Chi[1]
      + conj(DeltaChi[0])*DeltaChi[0] + conj(DeltaChi[1])*DeltaChi[1])
      +2.0*V*Nfourier*(
      - conj(Nc[0])*Nf[0] - conj(Nc[0])*Nf[1]
      - conj(Nc[1])*Nf[0] - conj(Nc[1])*Nf[1]
      + conj(DeltaNc[0])*DeltaNf[0] + conj(DeltaNc[0])*DeltaNf[1]
      + conj(DeltaNc[1])*DeltaNf[0] + conj(DeltaNc[1])*DeltaNf[1]
      + conj(Eta[0])*Eta[0] + conj(Eta[1])*Eta[1]
      + conj(DeltaEta[0])*DeltaEta[0] + conj(DeltaEta[1])*DeltaEta[1]);
    Ene = Ene_0 + Ene_U;

    if(step%10==0){
      printf("%6d %13.10f ",step,Ene/Nfourier);
      printf("%13.10f %13.10f ",Ene_0/Nfourier,Ene_U/Nfourier);
      printf("%13.10f ",creal(Nc[0]+Nc[1]+Nf[0]+Nf[1])*0.25);
      printf("%13.10f ",creal(Nc[0]));
      printf("%13.10f ",creal(Nc[1]));
      printf("%13.10f ",creal(Nf[0]));
      printf("%13.10f ",creal(Nf[1]));
      printf("%13.10f ",creal(DeltaNc[0]+DeltaNc[1]+DeltaNf[0]+DeltaNf[1])*0.25);
      printf("%13.10f ",creal(DeltaNc[0]));
      printf("%13.10f ",creal(DeltaNc[1]));
      printf("%13.10f ",creal(DeltaNf[0]));
      printf("%13.10f ",creal(DeltaNf[1]));
      printf("%13.10f ",creal(Chi[0]));
      printf("%13.10f ",creal(Chi[1]));
      printf("%13.10f ",creal(DeltaChi[0]));
      printf("%13.10f ",creal(DeltaChi[1]));
      printf("%13.10f ",creal(Eta[0]));
      printf("%13.10f ",creal(Eta[1]));
      printf("%13.10f ",creal(DeltaEta[0]));
      printf("%13.10f ",creal(DeltaEta[1]));
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
      GuijExpMIKx[i]=0.0;
      GuijExpMIKy[i]=0.0;
      GdijExpPIKx[i]=0.0;
      GdijExpPIKy[i]=0.0;
      GdijExpMIKx[i]=0.0;
      GdijExpMIKy[i]=0.0;
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
          tmp = conj(array[k].vec[i]) * cexp(-I*array[k].px);
          for(j=0; j<N; j++){
            GuijExpMIKx[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(-I*array[k].py);
          for(j=0; j<N; j++){
            GuijExpMIKy[i*N+j] += tmp * array[k].vec[j];
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
          tmp = conj(array[k].vec[i]) * cexp(-I*array[k].px);
          for(j=0; j<N; j++){
            GdijExpMIKx[i*N+j] += tmp * array[k].vec[j];
          }
          tmp = conj(array[k].vec[i]) * cexp(-I*array[k].py);
          for(j=0; j<N; j++){
            GdijExpMIKy[i*N+j] += tmp * array[k].vec[j];
          }
        }
      }
    }

    Nc[0] = (Guij[0*N+0]+Guij[1*N+1])/(1.0*Nfourier);
    Nc[1] = (Gdij[0*N+0]+Gdij[1*N+1])/(1.0*Nfourier);
    Nf[0] = (Guij[2*N+2]+Guij[3*N+3])/(1.0*Nfourier);
    Nf[1] = (Gdij[2*N+2]+Gdij[3*N+3])/(1.0*Nfourier);
    DeltaNc[0] = (Guij[0*N+1]+Guij[1*N+0])/(1.0*Nfourier);
    DeltaNc[1] = (Gdij[0*N+1]+Gdij[1*N+0])/(1.0*Nfourier);
    DeltaNf[0] = (Guij[2*N+3]+Guij[3*N+2])/(1.0*Nfourier);
    DeltaNf[1] = (Gdij[2*N+3]+Gdij[3*N+2])/(1.0*Nfourier);
    Chi[0] = (Guij[0*N+2]+Guij[2*N+0]+Guij[1*N+3]+Guij[3*N+1])/(2.0*Nfourier);
    Chi[1] = (Gdij[0*N+2]+Gdij[2*N+0]+Gdij[1*N+3]+Gdij[3*N+1])/(2.0*Nfourier);
    DeltaChi[0] = (Guij[0*N+3]+Guij[3*N+0]+Guij[1*N+2]+Guij[2*N+1])/(2.0*Nfourier);
    DeltaChi[1] = (Gdij[0*N+3]+Gdij[3*N+0]+Gdij[1*N+2]+Gdij[2*N+1])/(2.0*Nfourier);

    Eta[0] = (
      +GuijExpPIKx[0*N+2]+cexp(-I*Qx)*GuijExpPIKx[1*N+3]
      +GuijExpMIKx[2*N+0]+cexp(+I*Qx)*GuijExpMIKx[3*N+1]
      )/(2.0*Nfourier);
    Eta[1] = (
      +GdijExpPIKx[0*N+2]+cexp(-I*Qx)*GdijExpPIKx[1*N+3]
      +GdijExpMIKx[2*N+0]+cexp(+I*Qx)*GdijExpMIKx[3*N+1]
      )/(2.0*Nfourier);
    DeltaEta[0] = (
      +GuijExpPIKx[0*N+3]+cexp(-I*Qx)*GuijExpPIKx[1*N+2]
      +GuijExpMIKx[3*N+0]+cexp(+I*Qx)*GuijExpMIKx[2*N+1]
      )/(2.0*Nfourier);
    DeltaEta[1] = (
      +GdijExpPIKx[0*N+3]+cexp(-I*Qx)*GdijExpPIKx[1*N+2]
      +GdijExpMIKx[3*N+0]+cexp(+I*Qx)*GdijExpMIKx[2*N+1]
      )/(2.0*Nfourier);

    tmp_err = 0.0;
    for(i=0; i<N*N; i++){
      tmp_err += cabs(oldGuij[i] - Guij[i])*cabs(oldGuij[i] - Guij[i]);
      tmp_err += cabs(oldGdij[i] - Gdij[i])*cabs(oldGdij[i] - Gdij[i]);
    }
    if(sqrt(tmp_err)/Nfourier < orderpara_eps
      && fabs(creal(Nc[0] - oldNc[0])) < orderpara_eps
      && fabs(creal(Nc[1] - oldNc[1])) < orderpara_eps
      && fabs(creal(Nf[0] - oldNf[0])) < orderpara_eps
      && fabs(creal(Nf[1] - oldNf[1])) < orderpara_eps
      && fabs(creal(DeltaNc[0] - oldDeltaNc[0])) < orderpara_eps
      && fabs(creal(DeltaNc[1] - oldDeltaNc[1])) < orderpara_eps
      && fabs(creal(DeltaNf[0] - oldDeltaNf[0])) < orderpara_eps
      && fabs(creal(DeltaNf[1] - oldDeltaNf[1])) < orderpara_eps
      && fabs(creal(Chi[0] - oldChi[0])) < orderpara_eps
      && fabs(creal(Chi[1] - oldChi[1])) < orderpara_eps
      && fabs(creal(DeltaChi[0] - oldDeltaChi[0])) < orderpara_eps
      && fabs(creal(DeltaChi[1] - oldDeltaChi[1])) < orderpara_eps
      && fabs(creal(Eta[0] - oldEta[0])) < orderpara_eps
      && fabs(creal(Eta[1] - oldEta[1])) < orderpara_eps
      && fabs(creal(DeltaEta[0] - oldDeltaEta[0])) < orderpara_eps
      && fabs(creal(DeltaEta[1] - oldDeltaEta[1])) < orderpara_eps
    ){
      break;
    }

    for(spin=0; spin<2; spin++){
      Nc[spin] = Nc[spin]*mixing + oldNc[spin]*(1.0-mixing);
      Nf[spin] = Nf[spin]*mixing + oldNf[spin]*(1.0-mixing);
      DeltaNc[spin] = DeltaNc[spin]*mixing + oldDeltaNc[spin]*(1.0-mixing);
      DeltaNf[spin] = DeltaNf[spin]*mixing + oldDeltaNf[spin]*(1.0-mixing);
      Chi[spin] = Chi[spin]*mixing + Chi[spin]*(1.0-mixing);
      DeltaChi[spin] = DeltaChi[spin]*mixing + DeltaChi[spin]*(1.0-mixing);
      Eta[spin] = Eta[spin]*mixing + Eta[spin]*(1.0-mixing);
      DeltaEta[spin] = DeltaEta[spin]*mixing + DeltaEta[spin]*(1.0-mixing);
    }

    for(i=0; i<N*N; i++){
      oldGuij[i] = Guij[i];
      oldGdij[i] = Gdij[i];
    }
    for(spin=0; spin<2; spin++){
      oldNc[spin] = Nc[spin];
      oldNf[spin] = Nf[spin];
      oldDeltaNc[spin] = DeltaNc[spin];
      oldDeltaNf[spin] = DeltaNf[spin];
      oldChi[spin] = Chi[spin];
      oldDeltaChi[spin] = DeltaChi[spin];
      oldEta[spin] = Eta[spin];
      oldDeltaEta[spin] = DeltaEta[spin];
    }

  }// end step loop

  if(step==Niter){
    printf("\n\n");
    printf("# not converged within %d steps !!!\n",Niter);
    printf("# ");
  }else if(
    creal(Nc[0]+DeltaNc[0])>1.0+1.0e10 ||
    creal(Nc[0]+DeltaNc[0])<0.0-1.0e10 ||
    creal(Nc[0]-DeltaNc[0])>1.0+1.0e10 ||
    creal(Nc[0]-DeltaNc[0])<0.0-1.0e10 ||
    creal(Nc[1]+DeltaNc[1])>1.0+1.0e10 ||
    creal(Nc[1]+DeltaNc[1])<0.0-1.0e10 ||
    creal(Nc[1]-DeltaNc[1])>1.0+1.0e10 ||
    creal(Nc[1]-DeltaNc[1])<0.0-1.0e10 ||
    creal(Nf[0]+DeltaNf[0])>1.0+1.0e10 ||
    creal(Nf[0]+DeltaNf[0])<0.0-1.0e10 ||
    creal(Nf[0]-DeltaNf[0])>1.0+1.0e10 ||
    creal(Nf[0]-DeltaNf[0])<0.0-1.0e10 ||
    creal(Nf[1]+DeltaNf[1])>1.0+1.0e10 ||
    creal(Nf[1]+DeltaNf[1])<0.0-1.0e10 ||
    creal(Nf[1]-DeltaNf[1])>1.0+1.0e10 ||
    creal(Nf[1]-DeltaNf[1])<0.0-1.0e10
  ){
    printf("\n\n");
    printf("# converged to <n> < 0 or <n> >1 (%d steps) !!! \n",step);
    printf("# ");
  }else{
    printf("\n\n");
    printf("# converged in %d steps\n",step);
  }
  printf("%6d %13.10f ",step,Ene/Nfourier);
  printf("%13.10f %13.10f ",Ene_0/Nfourier,Ene_U/Nfourier);
  printf("%13.10f ",creal(Nc[0]+Nc[1]+Nf[0]+Nf[1])*0.25);
  printf("%13.10f ",creal(Nc[0]));
  printf("%13.10f ",creal(Nc[1]));
  printf("%13.10f ",creal(Nf[0]));
  printf("%13.10f ",creal(Nf[1]));
  printf("%13.10f ",creal(DeltaNc[0]+DeltaNc[1]+DeltaNf[0]+DeltaNf[1])*0.25);
  printf("%13.10f ",creal(DeltaNc[0]));
  printf("%13.10f ",creal(DeltaNc[1]));
  printf("%13.10f ",creal(DeltaNf[0]));
  printf("%13.10f ",creal(DeltaNf[1]));
  printf("%13.10f ",creal(Chi[0]));
  printf("%13.10f ",creal(Chi[1]));
  printf("%13.10f ",creal(DeltaChi[0]));
  printf("%13.10f ",creal(DeltaChi[1]));
  printf("%13.10f ",creal(Eta[0]));
  printf("%13.10f ",creal(Eta[1]));
  printf("%13.10f ",creal(DeltaEta[0]));
  printf("%13.10f ",creal(DeltaEta[1]));
  printf("\n");

  printf("\n\n");
  printf("# SzcA SzfA SzcB SzfB NcA NfA NcB NfB\n");
  printf("%13.10f ",creal(Nc[0]+DeltaNc[0]-Nc[1]-DeltaNc[1])*0.5);// SzcA
  printf("%13.10f ",creal(Nf[0]+DeltaNf[0]-Nf[1]-DeltaNf[1])*0.5);// SzfA
  printf("%13.10f ",creal(Nc[0]-DeltaNc[0]-Nc[1]+DeltaNc[1])*0.5);// SzcB
  printf("%13.10f ",creal(Nf[0]-DeltaNf[0]-Nf[1]+DeltaNf[1])*0.5);// SzfB
  printf("%13.10f ",creal(Nc[0]+DeltaNc[0]+Nc[1]+DeltaNc[1]));// NcA
  printf("%13.10f ",creal(Nf[0]+DeltaNf[0]+Nf[1]+DeltaNf[1]));// NfA
  printf("%13.10f ",creal(Nc[0]-DeltaNc[0]+Nc[1]-DeltaNc[1]));// NcB
  printf("%13.10f ",creal(Nf[0]-DeltaNf[0]+Nf[1]-DeltaNf[1]));// NfB
  printf("\n");

  if(fabs(array[Ne-1].data-array[Ne].data) < 1.0e-16){
    printf("\n\n");
    printf("# !!! OPEN SHELL !!!\n");
  }

  if(show_band==1){
/*
    printf("# tb1=%f tb2=%f tp=%f tq=%f\n",tb1,tb2,tp,tq);
    printf("\n\n");
    printf("# 1: band\n");
    printf("# number spin intkx intky kx ky band energy orig_posit\n");
    numk=0;
    for(spin=0; spin<2; spin++){
      for(intkx=0; intkx<intkmax; intkx++){
        for(intky=0; intky<intkmax; intky++){
          kx = doublekmax*((intkx+0.5)*invintkmax-0.5);
          ky = doublekmax*((intky+0.0)*invintkmax-0.5);
          if(fabs(kx)+fabs(ky) < doublekmax*0.5+1.0e-10){// reduced BZ for Q=(pi,pi)
            for(i=0; i<N; i++){
              printf("%4d %2d %3d %3d %13.10f %13.10f %2d ",i+numk*N,spin,intkx,intky,kx,ky,i);
              printf("%13.10f %4d ",array[i+numk*N].data,array[i+numk*N].orig_pos);
              printf("[");
              for(j=0; j<N; j++){
                printf("%13.10f+i%13.10f ",creal(array[i+numk*N].vec[j]),cimag(array[i+numk*N].vec[j]));
              }
              printf("]\n");
            }
            numk++;
          }
        }
      }
    }
*/
/*
    printf("\n\n");
    printf("# 2: energy\n");
    printf("# i Eigenvalue orig_posit\n");
    for(i=0; i<Norbital; i++){
      printf("%4d %13.10f %4d ",i,array[i].data,array[i].orig_pos);
      printf("[");
      for(j=0; j<N; j++){
        printf("%13.10f+i%13.10f ",creal(array[i].vec[j]),cimag(array[i].vec[j]));
      }
      printf("]\n");
    }
*/
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
  free(Epsilon[0]);
  free(Epsilon);
  free(EpsilonQ[0]);
  free(EpsilonQ);
  free(Eks[0][0]);
  free(Eks[0]);
  free(Eks);
  free(EQks[0][0]);
  free(EQks[0]);
  free(EQks);
  free(Yks[0][0]);
  free(Yks[0]);
  free(Yks);
  free(YQks[0][0]);
  free(YQks[0]);
  free(YQks);

  return 0;
}
