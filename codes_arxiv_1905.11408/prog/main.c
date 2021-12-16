#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include <complex.h>
#include "./dSFMT/dSFMT.h"
#include "mfmemory.h"
#include "func.h"
#include "lattice.h"

int main(int argc, char *argv[])
{
  char opt;
  double time_count[1];
  int seed;
  int sizeL;// length
  int sizeOrb;// # orbitals
  int sizeCooNum;// # coordination number
  int sizeNs;// system size
  double **spin;
  double **flux;
  int ***listSite;
  int **listLocSite;
  double ****listLocJ;
  double parJx;// Heisenberg
  double parJy;// Heisenberg
  double parJz;// Heisenberg
  double parKx;// Kitaev
  double parKy;// Kitaev
  double parKz;// Kitaev
  double parGx;// Gamma
  double parGy;// Gamma
  double parGz;// Gamma
  double parGxP;// GammaPrime
  double parGyP;// GammaPrime
  double parGzP;// GammaPrime
  double parH[3];// magnetic field
  double parHDir[3];// normalized magnetic field direction
  double parHSize;// size of magnetic field
  double parGK;// Gamma/|Kitaev|
  double ****parJ;
  int i;
  int j;
  int k;
  int site;
  double dist[2];
  int ix;
  int iy;
  int iorb;
  int ixy;
  int ixyorb;
  int step;
  int numStepRnd;
  int numLoopRnd;
  int numStepDet;
  int numLoopDet;
  int option;// init type
  int mode;// rnd or deterministic
  double ene;
  double *locField;
  double parGP;
  double parTmp;
  double parTmpRatio;
  double parInvTmp;
  double spinCand[3];
  double parBF;
  double parTheta;
  double ****listDist;

  time_count[0] = gettimeofday_sec();
  seed = (int)(time_count[0]
    + (time_count[0]-(int)time_count[0])*100000000);
//  seed=time(NULL);
//  seed=1234567;

  sizeL = 12;
  sizeOrb = 2;
  sizeCooNum = 3;
  sizeNs = sizeOrb*sizeL*sizeL;
//  numStepRnd = 26;
  numStepRnd = 176;
//  numLoopRnd = 100000;
  numLoopRnd = 1000000;
  numStepDet = 20;
  numLoopDet = 10000;

  parJx = 0.0;
  parJy = 0.0;
  parJz = 0.0;
  parKx = -1.0;
  parKy = -1.0;
  parKz = -1.0;
  parGx = 0.0;
  parGy = 0.0;
  parGz = 0.0;
  parGxP = 0.0;
  parGyP = 0.0;
  parGzP = 0.0;
  parGK = 0.0;

  // (1,1,1) direction
  parHDir[0] = 1.0/sqrt(3.0);
  parHDir[1] = 1.0/sqrt(3.0);
  parHDir[2] = 1.0/sqrt(3.0);
  // (1,1,-2) direction
//  parHDir[0] = 1.0/sqrt(6.0);
//  parHDir[1] = 1.0/sqrt(6.0);
//  parHDir[2] = -2.0/sqrt(6.0);
  // magnetic field: 0
  parHSize = 0.0;
  parH[0] = 0.0;
  parH[1] = 0.0;
  parH[2] = 0.0;

  parTheta = 0.0;
//  parTmpRatio = 0.5;
  parTmpRatio = 0.9;

  option = -2;// init: rnd
//  mode = 0;// optimization: deterministic
  mode = 1;// optimization: rnd

  // overwrite parameter settings
  for(i = 0; i < argc; ++i){
    if(*argv[i] == '-'){
      opt = *(argv[i]+1);
      if(opt == 'l'){
        sizeL = atoi(argv[i+1]);
        sizeNs = sizeOrb*sizeL*sizeL;
      }else if(opt == 's'){
        numStepRnd = atoi(argv[i+1]);
//      }else if(opt == 't'){
//        parTheta = atof(argv[i+1]);
//        parTheta *= M_PI;
//        parKx = -cos(parTheta);
//        parKy = -cos(parTheta);
//        parKz = -cos(parTheta);
//        parGx = sin(parTheta);
//        parGy = sin(parTheta);
//        parGz = sin(parTheta);
      }else if(opt == 'g'){
        parGK = atof(argv[i+1]);
        parGx = parGK;
        parGy = parGK;
        parGz = parGK;
        parTheta = atan(parGK);
      }else if(opt == 'p'){
        parGP = atof(argv[i+1]);
        parGxP = parGP;
        parGyP = parGP;
        parGzP = parGP;
      }else if(opt == 'H'){// z component of mag field
        parHSize = atof(argv[i+1]);
        parH[0] = - parHSize * parHDir[0];
        parH[1] = - parHSize * parHDir[1];
        parH[2] = - parHSize * parHDir[2];
      }else if(opt == 'r'){
        seed = atoi(argv[i+1]);
      }else if(opt == 'o'){
        option = atoi(argv[i+1]);
      }else if(opt == 'm'){
        mode = atoi(argv[i+1]);
      }
    }
  }

/*
0 initSpinConfGuessFMPMM(spin,sizeL,sizeOrb);
1 initSpinConfGuessZigzag(spin,sizeL,sizeOrb);
2 initSpinConfGuess12sub(spin,sizeL,sizeOrb);
3 initSpinConfGuess6sub(spin,sizeL,sizeOrb);
4 initSpinConfGuess18sub9dir(spin,sizeL,sizeOrb);
5 initSpinConfGuess18sub17dirC3Broken(spin,sizeL,sizeOrb);
6 initSpinConfGuess18sub17dirC3Sym(spin,sizeL,sizeOrb);
7 initSpinConfGuessFMPPP(spin,sizeL,sizeOrb);
8 initSpinConfGuess6sub6dir(spin,sizeL,sizeOrb);
9 initSpinConfGuess50sub(spin,sizeL,sizeOrb);
*/

  printf("# time rnd: %f %d\n",time_count[0],seed);
  printf("# SA cooling rate, num of steps: %16.12f %5d\n",parTmpRatio,numStepRnd);
  printf("# num of metropolis updates in each T: %5d\n",numLoopRnd);
  printf("# Deterministic optimization, num of steps: %16.12f %5d\n",parTmpRatio,numStepDet);
  printf("# num of updates: %5d\n",numLoopDet);
  printf("# mode (1:rnd, 0:deterministic) : %d\n",mode);
  printf("# init option: %d\n",option);
  printf("# Theta: %16.12f\n",parTheta);
  printf("# Theta/Pi: %16.12f\n",parTheta/M_PI);
  printf("# Jx: %16.12f\n",parJx);
  printf("# Jy: %16.12f\n",parJy);
  printf("# Jz: %16.12f\n",parJz);
  printf("# Kx: %16.12f\n",parKx);
  printf("# Ky: %16.12f\n",parKy);
  printf("# Kz: %16.12f\n",parKz);
  printf("# Gx: %16.12f\n",parGx);
  printf("# Gy: %16.12f\n",parGy);
  printf("# Gz: %16.12f\n",parGz);
  printf("# GxP: %16.12f\n",parGxP);
  printf("# GyP: %16.12f\n",parGyP);
  printf("# GzP: %16.12f\n",parGzP);
  printf("# mag field size: %16.12f\n",parHSize);
  printf("# mag field direction: %16.12f %16.12f %16.12f\n",parHDir[0],parHDir[1],parHDir[2]);
  printf("# mag field: %16.12f %16.12f %16.12f\n",parH[0],parH[1],parH[2]);
  dsfmt_init_gen_rand(&dsfmt,seed);

  d_malloc1(locField,3);
  d_malloc2(spin,sizeNs,3);
  d_malloc2(flux,sizeL,sizeL);
  i_malloc3(listSite,sizeL,sizeL,sizeOrb);
  i_malloc2(listLocSite,sizeNs,sizeCooNum);
  d_malloc4(listLocJ,sizeNs,sizeCooNum,3,3);
  d_malloc4(parJ,sizeOrb,sizeCooNum,3,3);
  d_malloc4(listDist,sizeL,sizeL,sizeOrb,2);

  setJ(parJx,parJy,parJz,parKx,parKy,parKz,parGx,parGy,parGz,parGxP,parGyP,parGzP,sizeOrb,parJ);
  makeListXYOrb2Site(sizeL,sizeOrb,listSite);
  makeListLocSite(sizeL,sizeOrb,sizeCooNum,listSite,listLocSite);
  makeListLocJ(sizeOrb,sizeCooNum,sizeNs,parJ,listLocJ);
  initSpinConfGuess(option,spin,sizeL,sizeOrb);

  if(mode==1){
    for(step=0; step<numStepRnd; step++){
      parTmp = pow(parTmpRatio,step);
//      parTmp = pow(parTmpRatio,step)*0.2;
      parInvTmp = 1.0/parTmp;
      for(i=0; i<numLoopRnd; i++){
        for(j=0; j<sizeNs; j++){
          getRndSite(sizeL,sizeOrb,sizeNs,&ix,&iy,&iorb);
          makeSpinCandidate(spinCand);
          calcLocField(sizeL,sizeOrb,sizeCooNum,spin,listLocSite,listLocJ,spinCand,ix,iy,iorb,parH,locField);
          calcBoltzFac(spinCand,spin,listSite,ix,iy,iorb,locField,parInvTmp,&parBF);
          if(parBF > dsfmt_genrand_close_open(&dsfmt)){
            updateSpinCand(spinCand,spin,listSite,ix,iy,iorb);
          }
        }
      }
      calcEne(sizeCooNum,sizeNs,spin,listLocSite,listLocJ,parH,&ene);
      printf("step tmp hz ene: %8d %18.12e %16.12f %16.12f\n",step,parTmp,parHSize,ene);
    }
  }

  printf("\n\n");
  step = 0;
  parTmp = 0.0;
  calcEne(sizeCooNum,sizeNs,spin,listLocSite,listLocJ,parH,&ene);
  printf("step tmp hz ene: %8d %18.12e %16.12f %16.12f\n",step,parTmp,parHSize,ene);
  for(step=1; step<=numStepDet; step++){
    for(i=0; i<numLoopDet; i++){
      for(ixyorb=0; ixyorb<sizeNs; ixyorb++){
        iorb = ixyorb%sizeOrb;
        ixy = ixyorb/sizeOrb;
        ix = ixy%sizeL;
        iy = ixy/sizeL;
        makeSpinCandidateDeterministic(sizeL,sizeOrb,sizeCooNum,spin,listLocSite,listLocJ,spinCand,ix,iy,iorb,parH,locField);
        updateSpinCand(spinCand,spin,listSite,ix,iy,iorb);
      }
    }
    calcEne(sizeCooNum,sizeNs,spin,listLocSite,listLocJ,parH,&ene);
    printf("step tmp hz ene: %8d %18.12e %16.12f %16.12f\n",step,parTmp,parHSize,ene);
  }

  printf("\n\n");
  makeListDist(sizeL,sizeOrb,listDist);
  for(i=0; i<sizeL; i++){
    for(j=0; j<sizeL; j++){
      for(k=0; k<sizeOrb; k++){
        site = listSite[i][j][k];
        dist[0] = listDist[i][j][k][0];
        dist[1] = listDist[i][j][k][1];
        printf("%4d %4d %d  %16.12f %16.12f  %16.12f %16.12f %16.12f\n",
          i,j,k,
          dist[0],dist[1],
          spin[site][0],spin[site][1],spin[site][2]);
      }
    }
  }

  printf("\n\n");
  calcFlux(sizeL,spin,flux);
  for(i=0; i<sizeL; i++){
    for(j=0; j<sizeL; j++){
      printf("%4d %4d  %16.12f\n",
        i,j,flux[i][j]);
    }
  }

  d_free1(locField,3);
  d_free2(spin,sizeNs,3);
  d_free2(flux,sizeL,sizeL);
  i_free3(listSite,sizeL,sizeL,sizeOrb);
  i_free2(listLocSite,sizeNs,sizeCooNum);
  d_free4(listLocJ,sizeNs,sizeCooNum,3,3);
  d_free4(parJ,sizeOrb,sizeCooNum,3,3);
  d_free4(listDist,sizeL,sizeL,sizeOrb,2);

  return 0;
}
