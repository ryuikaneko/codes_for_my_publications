#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include <complex.h>
#include "./dSFMT/dSFMT.h"
#include "mfmemory.h"
#include "lattice.h"

int setJ(
  double parJx,
  double parJy,
  double parJz,
  double parKx,
  double parKy,
  double parKz,
  double parGx,
  double parGy,
  double parGz,
  double parGxP,
  double parGyP,
  double parGzP,
  int sizeOrb,
  double ****parJ
){
  int i,j,k;
  for(i=0; i<sizeOrb; i++){
    parJ[i][0][0][0] = parKx + parJx;
    parJ[i][0][0][1] = parGxP;
    parJ[i][0][0][2] = parGxP;
    parJ[i][0][1][0] = parGxP;
    parJ[i][0][1][1] = parJy;
    parJ[i][0][1][2] = parGx;
    parJ[i][0][2][0] = parGxP;
    parJ[i][0][2][1] = parGx;
    parJ[i][0][2][2] = parJz;

    parJ[i][1][0][0] = parJx;
    parJ[i][1][0][1] = parGyP;
    parJ[i][1][0][2] = parGy;
    parJ[i][1][1][0] = parGyP;
    parJ[i][1][1][1] = parKy + parJy;
    parJ[i][1][1][2] = parGyP;
    parJ[i][1][2][0] = parGy;
    parJ[i][1][2][1] = parGyP;
    parJ[i][1][2][2] = parJz;

    parJ[i][2][0][0] = parJx;
    parJ[i][2][0][1] = parGz;
    parJ[i][2][0][2] = parGzP;
    parJ[i][2][1][0] = parGz;
    parJ[i][2][1][1] = parJy;
    parJ[i][2][1][2] = parGzP;
    parJ[i][2][2][0] = parGzP;
    parJ[i][2][2][1] = parGzP;
    parJ[i][2][2][2] = parKz + parJz;
  }
  return 0;
}

int makeListLocSite(
  int sizeL,
  int sizeOrb,
  int sizeCooNum,
  int ***listSite,
  int **listLocSite
){
  int i,j,k,l;
  int ip,im,jp,jm;
  int site[sizeOrb][sizeCooNum+1];
  int site0;
  for(i=0; i<sizeL; i++){
    ip = (i+1+sizeL)%sizeL;
    im = (i-1+sizeL)%sizeL;
    for(j=0; j<sizeL; j++){
      jp = (j+1+sizeL)%sizeL;
      jm = (j-1+sizeL)%sizeL;

      site[0][0] = listSite[i][j][0];
      site[0][1] = listSite[im][j][1];
      site[0][2] = listSite[i][jm][1];
      site[0][3] = listSite[i][j][1];

      site[1][0] = listSite[i][j][1];
      site[1][1] = listSite[ip][j][0];
      site[1][2] = listSite[i][jp][0];
      site[1][3] = listSite[i][j][0];

      for(k=0; k<sizeOrb; k++){
        site0 = site[k][0];
        for(l=0; l<sizeCooNum; l++){
          listLocSite[site0][l] = site[k][l+1];
        }
      }
    }
  }
  return 0;
}

int makeListDist(
  int sizeL,
  int sizeOrb,
  double ****listDist
){
  int i,j,k,l;
  double Rintra[sizeOrb][2];
  double Rinter[2][2];

  Rintra[0][0] = 0.0;
  Rintra[0][1] = 0.0;
  Rintra[1][0] = 0.5*sqrt(3.0);
  Rintra[1][1] = 0.5;

  Rinter[0][0] = sqrt(3.0);
  Rinter[0][1] = 0.0;
  Rinter[1][0] = 0.5*sqrt(3.0);
  Rinter[1][1] = 1.5;

  for(i=0; i<sizeL; i++){
    for(j=0; j<sizeL; j++){
      for(k=0; k<sizeOrb; k++){
        for(l=0; l<2; l++){
          listDist[i][j][k][l] = i * Rinter[0][l] + j * Rinter[1][l] + Rintra[k][l];
        }
      }
    }
  }
  return 0;
}
