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

double gettimeofday_sec()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

int initSpinConf(
  double **spin,
  int sizeNs
)
{
  int i;
  double theta;
  double phi;
  for(i=0; i<sizeNs; i++){
    theta = acos(1.0 - 2.0 * dsfmt_genrand_close_open(&dsfmt));
    phi = 2.0 * M_PI * dsfmt_genrand_close_open(&dsfmt);
    spin[i][0] = sin(theta)*cos(phi);
    spin[i][1] = sin(theta)*sin(phi);
    spin[i][2] = cos(theta);
  }
  return 0;
}

int initSpinConfMarsaglia(
  double **spin,
  int sizeNs
)
{
  int i;
  double rnd1,rnd2,r2,coef;
  double z;
  for(i=0; i<sizeNs; i++){
    do{
      rnd1 = 2.0 * dsfmt_genrand_close_open(&dsfmt) - 1.0;// uniform in [-1,1)
      rnd2 = 2.0 * dsfmt_genrand_close_open(&dsfmt) - 1.0;// uniform in [-1,1)
      r2 = rnd1*rnd1 + rnd2*rnd2;// uniform in [0,1]
    }while(r2>1.0);
    z = 1.0 - 2.0*r2;// uniform in [-1,1] (d*(r2-0.5) (0<d<2) will be in [-d,d])
    coef = sqrt(2.0*(1.0+z));
    spin[i][0] = coef*rnd1;
    spin[i][1] = coef*rnd2;
    spin[i][2] = z;
  }
  return 0;
}

int initSpinConfGuessFMPMM(
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  int i;
  int sizeNs;
  double tmp;
  sizeNs = sizeOrb*sizeL*sizeL;
  tmp = 1.0/sqrt(3.0);
  for(i=0; i<sizeNs; i++){
    spin[i][0] = tmp;
    spin[i][1] = -tmp;
    spin[i][2] = -tmp;
  }
  return 0;
}

int initSpinConfGuessZigzag(
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  int i,j,sub;
  int site;
  int sizeNs;
  double tmp0,tmp1;
  sizeNs = sizeOrb*sizeL*sizeL;
  tmp0 = 1.0/sqrt(2.0);
  for(i=0; i<sizeL; i++){
  for(j=0; j<sizeL; j++){
  for(sub=0; sub<sizeOrb; sub++){
    site = sizeOrb*(sizeL*j + i) + sub;
    if(j%2==0){
      tmp1 = tmp0;
    }else{
      tmp1 = -tmp0;
    }
    spin[site][0] = tmp1;
    spin[site][1] = 0.0;
    spin[site][2] = tmp1;
  }
  }
  }
//
//  tmp0 = 1.0/sqrt(3.0);
//  for(i=0; i<sizeL; i++){
//  for(j=0; j<sizeL; j++){
//  for(sub=0; sub<sizeOrb; sub++){
//    site = sizeOrb*(sizeL*j + i) + sub;
//    if(j%2==0){
//      tmp1 = tmp0;
//    }else{
//      tmp1 = -tmp0;
//    }
//    spin[site][0] = -tmp1;
//    spin[site][1] = -tmp1;
//    spin[site][2] = +tmp1;
//  }
//  }
//  }
//
  return 0;
}

int initSpinConfGuess12sub(
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  int i,j,sub;
  int site0,site1;
  for(i=0; i<sizeL; i++){
  for(j=0; j<sizeL; j++){
    site0 =  sizeOrb*(sizeL*j + i) + 0;
    site1 =  sizeOrb*(sizeL*j + i) + 1;
    if((i*8+j)%12 == 0){
      spin[site0][0] = 0.706937777074;
      spin[site0][1] = 0.021862266353;
      spin[site0][2] = 0.706937777074;
      spin[site1][0] = -0.706937777074;
      spin[site1][1] = -0.021862266353;
      spin[site1][2] = 0.706937777074;
    }else if((i*8+j)%12 == 1){
      spin[site0][0] = -0.706937777074;
      spin[site0][1] = 0.021862266353;
      spin[site0][2] = 0.706937777074;
      spin[site1][0] = -0.706937777074;
      spin[site1][1] = 0.021862266353;
      spin[site1][2] = 0.706937777074;
    }else if((i*8+j)%12 == 2){
      spin[site0][0] = -0.706937777074;
      spin[site0][1] = -0.021862266353;
      spin[site0][2] = 0.706937777074;
      spin[site1][0] = 0.706937777074;
      spin[site1][1] = 0.021862266353;
      spin[site1][2] = 0.706937777074;
    }else if((i*8+j)%12 == 3){
      spin[site0][0] = -0.706937777074;
      spin[site0][1] = -0.021862266353;
      spin[site0][2] = -0.706937777074;
      spin[site1][0] = 0.706937777074;
      spin[site1][1] = 0.021862266353;
      spin[site1][2] = -0.706937777074;
    }else if((i*8+j)%12 == 4){
      spin[site0][0] = 0.706937777074;
      spin[site0][1] = -0.021862266353;
      spin[site0][2] = -0.706937777074;
      spin[site1][0] = 0.706937777074;
      spin[site1][1] = -0.021862266353;
      spin[site1][2] = -0.706937777074;
    }else if((i*8+j)%12 == 5){
      spin[site0][0] = 0.706937777074;
      spin[site0][1] = 0.021862266353;
      spin[site0][2] = -0.706937777074;
      spin[site1][0] = -0.706937777074;
      spin[site1][1] = -0.021862266353;
      spin[site1][2] = -0.706937777074;
    }else if((i*8+j)%12 == 6){
      spin[site0][0] = 0.706937777074;
      spin[site0][1] = 0.021862266353;
      spin[site0][2] = 0.706937777074;
      spin[site1][0] = -0.706937777074;
      spin[site1][1] = -0.021862266353;
      spin[site1][2] = 0.706937777074;
    }else if((i*8+j)%12 == 7){
      spin[site0][0] = -0.706937777074;
      spin[site0][1] = 0.021862266353;
      spin[site0][2] = 0.706937777074;
      spin[site1][0] = -0.706937777074;
      spin[site1][1] = 0.021862266353;
      spin[site1][2] = 0.706937777074;
    }else if((i*8+j)%12 == 8){
      spin[site0][0] = -0.706937777074;
      spin[site0][1] = -0.021862266353;
      spin[site0][2] = 0.706937777074;
      spin[site1][0] = 0.706937777074;
      spin[site1][1] = 0.021862266353;
      spin[site1][2] = 0.706937777074;
    }else if((i*8+j)%12 == 9){
      spin[site0][0] = -0.706937777074;
      spin[site0][1] = -0.021862266353;
      spin[site0][2] = -0.706937777074;
      spin[site1][0] = 0.706937777074;
      spin[site1][1] = 0.021862266353;
      spin[site1][2] = -0.706937777074;
    }else if((i*8+j)%12 == 10){
      spin[site0][0] = 0.706937777074;
      spin[site0][1] = -0.021862266353;
      spin[site0][2] = -0.706937777074;
      spin[site1][0] = 0.706937777074;
      spin[site1][1] = -0.021862266353;
      spin[site1][2] = -0.706937777074;
    }else if((i*8+j)%12 == 11){
      spin[site0][0] = 0.706937777074;
      spin[site0][1] = 0.021862266353;
      spin[site0][2] = -0.706937777074;
      spin[site1][0] = -0.706937777074;
      spin[site1][1] = -0.021862266353;
      spin[site1][2] = -0.706937777074;
    }
  }
  }
  return 0;
}

int initSpinConfGuess6sub(
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  int i,j,sub;
  int site0,site1;
  for(i=0; i<sizeL; i++){
  for(j=0; j<sizeL; j++){
    site0 =  sizeOrb*(sizeL*j + i) + 0;
    site1 =  sizeOrb*(sizeL*j + i) + 1;
    if((i+j)%3 == 0){
      spin[site0][0] = 0.640374263316;
      spin[site0][1] = 0.640374263316;
      spin[site0][2] = -0.424077358232;
      spin[site1][0] = -0.706607357262;
      spin[site1][1] = -0.706607357262;
      spin[site1][2] = -0.037577723802;
    }else if((i+j)%3 == 1){
      spin[site0][0] = -0.241000101375;
      spin[site0][1] = -0.241000101375;
      spin[site0][2] = 0.940126535247;
      spin[site1][0] = -0.241000101375;
      spin[site1][1] = -0.241000101375;
      spin[site1][2] = 0.940126535247;
    }else if((i+j)%3 == 2){
      spin[site0][0] = -0.706607357262;
      spin[site0][1] = -0.706607357262;
      spin[site0][2] = -0.037577723802;
      spin[site1][0] = 0.640374263316;
      spin[site1][1] = 0.640374263316;
      spin[site1][2] = -0.424077358232;
    }
  }
  }
  return 0;
}

int initSpinConfGuess18sub9dir(
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  int i,j,sub;
  int site0,site1;
  for(i=0; i<sizeL; i++){
  for(j=0; j<sizeL; j++){
    site0 =  sizeOrb*(sizeL*j + i) + 0;
    site1 =  sizeOrb*(sizeL*j + i) + 1;
    if(i%3 == 0 && j%3 == 0){
      spin[site0][0] = 0.664092285024;
      spin[site0][1] = 0.343457237431;
      spin[site0][2] = 0.664092285024;
      spin[site1][0] = -0.664092285024;
      spin[site1][1] = -0.343457237431;
      spin[site1][2] = 0.664092285024;
    }else if(i%3 == 0 && j%3 == 1){
      spin[site0][0] = -0.704364605876;
      spin[site0][1] = -0.087982975505;
      spin[site0][2] = 0.704364605876;
      spin[site1][0] = -0.363594271919;
      spin[site1][1] = 0.857670339266;
      spin[site1][2] = 0.363594271919;
    }else if(i%3 == 0 && j%3 == 2){
      spin[site0][0] = 0.363594271919;
      spin[site0][1] = 0.857670339266;
      spin[site0][2] = -0.363594271919;
      spin[site1][0] = -0.704364605876;
      spin[site1][1] = 0.087982975505;
      spin[site1][2] = -0.704364605876;
    }else if(i%3 == 1 && j%3 == 0){
      spin[site0][0] = -0.664092285024;
      spin[site0][1] = -0.343457237431;
      spin[site0][2] = 0.664092285024;
      spin[site1][0] = 0.664092285024;
      spin[site1][1] = 0.343457237431;
      spin[site1][2] = 0.664092285024;
    }else if(i%3 == 1 && j%3 == 1){
      spin[site0][0] = -0.704364605876;
      spin[site0][1] = 0.087982975505;
      spin[site0][2] = -0.704364605876;
      spin[site1][0] = 0.363594271919;
      spin[site1][1] = 0.857670339266;
      spin[site1][2] = -0.363594271919;
    }else if(i%3 == 1 && j%3 == 2){
      spin[site0][0] = -0.363594271919;
      spin[site0][1] = 0.857670339266;
      spin[site0][2] = 0.363594271919;
      spin[site1][0] = -0.704364605876;
      spin[site1][1] = -0.087982975505;
      spin[site1][2] = 0.704364605876;
    }else if(i%3 == 2 && j%3 == 0){
      spin[site0][0] = 0.664092285024;
      spin[site0][1] = -0.343457237431;
      spin[site0][2] = -0.664092285024;
      spin[site1][0] = 0.664092285024;
      spin[site1][1] = -0.343457237431;
      spin[site1][2] = -0.664092285024;
    }else if(i%3 == 2 && j%3 == 1){
      spin[site0][0] = 0.704364605876;
      spin[site0][1] = -0.087982975505;
      spin[site0][2] = -0.704364605876;
      spin[site1][0] = -0.363594271919;
      spin[site1][1] = -0.857670339266;
      spin[site1][2] = -0.363594271919;
    }else if(i%3 == 2 && j%3 == 2){
      spin[site0][0] = -0.363594271919;
      spin[site0][1] = -0.857670339266;
      spin[site0][2] = -0.363594271919;
      spin[site1][0] = 0.704364605876;
      spin[site1][1] = -0.087982975505;
      spin[site1][2] = -0.704364605876;
    }
  }
  }
  return 0;
}

int initSpinConfGuess18sub17dirC3Broken(
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  int i,j,sub;
  int site0,site1;
  for(i=0; i<sizeL; i++){
  for(j=0; j<sizeL; j++){
    site0 =  sizeOrb*(sizeL*j + i) + 0;
    site1 =  sizeOrb*(sizeL*j + i) + 1;
    if(i%3 == 0 && j%3 == 0){
      spin[site0][0] = -0.632286554682;
      spin[site0][1] = -0.447691216729;
      spin[site0][2] = -0.632286554682;
      spin[site1][0] = 0.387903544014;
      spin[site1][1] = 0.847454373375;
      spin[site1][2] = -0.362425062031;
    }else if(i%3 == 0 && j%3 == 1){
      spin[site0][0] = 0.135328848627;
      spin[site0][1] = 0.700601920754;
      spin[site0][2] = -0.700601920754;
      spin[site1][0] = -0.387903544014;
      spin[site1][1] = -0.362425062031;
      spin[site1][2] = -0.847454373375;
    }else if(i%3 == 0 && j%3 == 2){
      spin[site0][0] = 0.632286554682;
      spin[site0][1] = -0.632286554682;
      spin[site0][2] = 0.447691216729;
      spin[site1][0] = 0.577350269190;
      spin[site1][1] = -0.577350269190;
      spin[site1][2] = 0.577350269190;
    }else if(i%3 == 1 && j%3 == 0){
      spin[site0][0] = -0.284817929625;
      spin[site0][1] = 0.284817929625;
      spin[site0][2] = -0.915290934036;
      spin[site1][0] = -0.577350269190;
      spin[site1][1] = 0.577350269190;
      spin[site1][2] = -0.577350269190;
    }else if(i%3 == 1 && j%3 == 1){
      spin[site0][0] = 0.284817929625;
      spin[site0][1] = 0.915290934036;
      spin[site0][2] = 0.284817929625;
      spin[site1][0] = -0.847454373375;
      spin[site1][1] = -0.362425062031;
      spin[site1][2] = -0.387903544014;
    }else if(i%3 == 1 && j%3 == 2){
      spin[site0][0] = 0.447691216729;
      spin[site0][1] = -0.632286554682;
      spin[site0][2] = 0.632286554682;
      spin[site1][0] = 0.847454373375;
      spin[site1][1] = -0.387903544014;
      spin[site1][2] = 0.362425062031;
    }else if(i%3 == 2 && j%3 == 0){
      spin[site0][0] = -0.915290934036;
      spin[site0][1] = 0.284817929625;
      spin[site0][2] = -0.284817929625;
      spin[site1][0] = -0.362425062031;
      spin[site1][1] = 0.847454373375;
      spin[site1][2] = 0.387903544014;
    }else if(i%3 == 2 && j%3 == 1){
      spin[site0][0] = -0.700601920754;
      spin[site0][1] = 0.700601920754;
      spin[site0][2] = 0.135328848627;
      spin[site1][0] = -0.577350269190;
      spin[site1][1] = 0.577350269190;
      spin[site1][2] = -0.577350269190;
    }else if(i%3 == 2 && j%3 == 2){
      spin[site0][0] = 0.700601920754;
      spin[site0][1] = -0.135328848627;
      spin[site0][2] = 0.700601920754;
      spin[site1][0] = 0.362425062031;
      spin[site1][1] = -0.387903544014;
      spin[site1][2] = 0.847454373375;
    }
  }
  }
  return 0;
}

int initSpinConfGuess18sub17dirC3Sym(
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  int i,j,sub;
  int site0,site1;
  for(i=0; i<sizeL; i++){
  for(j=0; j<sizeL; j++){
    site0 =  sizeOrb*(sizeL*j + i) + 0;
    site1 =  sizeOrb*(sizeL*j + i) + 1;
    if(i%3 == 0 && j%3 == 0){
      spin[site0][0] = 0.577350269190;
      spin[site0][1] = 0.577350269190;
      spin[site0][2] = 0.577350269190;
      spin[site1][0] = -0.680104516679;
      spin[site1][1] = -0.680104516679;
      spin[site1][2] = -0.273707312262;
    }else if(i%3 == 0 && j%3 == 1){
      spin[site0][0] = 0.413991838943;
      spin[site0][1] = -0.804163601407;
      spin[site0][2] = 0.426534476286;
      spin[site1][0] = 0.859480130649;
      spin[site1][1] = -0.361451175831;
      spin[site1][2] = -0.361451175831;
    }else if(i%3 == 0 && j%3 == 2){
      spin[site0][0] = 0.413991838943;
      spin[site0][1] = 0.426534476286;
      spin[site0][2] = -0.804163601407;
      spin[site1][0] = -0.680104516679;
      spin[site1][1] = -0.273707312262;
      spin[site1][2] = -0.680104516679;
    }else if(i%3 == 1 && j%3 == 0){
      spin[site0][0] = -0.804163601407;
      spin[site0][1] = 0.413991838943;
      spin[site0][2] = 0.426534476286;
      spin[site1][0] = -0.361451175831;
      spin[site1][1] = 0.859480130649;
      spin[site1][2] = -0.361451175831;
    }else if(i%3 == 1 && j%3 == 1){
      spin[site0][0] = 0.577350269190;
      spin[site0][1] = 0.577350269190;
      spin[site0][2] = 0.577350269190;
      spin[site1][0] = -0.361451175831;
      spin[site1][1] = -0.361451175831;
      spin[site1][2] = 0.859480130649;
    }else if(i%3 == 1 && j%3 == 2){
      spin[site0][0] = -0.804163601407;
      spin[site0][1] = 0.426534476286;
      spin[site0][2] = 0.413991838943;
      spin[site1][0] = -0.450453113688;
      spin[site1][1] = 0.631304994582;
      spin[site1][2] = 0.631304994582;
    }else if(i%3 == 2 && j%3 == 0){
      spin[site0][0] = 0.426534476286;
      spin[site0][1] = 0.413991838943;
      spin[site0][2] = -0.804163601407;
      spin[site1][0] = -0.273707312262;
      spin[site1][1] = -0.680104516679;
      spin[site1][2] = -0.680104516679;
    }else if(i%3 == 2 && j%3 == 1){
      spin[site0][0] = 0.426534476286;
      spin[site0][1] = -0.804163601407;
      spin[site0][2] = 0.413991838943;
      spin[site1][0] = 0.631304994582;
      spin[site1][1] = -0.450453113688;
      spin[site1][2] = 0.631304994582;
    }else if(i%3 == 2 && j%3 == 2){
      spin[site0][0] = -0.577350269190;
      spin[site0][1] = -0.577350269190;
      spin[site0][2] = -0.577350269190;
      spin[site1][0] = 0.631304994582;
      spin[site1][1] = 0.631304994582;
      spin[site1][2] = -0.450453113688;
    }
  }
  }
  return 0;
}

int initSpinConfGuessFMPPP(
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  int i;
  int sizeNs;
  double tmp;
  sizeNs = sizeOrb*sizeL*sizeL;
  tmp = 1.0/sqrt(3.0);
  for(i=0; i<sizeNs; i++){
    spin[i][0] = tmp;
    spin[i][1] = tmp;
    spin[i][2] = tmp;
  }
  return 0;
}

int initSpinConfGuess6sub6dir(
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  int i,j,sub;
  int site0,site1;
  for(i=0; i<sizeL; i++){
  for(j=0; j<sizeL; j++){
    site0 =  sizeOrb*(sizeL*j + i) + 0;
    site1 =  sizeOrb*(sizeL*j + i) + 1;
    if((i-j+sizeL)%3 == 0){
      spin[site0][0] = 0.674453194025;
      spin[site0][1] = -0.735767268466;
      spin[site0][2] = 0.061314074441;
      spin[site1][0] = 0.735767268466;
      spin[site1][1] = -0.674453194025;
      spin[site1][2] = -0.061314074441;
    }else if((i-j+sizeL)%3 == 1){
      spin[site0][0] = -0.735767268466;
      spin[site0][1] = 0.061314074441;
      spin[site0][2] = 0.674453194025;
      spin[site1][0] = -0.061314074441;
      spin[site1][1] = 0.735767268466;
      spin[site1][2] = -0.674453194025;
    }else if((i-j+sizeL)%3 == 2){
      spin[site0][0] = 0.061314074441;
      spin[site0][1] = 0.674453194025;
      spin[site0][2] = -0.735767268466;
      spin[site1][0] = -0.674453194025;
      spin[site1][1] = -0.061314074441;
      spin[site1][2] = 0.735767268466;
    }
  }
  }
  return 0;
}

int initSpinConfGuess50sub(
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  int i,j,sub;
  int site0,site1;
  for(i=0; i<sizeL; i++){
  for(j=0; j<sizeL; j++){
    site0 =  sizeOrb*(sizeL*j + i) + 0;
    site1 =  sizeOrb*(sizeL*j + i) + 1;
    if(i%5 == 0 && j%5 == 0){
      spin[site0][0] = 0.820540142541;
      spin[site0][1] = 0.425602884860;
      spin[site0][2] = -0.381544307883;
      spin[site1][0] = -0.242365903607;
      spin[site1][1] = -0.928906446819;
      spin[site1][2] = -0.279984967144;
    }else if(i%5 == 0 && j%5 == 1){
      spin[site0][0] = 0.394218076752;
      spin[site0][1] = -0.842783774092;
      spin[site0][2] = 0.366480037776;
      spin[site1][0] = 0.976626059938;
      spin[site1][1] = 0.032651925549;
      spin[site1][2] = 0.212450913879;
    }else if(i%5 == 0 && j%5 == 2){
      spin[site0][0] = 0.158365103778;
      spin[site0][1] = 0.158365103778;
      spin[site0][2] = -0.974597859535;
      spin[site1][0] = 0.158365103778;
      spin[site1][1] = 0.158365103778;
      spin[site1][2] = -0.974597859535;
    }else if(i%5 == 0 && j%5 == 3){
      spin[site0][0] = 0.976626059938;
      spin[site0][1] = 0.032651925549;
      spin[site0][2] = 0.212450913879;
      spin[site1][0] = 0.394218076752;
      spin[site1][1] = -0.842783774092;
      spin[site1][2] = 0.366480037776;
    }else if(i%5 == 0 && j%5 == 4){
      spin[site0][0] = -0.242365903607;
      spin[site0][1] = -0.928906446819;
      spin[site0][2] = -0.279984967144;
      spin[site1][0] = 0.820540142541;
      spin[site1][1] = 0.425602884860;
      spin[site1][2] = -0.381544307883;
    }else if(i%5 == 1 && j%5 == 0){
      spin[site0][0] = -0.339592233782;
      spin[site0][1] = 0.416865330789;
      spin[site0][2] = 0.843149103505;
      spin[site1][0] = 0.039785945936;
      spin[site1][1] = 0.904254101681;
      spin[site1][2] = 0.425137152105;
    }else if(i%5 == 1 && j%5 == 1){
      spin[site0][0] = 0.554388357407;
      spin[site0][1] = 0.570509888809;
      spin[site0][2] = 0.605947205574;
      spin[site1][0] = 0.570509888809;
      spin[site1][1] = 0.554388357407;
      spin[site1][2] = 0.605947205574;
    }else if(i%5 == 1 && j%5 == 2){
      spin[site0][0] = 0.032651925549;
      spin[site0][1] = 0.976626059938;
      spin[site0][2] = 0.212450913879;
      spin[site1][0] = -0.842783774092;
      spin[site1][1] = 0.394218076752;
      spin[site1][2] = 0.366480037776;
    }else if(i%5 == 1 && j%5 == 3){
      spin[site0][0] = 0.172844699739;
      spin[site0][1] = 0.172844699739;
      spin[site0][2] = 0.969664591260;
      spin[site1][0] = 0.586769221104;
      spin[site1][1] = 0.586769221104;
      spin[site1][2] = 0.558035628190;
    }else if(i%5 == 1 && j%5 == 4){
      spin[site0][0] = 0.446090669433;
      spin[site0][1] = 0.827834082308;
      spin[site0][2] = 0.340137982023;
      spin[site1][0] = -0.639125996627;
      spin[site1][1] = 0.273393863376;
      spin[site1][2] = 0.718869776736;
    }else if(i%5 == 2 && j%5 == 0){
      spin[site0][0] = 0.176834526881;
      spin[site0][1] = 0.453330605560;
      spin[site0][2] = -0.873625155410;
      spin[site1][0] = 0.453330605560;
      spin[site1][1] = 0.176834526881;
      spin[site1][2] = -0.873625155410;
    }else if(i%5 == 2 && j%5 == 1){
      spin[site0][0] = 0.904254101681;
      spin[site0][1] = 0.039785945936;
      spin[site0][2] = 0.425137152105;
      spin[site1][0] = 0.416865330789;
      spin[site1][1] = -0.339592233782;
      spin[site1][2] = 0.843149103505;
    }else if(i%5 == 2 && j%5 == 2){
      spin[site0][0] = -0.928906446819;
      spin[site0][1] = -0.242365903607;
      spin[site0][2] = -0.279984967144;
      spin[site1][0] = 0.425602884860;
      spin[site1][1] = 0.820540142541;
      spin[site1][2] = -0.381544307883;
    }else if(i%5 == 2 && j%5 == 3){
      spin[site0][0] = 0.827834082308;
      spin[site0][1] = 0.446090669433;
      spin[site0][2] = 0.340137982023;
      spin[site1][0] = 0.273393863376;
      spin[site1][1] = -0.639125996627;
      spin[site1][2] = 0.718869776736;
    }else if(i%5 == 2 && j%5 == 4){
      spin[site0][0] = -0.686514621978;
      spin[site0][1] = -0.686514621978;
      spin[site0][2] = 0.239573261487;
      spin[site1][0] = 0.693160523417;
      spin[site1][1] = 0.693160523417;
      spin[site1][2] = 0.197628382457;
    }else if(i%5 == 3 && j%5 == 0){
      spin[site0][0] = 0.693160523417;
      spin[site0][1] = 0.693160523417;
      spin[site0][2] = 0.197628382457;
      spin[site1][0] = -0.686514621978;
      spin[site1][1] = -0.686514621978;
      spin[site1][2] = 0.239573261487;
    }else if(i%5 == 3 && j%5 == 1){
      spin[site0][0] = 0.273393863376;
      spin[site0][1] = -0.639125996627;
      spin[site0][2] = 0.718869776736;
      spin[site1][0] = 0.827834082308;
      spin[site1][1] = 0.446090669433;
      spin[site1][2] = 0.340137982023;
    }else if(i%5 == 3 && j%5 == 2){
      spin[site0][0] = 0.425602884860;
      spin[site0][1] = 0.820540142541;
      spin[site0][2] = -0.381544307883;
      spin[site1][0] = -0.928906446819;
      spin[site1][1] = -0.242365903607;
      spin[site1][2] = -0.279984967144;
    }else if(i%5 == 3 && j%5 == 3){
      spin[site0][0] = 0.416865330789;
      spin[site0][1] = -0.339592233782;
      spin[site0][2] = 0.843149103505;
      spin[site1][0] = 0.904254101681;
      spin[site1][1] = 0.039785945936;
      spin[site1][2] = 0.425137152105;
    }else if(i%5 == 3 && j%5 == 4){
      spin[site0][0] = 0.453330605560;
      spin[site0][1] = 0.176834526881;
      spin[site0][2] = -0.873625155410;
      spin[site1][0] = 0.176834526881;
      spin[site1][1] = 0.453330605560;
      spin[site1][2] = -0.873625155410;
    }else if(i%5 == 4 && j%5 == 0){
      spin[site0][0] = -0.639125996627;
      spin[site0][1] = 0.273393863376;
      spin[site0][2] = 0.718869776736;
      spin[site1][0] = 0.446090669433;
      spin[site1][1] = 0.827834082308;
      spin[site1][2] = 0.340137982023;
    }else if(i%5 == 4 && j%5 == 1){
      spin[site0][0] = 0.586769221104;
      spin[site0][1] = 0.586769221104;
      spin[site0][2] = 0.558035628190;
      spin[site1][0] = 0.172844699739;
      spin[site1][1] = 0.172844699739;
      spin[site1][2] = 0.969664591260;
    }else if(i%5 == 4 && j%5 == 2){
      spin[site0][0] = -0.842783774092;
      spin[site0][1] = 0.394218076752;
      spin[site0][2] = 0.366480037776;
      spin[site1][0] = 0.032651925549;
      spin[site1][1] = 0.976626059938;
      spin[site1][2] = 0.212450913879;
    }else if(i%5 == 4 && j%5 == 3){
      spin[site0][0] = 0.570509888809;
      spin[site0][1] = 0.554388357407;
      spin[site0][2] = 0.605947205574;
      spin[site1][0] = 0.554388357407;
      spin[site1][1] = 0.570509888809;
      spin[site1][2] = 0.605947205574;
    }else if(i%5 == 4 && j%5 == 4){
      spin[site0][0] = 0.039785945936;
      spin[site0][1] = 0.904254101681;
      spin[site0][2] = 0.425137152105;
      spin[site1][0] = -0.339592233782;
      spin[site1][1] = 0.416865330789;
      spin[site1][2] = 0.843149103505;
    }
  }
  }
  return 0;
}


int initSpinConfFlipAll(
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  int i;
  int sizeNs;
  sizeNs = sizeOrb*sizeL*sizeL;
  for(i=0; i<sizeNs; i++){
    spin[i][0] *= -1.0;
    spin[i][1] *= -1.0;
    spin[i][2] *= -1.0;
  }
  return 0;
}

int initSpinConfGuess(
  int option,
  double **spin,
  int sizeL,
  int sizeOrb
)
{
  if(option==0){
    initSpinConfGuessFMPMM(spin,sizeL,sizeOrb);
  }else if(option==1){
    initSpinConfGuessZigzag(spin,sizeL,sizeOrb);
  }else if(option==2){
    initSpinConfGuess12sub(spin,sizeL,sizeOrb);
  }else if(option==3){
    initSpinConfGuess6sub(spin,sizeL,sizeOrb);
  }else if(option==4){
    initSpinConfGuess18sub9dir(spin,sizeL,sizeOrb);
  }else if(option==5){
    initSpinConfGuess18sub17dirC3Broken(spin,sizeL,sizeOrb);
  }else if(option==6){
    initSpinConfGuess18sub17dirC3Sym(spin,sizeL,sizeOrb);
  }else if(option==7){
    initSpinConfGuessFMPPP(spin,sizeL,sizeOrb);
  }else if(option==8){
    initSpinConfGuess6sub6dir(spin,sizeL,sizeOrb);
  }else if(option==9){
    initSpinConfGuess50sub(spin,sizeL,sizeOrb);
//
  }else if(option==100){
    initSpinConfGuessFMPMM(spin,sizeL,sizeOrb);
    initSpinConfFlipAll(spin,sizeL,sizeOrb);
  }else if(option==101){
    initSpinConfGuessZigzag(spin,sizeL,sizeOrb);
    initSpinConfFlipAll(spin,sizeL,sizeOrb);
  }else if(option==102){
    initSpinConfGuess12sub(spin,sizeL,sizeOrb);
    initSpinConfFlipAll(spin,sizeL,sizeOrb);
  }else if(option==103){
    initSpinConfGuess6sub(spin,sizeL,sizeOrb);
    initSpinConfFlipAll(spin,sizeL,sizeOrb);
  }else if(option==104){
    initSpinConfGuess18sub9dir(spin,sizeL,sizeOrb);
    initSpinConfFlipAll(spin,sizeL,sizeOrb);
  }else if(option==105){
    initSpinConfGuess18sub17dirC3Broken(spin,sizeL,sizeOrb);
    initSpinConfFlipAll(spin,sizeL,sizeOrb);
  }else if(option==106){
    initSpinConfGuess18sub17dirC3Sym(spin,sizeL,sizeOrb);
    initSpinConfFlipAll(spin,sizeL,sizeOrb);
  }else if(option==107){
    initSpinConfGuessFMPPP(spin,sizeL,sizeOrb);
    initSpinConfFlipAll(spin,sizeL,sizeOrb);
  }else if(option==108){
    initSpinConfGuess6sub6dir(spin,sizeL,sizeOrb);
    initSpinConfFlipAll(spin,sizeL,sizeOrb);
  }else if(option==108){
    initSpinConfGuess6sub6dir(spin,sizeL,sizeOrb);
    initSpinConfFlipAll(spin,sizeL,sizeOrb);
  }else if(option==109){
    initSpinConfGuess50sub(spin,sizeL,sizeOrb);
    initSpinConfFlipAll(spin,sizeL,sizeOrb);
//
  }else{
    int sizeNs;
    sizeNs = sizeOrb*sizeL*sizeL;
    initSpinConfMarsaglia(spin,sizeNs);
  }
  return 0;
}

int makeSpinCandidate(
  double *spin_site
){
  double rnd1,rnd2,r2,coef;
  double z;
  do{
    rnd1 = 2.0 * dsfmt_genrand_close_open(&dsfmt) - 1.0;// uniform in [-1,1)
    rnd2 = 2.0 * dsfmt_genrand_close_open(&dsfmt) - 1.0;// uniform in [-1,1)
    r2 = rnd1*rnd1 + rnd2*rnd2;// uniform in [0,1]
  }while(r2>1.0);
  z = 1.0 - 2.0*r2;// uniform in [-1,1] (d*(r2-0.5) (0<d<2) will be in [-d,d])
  coef = sqrt(2.0*(1.0+z));
  spin_site[0] = coef*rnd1;
  spin_site[1] = coef*rnd2;
  spin_site[2] = z;
  return 0;
}

int getRndSite(
  int sizeL,
  int sizeOrb,
  int sizeNs,
  int *ix,
  int *iy,
  int *iorb
)
{
  int ixyorb;
  int ixy;
  ixyorb = (int) (sizeNs * dsfmt_genrand_close_open(&dsfmt));
  iorb[0] = ixyorb%sizeOrb;
  ixy = ixyorb/sizeOrb;
  ix[0] = ixy%sizeL;
  iy[0] = ixy/sizeL;
  return 0;
}

int calcBoltzFac(
  double *spinCand,
  double **spin,
  int ***listSite,
  int ix,
  int iy,
  int iorb,
  double *locField,
  double parInvTmp,
  double *parBF
){
  int site;
  int i,j;
  double eneDiff;
  double spinDiff;
  site = listSite[ix][iy][iorb];
  eneDiff = 0.0;
  for(i=0; i<3; i++){
    eneDiff += (spin[site][i] - spinCand[i]) * locField[i];
  }
  parBF[0] = exp(parInvTmp * eneDiff);
  return 0;
}

int makeListXYOrb2Site(
  int sizeL,
  int sizeOrb,
  int ***listSite
){
  int i,j,k;
  for(i=0; i<sizeL; i++){
    for(j=0; j<sizeL; j++){
      for(k=0; k<sizeOrb; k++){
        listSite[i][j][k] = sizeOrb*(sizeL*j + i) + k;
      }
    }
  }
  return 0;
}

int makeListLocJ(
  int sizeOrb,
  int sizeCooNum,
  int sizeNs,
  double ****parJ,
  double ****listLocJ
){
  int i,j,k,l,orb;
  for(i=0; i<sizeNs; i++){
    orb = i%sizeOrb;
    for(j=0; j<sizeCooNum; j++){
      for(k=0; k<3; k++){
        for(l=0; l<3; l++){
          listLocJ[i][j][k][l] = parJ[orb][j][k][l];
        }
      }
    }
  }
  return 0;
}

int calcLocField(
  int sizeL,
  int sizeOrb,
  int sizeCooNum,
  double **spin,
  int **listLocSite,
  double ****listLocJ,
  double *spinCand,
  int ix,
  int iy,
  int iorb,
  double *parH,
  double *locField
){
  int i,j,k,site0,site1;
  site0 = sizeOrb*(sizeL*iy + ix) + iorb;
  for(i=0; i<3; i++){
    locField[i] = parH[i];
  }
  for(i=0; i<sizeCooNum; i++){
    site1 = listLocSite[site0][i];
    for(j=0; j<3; j++){
      for(k=0; k<3; k++){
        locField[j] += spin[site1][k] * listLocJ[site0][i][j][k];
      }
    }
  }
  return 0;
}

int makeSpinCandidateDeterministic(
  int sizeL,
  int sizeOrb,
  int sizeCooNum,
  double **spin,
  int **listLocSite,
  double ****listLocJ,
  double *spinCand,
  int ix,
  int iy,
  int iorb,
  double *parH,
  double *locField
){
  int i,j,k,site0,site1;
  double norm;
  site0 = sizeOrb*(sizeL*iy + ix) + iorb;
  for(i=0; i<3; i++){
    locField[i] = parH[i];
  }
  for(i=0; i<sizeCooNum; i++){
    site1 = listLocSite[site0][i];
    for(j=0; j<3; j++){
      for(k=0; k<3; k++){
        locField[j] += spin[site1][k] * listLocJ[site0][i][j][k];
      }
    }
  }
  norm = 0.0;
  for(j=0; j<3; j++){
    norm += locField[j]*locField[j];
  }
  norm = - 1.0/sqrt(norm);
  for(j=0; j<3; j++){
    spinCand[j] = locField[j] * norm;
  }
  return 0;
}

int updateSpinCand(
  double *spinCand,
  double **spin,
  int ***listSite,
  int ix,
  int iy,
  int iorb
){
  int site;
  int i;
  site = listSite[ix][iy][iorb];
  for(i=0; i<3; i++){
    spin[site][i] = spinCand[i];
  }
  return 0;
}

int calcEne(
  int sizeCooNum,
  int sizeNs,
  double **spin,
  int **listLocSite,
  double ****listLocJ,
  double *parH,
  double *ene
){
  int i,j,k,l,site;
  double e0,e1;
  double ssk;
  e0 = 0.0;
  e1 = 0.0;
  for(i=0; i<sizeNs; i++){
    for(k=0; k<3; k++){
      e0 += spin[i][k] * parH[k];
    }
    for(j=0; j<sizeCooNum; j++){
      site = listLocSite[i][j];
      for(k=0; k<3; k++){
        for(l=0; l<3; l++){
          e1 += spin[i][k] * spin[site][l] * listLocJ[i][j][k][l];
        }
      }
    }
  }
  e0 /= (1.0*sizeNs);
  e1 /= (2.0*sizeNs);
  ene[0] = e0 + e1;
  return 0;
}

int calcFlux(
  int sizeL,
  double **spin,
  double **flux
){
  int i,j;
  int ip,jp;
  double s[6];
  for(i=0; i<sizeL; i++){
    ip = (i+1+sizeL)%sizeL;
    for(j=0; j<sizeL; j++){
      jp = (j+1+sizeL)%sizeL;
      s[0] = spin[2*(sizeL*j+i)+1][2];
      s[1] = spin[2*(sizeL*jp+i)+0][0];
      s[2] = spin[2*(sizeL*jp+i)+1][1];
      s[3] = spin[2*(sizeL*jp+ip)+0][2];
      s[4] = spin[2*(sizeL*j+ip)+1][0];
      s[5] = spin[2*(sizeL*j+ip)+0][1];
      flux[i][j] = s[0]*s[1]*s[2]*s[3]*s[4]*s[5];
    }
  }
  return 0;
}
