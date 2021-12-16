#ifndef _FUNC_H
#define _FUNC_H

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

dsfmt_t dsfmt;

double gettimeofday_sec();
int initSpinConf(double **spin, int sizeNs);
int initSpinConfMarsaglia(double **spin, int sizeNs);
int initSpinConfGuessFMPMM(double **spin, int sizeL, int sizeOrb);
int initSpinConfGuessZigzag(double **spin, int sizeL, int sizeOrb);
int initSpinConfGuess12sub(double **spin, int sizeL, int sizeOrb);
int initSpinConfGuess6sub(double **spin, int sizeL, int sizeOrb);
int initSpinConfGuess18sub9dir(double **spin, int sizeL, int sizeOrb);
int initSpinConfGuess18sub17dirC3Broken(double **spin, int sizeL, int sizeOrb);
int initSpinConfGuess18sub17dirC3Sym(double **spin, int sizeL, int sizeOrb);
int initSpinConfGuessFMPPP(double **spin, int sizeL, int sizeOrb);
int initSpinConfGuess6sub6dir(double **spin, int sizeL, int sizeOrb);
int initSpinConfGuess50sub(double **spin, int sizeL, int sizeOrb);
int initSpinConfFlipAll(double **spin, int sizeL, int sizeOrb);
int initSpinConfGuess(int option, double **spin, int sizeL, int sizeOrb);
int makeSpinCandidate(double *spin);
int getRndSite(int sizeL, int sizeOrb, int sizeNs, int *ix, int *iy, int *iorb);
int calcBoltzFac(double *spinCand, double **spin, int ***listSite, int ix, int iy, int iorb, double *locField, double parInvTmp, double *parBF);
int makeListXYOrb2Site(int sizeL, int sizeOrb, int ***listSite);
int makeListLocJ(int sizeOrb, int sizeCooNum, int sizeNs, double ****parJ, double ****listLocJ);
int calcLocField(int sizeL, int sizeOrb, int sizeCooNum, double **spin, int **listLocSite, double ****listLocJ, double *spinCand, int ix, int iy, int iorb, double *parH, double *locField);
int makeSpinCandidateDeterministic(int sizeL, int sizeOrb, int sizeCooNum, double **spin, int **listLocSite, double ****listLocJ, double *spinCand, int ix, int iy, int iorb, double *parH, double *locField);
int updateSpinCand(double *spinCand, double **spin, int ***listSite, int ix, int iy, int iorb);
int calcEne(int sizeCooNum, int sizeNs, double **spin, int **listLocSite, double ****listLocJ, double *parH, double *ene);
int calcFlux(int sizeL, double **spin, double **flux);

#endif
