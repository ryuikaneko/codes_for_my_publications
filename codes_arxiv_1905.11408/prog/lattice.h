#ifndef _LATTICE_H
#define _LATTICE_H

int setJ(double parJx, double parJy, double parJz, double parKx, double parKy, double parKz, double parGx, double parGy, double parGz, double parGxP, double parGyP, double parGzP, int sizeOrb, double ****parJ);
int makeListLocSite(int sizeL, int sizeOrb, int sizeCooNum, int ***listSite, int **listLocSite);
int makeListDist(int sizeL, int sizeOrb, double ****listDist);

#endif
