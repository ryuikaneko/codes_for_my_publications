#ifndef _MFMEMORY_H
#define _MFMEMORY_H

int mfint[7];

/* int type ============================================================== */
#define i_malloc1(X, N1) \
  X = (int*)malloc((N1)*sizeof(int));

#define i_malloc2(X, N1, N2) \
  X = (int**)malloc((N1)*sizeof(int*));\
  X[0] = (int*)malloc((N1)*(N2)*sizeof(int));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
  }

#define i_malloc3(X, N1, N2, N3) \
  X = (int***)malloc((N1)*sizeof(int**));\
  X[0] = (int**)malloc((N1)*(N2)*sizeof(int*));\
  X[0][0] = (int*)malloc((N1)*(N2)*(N3)*sizeof(int));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
    }\
  }

#define i_malloc4(X, N1, N2, N3, N4) \
  X = (int****)malloc((N1)*sizeof(int***));\
  X[0] = (int***)malloc((N1)*(N2)*sizeof(int**));\
  X[0][0] = (int**)malloc((N1)*(N2)*(N3)*sizeof(int*));\
  X[0][0][0] = (int*)malloc((N1)*(N2)*(N3)*(N4)*sizeof(int));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
      for(mfint[2]=0; mfint[2]<(N3); mfint[2]++){\
        X[mfint[0]][mfint[1]][mfint[2]] = X[0][0][0] + ((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4);\
      }\
    }\
  }

#define i_malloc5(X, N1, N2, N3, N4, N5) \
  X = (int*****)malloc((N1)*sizeof(int****));\
  X[0] = (int****)malloc((N1)*(N2)*sizeof(int***));\
  X[0][0] = (int***)malloc((N1)*(N2)*(N3)*sizeof(int**));\
  X[0][0][0] = (int**)malloc((N1)*(N2)*(N3)*(N4)*sizeof(int*));\
  X[0][0][0][0] = (int*)malloc((N1)*(N2)*(N3)*(N4)*(N5)*sizeof(int));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
      for(mfint[2]=0; mfint[2]<(N3); mfint[2]++){\
        X[mfint[0]][mfint[1]][mfint[2]] = X[0][0][0] + ((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4);\
        for(mfint[3]=0; mfint[3]<(N4); mfint[3]++){\
          X[mfint[0]][mfint[1]][mfint[2]][mfint[3]] = X[0][0][0][0] + (((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4) + mfint[3])*(N5);\
        }\
      }\
    }\
  }

#define i_malloc6(X, N1, N2, N3, N4, N5, N6) \
  X = (int******)malloc((N1)*sizeof(int*****));\
  X[0] = (int*****)malloc((N1)*(N2)*sizeof(int****));\
  X[0][0] = (int****)malloc((N1)*(N2)*(N3)*sizeof(int***));\
  X[0][0][0] = (int***)malloc((N1)*(N2)*(N3)*(N4)*sizeof(int**));\
  X[0][0][0][0] = (int**)malloc((N1)*(N2)*(N3)*(N4)*(N5)*sizeof(int*));\
  X[0][0][0][0][0] = (int*)malloc((N1)*(N2)*(N3)*(N4)*(N5)*(N6)*sizeof(int));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
      for(mfint[2]=0; mfint[2]<(N3); mfint[2]++){\
        X[mfint[0]][mfint[1]][mfint[2]] = X[0][0][0] + ((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4);\
        for(mfint[3]=0; mfint[3]<(N4); mfint[3]++){\
          X[mfint[0]][mfint[1]][mfint[2]][mfint[3]] = X[0][0][0][0] + (((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4) + mfint[3])*(N5);\
          for(mfint[4]=0; mfint[4]<(N5); mfint[4]++){\
            X[mfint[0]][mfint[1]][mfint[2]][mfint[3]][mfint[4]] = X[0][0][0][0][0] + ((((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4) + mfint[3])*(N5) + mfint[4])*(N6);\
          }\
        }\
      }\
    }\
  }

#define i_free1(X, N1) \
  free(X);

#define i_free2(X, N1, N2) \
  free(X[0]);\
  free(X);

#define i_free3(X, N1, N2, N3) \
  free(X[0][0]);\
  free(X[0]);\
  free(X);

#define i_free4(X, N1, N2, N3, N4) \
  free(X[0][0][0]);\
  free(X[0][0]);\
  free(X[0]);\
  free(X);

#define i_free5(X, N1, N2, N3, N4, N5) \
  free(X[0][0][0][0]);\
  free(X[0][0][0]);\
  free(X[0][0]);\
  free(X[0]);\
  free(X);

#define i_free6(X, N1, N2, N3, N4, N5, N6) \
  free(X[0][0][0][0][0]);\
  free(X[0][0][0][0]);\
  free(X[0][0][0]);\
  free(X[0][0]);\
  free(X[0]);\
  free(X);
/* end of int type ==============================================================*/

/* double type ==============================================================*/
#define d_malloc1(X, N1) \
  X = (double*)malloc((N1)*sizeof(double));

#define d_malloc2(X, N1, N2) \
  X = (double**)malloc((N1)*sizeof(double*));\
  X[0] = (double*)malloc((N1)*(N2)*sizeof(double));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
  }

#define d_malloc3(X, N1, N2, N3) \
  X = (double***)malloc((N1)*sizeof(double**));\
  X[0] = (double**)malloc((N1)*(N2)*sizeof(double*));\
  X[0][0] = (double*)malloc((N1)*(N2)*(N3)*sizeof(double));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
    }\
  }

#define d_malloc4(X, N1, N2, N3, N4) \
  X = (double****)malloc((N1)*sizeof(double***));\
  X[0] = (double***)malloc((N1)*(N2)*sizeof(double**));\
  X[0][0] = (double**)malloc((N1)*(N2)*(N3)*sizeof(double*));\
  X[0][0][0] = (double*)malloc((N1)*(N2)*(N3)*(N4)*sizeof(double));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
      for(mfint[2]=0; mfint[2]<(N3); mfint[2]++){\
        X[mfint[0]][mfint[1]][mfint[2]] = X[0][0][0] + ((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4);\
      }\
    }\
  }

#define d_malloc5(X, N1, N2, N3, N4, N5) \
  X = (double*****)malloc((N1)*sizeof(double****));\
  X[0] = (double****)malloc((N1)*(N2)*sizeof(double***));\
  X[0][0] = (double***)malloc((N1)*(N2)*(N3)*sizeof(double**));\
  X[0][0][0] = (double**)malloc((N1)*(N2)*(N3)*(N4)*sizeof(double*));\
  X[0][0][0][0] = (double*)malloc((N1)*(N2)*(N3)*(N4)*(N5)*sizeof(double));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
      for(mfint[2]=0; mfint[2]<(N3); mfint[2]++){\
        X[mfint[0]][mfint[1]][mfint[2]] = X[0][0][0] + ((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4);\
        for(mfint[3]=0; mfint[3]<(N4); mfint[3]++){\
          X[mfint[0]][mfint[1]][mfint[2]][mfint[3]] = X[0][0][0][0] + (((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4) + mfint[3])*(N5);\
        }\
      }\
    }\
  }

#define d_malloc6(X, N1, N2, N3, N4, N5, N6) \
  X = (double******)malloc((N1)*sizeof(double*****));\
  X[0] = (double*****)malloc((N1)*(N2)*sizeof(double****));\
  X[0][0] = (double****)malloc((N1)*(N2)*(N3)*sizeof(double***));\
  X[0][0][0] = (double***)malloc((N1)*(N2)*(N3)*(N4)*sizeof(double**));\
  X[0][0][0][0] = (double**)malloc((N1)*(N2)*(N3)*(N4)*(N5)*sizeof(double*));\
  X[0][0][0][0][0] = (double*)malloc((N1)*(N2)*(N3)*(N4)*(N5)*(N6)*sizeof(double));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
      for(mfint[2]=0; mfint[2]<(N3); mfint[2]++){\
        X[mfint[0]][mfint[1]][mfint[2]] = X[0][0][0] + ((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4);\
        for(mfint[3]=0; mfint[3]<(N4); mfint[3]++){\
          X[mfint[0]][mfint[1]][mfint[2]][mfint[3]] = X[0][0][0][0] + (((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4) + mfint[3])*(N5);\
          for(mfint[4]=0; mfint[4]<(N5); mfint[4]++){\
            X[mfint[0]][mfint[1]][mfint[2]][mfint[3]][mfint[4]] = X[0][0][0][0][0] + ((((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4) + mfint[3])*(N5) + mfint[4])*(N6);\
          }\
        }\
      }\
    }\
  }

#define d_free1(X, N1) \
  free(X);

#define d_free2(X, N1, N2) \
  free(X[0]);\
  free(X);

#define d_free3(X, N1, N2, N3) \
  free(X[0][0]);\
  free(X[0]);\
  free(X);

#define d_free4(X, N1, N2, N3, N4) \
  free(X[0][0][0]);\
  free(X[0][0]);\
  free(X[0]);\
  free(X);

#define d_free5(X, N1, N2, N3, N4, N5) \
  free(X[0][0][0][0]);\
  free(X[0][0][0]);\
  free(X[0][0]);\
  free(X[0]);\
  free(X);

#define d_free6(X, N1, N2, N3, N4, N5, N6) \
  free(X[0][0][0][0][0]);\
  free(X[0][0][0][0]);\
  free(X[0][0][0]);\
  free(X[0][0]);\
  free(X[0]);\
  free(X);
/* end of double type ==============================================================*/

/* double complex type ==============================================================*/
#define dc_malloc1(X, N1) \
  X = (double complex*)malloc((N1)*sizeof(double complex));

#define dc_malloc2(X, N1, N2) \
  X = (double complex**)malloc((N1)*sizeof(double complex*));\
  X[0] = (double complex*)malloc((N1)*(N2)*sizeof(double complex));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
  }

#define dc_malloc3(X, N1, N2, N3) \
  X = (double complex***)malloc((N1)*sizeof(double complex**));\
  X[0] = (double complex**)malloc((N1)*(N2)*sizeof(double complex*));\
  X[0][0] = (double complex*)malloc((N1)*(N2)*(N3)*sizeof(double complex));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
    }\
  }

#define dc_malloc4(X, N1, N2, N3, N4) \
  X = (double complex****)malloc((N1)*sizeof(double complex***));\
  X[0] = (double complex***)malloc((N1)*(N2)*sizeof(double complex**));\
  X[0][0] = (double complex**)malloc((N1)*(N2)*(N3)*sizeof(double complex*));\
  X[0][0][0] = (double complex*)malloc((N1)*(N2)*(N3)*(N4)*sizeof(double complex));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
      for(mfint[2]=0; mfint[2]<(N3); mfint[2]++){\
        X[mfint[0]][mfint[1]][mfint[2]] = X[0][0][0] + ((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4);\
      }\
    }\
  }

#define dc_malloc5(X, N1, N2, N3, N4, N5) \
  X = (double complex*****)malloc((N1)*sizeof(double complex****));\
  X[0] = (double complex****)malloc((N1)*(N2)*sizeof(double complex***));\
  X[0][0] = (double complex***)malloc((N1)*(N2)*(N3)*sizeof(double complex**));\
  X[0][0][0] = (double complex**)malloc((N1)*(N2)*(N3)*(N4)*sizeof(double complex*));\
  X[0][0][0][0] = (double complex*)malloc((N1)*(N2)*(N3)*(N4)*(N5)*sizeof(double complex));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
      for(mfint[2]=0; mfint[2]<(N3); mfint[2]++){\
        X[mfint[0]][mfint[1]][mfint[2]] = X[0][0][0] + ((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4);\
        for(mfint[3]=0; mfint[3]<(N4); mfint[3]++){\
          X[mfint[0]][mfint[1]][mfint[2]][mfint[3]] = X[0][0][0][0] + (((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4) + mfint[3])*(N5);\
        }\
      }\
    }\
  }

#define dc_malloc6(X, N1, N2, N3, N4, N5, N6) \
  X = (double complex******)malloc((N1)*sizeof(double complex*****));\
  X[0] = (double complex*****)malloc((N1)*(N2)*sizeof(double complex****));\
  X[0][0] = (double complex****)malloc((N1)*(N2)*(N3)*sizeof(double complex***));\
  X[0][0][0] = (double complex***)malloc((N1)*(N2)*(N3)*(N4)*sizeof(double complex**));\
  X[0][0][0][0] = (double complex**)malloc((N1)*(N2)*(N3)*(N4)*(N5)*sizeof(double complex*));\
  X[0][0][0][0][0] = (double complex*)malloc((N1)*(N2)*(N3)*(N4)*(N5)*(N6)*sizeof(double complex));\
  for(mfint[0]=0; mfint[0]<(N1); mfint[0]++){\
    X[mfint[0]] = X[0] + mfint[0]*(N2);\
    for(mfint[1]=0; mfint[1]<(N2); mfint[1]++){\
      X[mfint[0]][mfint[1]] = X[0][0] + (mfint[0]*(N2) + mfint[1])*(N3);\
      for(mfint[2]=0; mfint[2]<(N3); mfint[2]++){\
        X[mfint[0]][mfint[1]][mfint[2]] = X[0][0][0] + ((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4);\
        for(mfint[3]=0; mfint[3]<(N4); mfint[3]++){\
          X[mfint[0]][mfint[1]][mfint[2]][mfint[3]] = X[0][0][0][0] + (((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4) + mfint[3])*(N5);\
          for(mfint[4]=0; mfint[4]<(N5); mfint[4]++){\
            X[mfint[0]][mfint[1]][mfint[2]][mfint[3]][mfint[4]] = X[0][0][0][0][0] + ((((mfint[0]*(N2) + mfint[1])*(N3) + mfint[2])*(N4) + mfint[3])*(N5) + mfint[4])*(N6);\
          }\
        }\
      }\
    }\
  }

#define dc_free1(X, N1) \
  free(X);

#define dc_free2(X, N1, N2) \
  free(X[0]);\
  free(X);

#define dc_free3(X, N1, N2, N3) \
  free(X[0][0]);\
  free(X[0]);\
  free(X);

#define dc_free4(X, N1, N2, N3, N4) \
  free(X[0][0][0]);\
  free(X[0][0]);\
  free(X[0]);\
  free(X);

#define dc_free5(X, N1, N2, N3, N4, N5) \
  free(X[0][0][0][0]);\
  free(X[0][0][0]);\
  free(X[0][0]);\
  free(X[0]);\
  free(X);

#define dc_free6(X, N1, N2, N3, N4, N5, N6) \
  free(X[0][0][0][0][0]);\
  free(X[0][0][0][0]);\
  free(X[0][0][0]);\
  free(X[0][0]);\
  free(X[0]);\
  free(X);
/* end of double complex type ==============================================================*/

#endif
