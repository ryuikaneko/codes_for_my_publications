CC       = gcc
WARN     =
#WARN     = -Wmissing-prototypes -Wall -Winline 
CFLAGSMT = -O3 -finline-functions -fomit-frame-pointer -DNDEBUG -fno-strict-aliasing \
 --param max-inline-insns-single=1800 \
 --param inline-unit-growth=500 \
 --param large-function-growth=900 \
 -std=c99 -lm
CFLAGS   = -O3 -msse2 -DHAVE_SSE2 -DDSFMT_MEXP=19937 -std=c99 -lm
OBJS1    = triangle
OBJSMT   = ../dSFMT/dSFMT

$(OBJS1).out: $(OBJSMT).o $(OBJS1).o
	$(CC) $(CFLAGSMT) $(WARN) -o $@ $(OBJSMT).o $(OBJS1).o

$(OBJS1).o: $(OBJS1).c
	$(CC) $(CFLAGS) $(WARN) -c -o $(OBJS1).o $(OBJS1).c

SUFFIXES: .c .o

.c .o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o ../dSFMT/*.o *~ ../dSFMT/*~

remove:
	rm -f $(OBJS1).out *.o ../dSFMT/*.o *~ ../dSFMT/*~
