CC = gcc
WARN = -Wmissing-prototypes -Wall -Winline 
CFLAGSFMT = -O3 -finline-functions -fomit-frame-pointer -DNDEBUG -fno-strict-aliasing \
 --param max-inline-insns-single=1800 \
 --param inline-unit-growth=500 \
 --param large-function-growth=900 \
 -std=c99 -lm
CFLAGS = -O3 -msse2 -DHAVE_SSE2 -DDSFMT_MEXP=19937 -std=c99 -lm
OBJSFMT = ./dSFMT/dSFMT

main.out: $(OBJSFMT).o func.o lattice.o main.o
	$(CC) $(CFLAGSFMT) $(WARN) -o $@ $(OBJSFMT).o func.o lattice.o main.o

SUFFIXES: .c .o

.c .o:
	$(CC) $(CFLAGS) $(WARN) -c $<

clean:
	rm -f *.o ./dSFMT/*.o *~ ./dSFMT/*~

remove:
	rm -f main.out *.o ./dSFMT/*.o *~ ./dSFMT/*~
