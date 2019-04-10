# Compiler choice:

# Gcc
CC = gcc
CFLAGS = -O2 -fopenmp

# IBM p690
#CC = xlc_r
#CFLAGS = -O2 -q64 -qsmp=omp -qreport -fopenmp

MLIBS	=   -lm

QLIB = -L/users/bfalck/qhull-2012.1/lib
QINC = -I/users/bfalck/qhull-2012.1/src/libqhull/

###############

all: vozinit voz1b1 voztie origamitag origamihalo origamicen origamienv

origami: origamitag origamihalo origamicen origamienv

voboz: vozinit voz1b1 voztie

origamitag: origamitag.o readfiles.o
	$(CC) $(CFLAGS) -o origamitag origamitag.o readfiles.o
origamitag.o: origamitag.c
	$(CC) $(CFLAGS) -c origamitag.c

origamihalo: origamihalo.o readfiles.o
	$(CC) $(CFLAGS) -o origamihalo origamihalo.o readfiles.o
origamihalo.o: origamihalo.c
	$(CC) $(CFLAGS) -c origamihalo.c

origamienv: origamienv.o readfiles.o
	$(CC) $(CFLAGS) -o origamienv origamienv.o readfiles.o
origamienv.o: origamienv.c
	$(CC) $(CFLAGS) -c origamienv.c

origamicen: origamicen.o readfiles.o
	$(CC) $(CFLAGS) -o origamicen origamicen.o readfiles.o
origamicen.o: origamicen.c
	$(CC) $(CFLAGS) -c origamicen.c

readfiles.o: readfiles.c
	$(CC) $(CFLAGS) -c readfiles.c

voz1b1: voz1b1.o readfiles.o vozutil.o voz.h
	$(CC)  -o voz1b1 $(CFLAGS) voz1b1.o readfiles.o vozutil.o -L. $(QLIB) -lqhull $(MLIBS) 
voz1b1.o: voz1b1.c
	$(CC) $(CFLAGS) $(QINC) -c voz1b1.c

vozutil.o: vozutil.c
	$(CC) $(CFLAGS) $(QINC) -c vozutil.c

vozinit: vozinit.o readfiles.o voz.h
	$(CC)  -o vozinit $(CFLAGS) vozinit.o readfiles.o -L. $(MLIBS) 
vozinit.o: vozinit.c
	$(CC) $(CFLAGS) -c vozinit.c $(QINC)

voztie: voztie.o readfiles.o
	$(CC)  -o voztie $(CFLAGS) voztie.o readfiles.o
voztie.o: voztie.c
	$(CC) $(CFLAGS) -c voztie.c

clean:
	rm -f *.o

