OBJS = genutil.o inputparams.o readkey.o keys.o indexarrays.o waveletutil.o covarutil.o

#OBJS = waveletutil.o covarutil.o genutil.o inputparams.o readkey.o keys.o msdutil.o covaranqutil.o

FC = gfortran
FLAGS = -O3 
DEBUGFLAGS =
LIBS = -lblas -llapack
#LIBS = -L/usr/local/lib/atlas/lib -llapack -lf77blas -lcblas -latlas

# executable name
EXE = ../waveletCov.exe
TESTEXE = ../testrun.exe
TESTCOMPEXE = ../testcomprun.exe

.SUFFIXES: .o .f90 

# general suffix rules
%.o: %.f90
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

default: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(EXE) $(LIBS)

test: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o -o $(TESTEXE) $(LIBS)

testcomp: $(OBJS) testcomp.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testcomp.o -o $(TESTCOMPEXE) $(LIBS)

clean: 
	rm *.o *.mod

veryclean: clean
	rm *~ $(EXE)


# ------------
# individual file dependencies
# ------------
#covarutil.o: waveletutil.o genutil.o keys.o
readkey.o: keys.o inputparams.o
waveletutil.o: indexarrays.o
covarutil.o: indexarrays.o waveletutil.o
