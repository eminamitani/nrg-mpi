## for gfortran
F90 = mpif90 -fopenmp -ggdb
FPPFLAGS = -cpp
FFLAGS= -O2
LFLAGS= -llapack -lblas

## for ifort
#F90 = mpif90 -openmp
#FPPFLAGS = -fpp
#FFLAGS= -O2 -CB -traceback -mcmodel=large
#LFLAGS=  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core


## for pgf90
#F90 = mpif90 -mp
#FPPFLAGS =
#FFLAGS= -O2
#LFLAGS=  -llapack -lblas

TARGET = test

OBJECTS = \
	count.o module.o main.o read_initial_basis.o loadbalance.o \
	invariantMatrix.o nrgIteration.o prepareBasis.o \
	diagonalization.o postprocess.o initializeReducedDensity.o \
	reducedDensityMatrix.o invariantMatrixForSpectrum.o \
	calculateSpectrum.o initialState.o 

MODULE = module.F90

.SUFFIXES: .o .f90 .F90

.F90.o:
	${F90} -c ${FPPFLAGS} ${FFLAGS} $<

.f90.o:
	${F90} -c ${FFLAGS} $<

all:	$(TARGET)

$(TARGET): $(OBJECTS)
	$(F90) -o  $(TARGET) $(OBJECTS) $(LFLAGS)

${OBJECTS}: ${MODULE}			

clean:
	rm -f *.o *.mod $(TARGET)
