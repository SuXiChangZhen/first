# Makefile for compiling with personal F90 and
# link with MPICH and FFTW libraries
# To compile and link write simply "make"
#
#module load fftw/2.1.5-double
SHELL =/bin/bash
FC=mpif90
MPIDIR=/nvme1/mpich-pgi
#FFLAG = -O5 -w -r8 -bnoquiet -bloadmap -mcmodel=large -shared-intel
FFLAG =  -r8 -O4 -Minfo=accel -ta=tesla:ccall -DMPI_DEBUG
IDIR =  -I/nvme1/fftw3-pgi/include
LDIR =  -L/nvme1/fftw3-pgi/lib
#IDIR =   -I/PARA/blsc070/software/fftw_2.1.5_gnu/include
#LDIR =   -L/PARA/blsc070/software/fftw_2.1.5_gnu/lib
#FCOMP = mpif90 -g -qfree=f90 -c -qmaxmem=-1 ${FFLAG}
#FCOMP = ${FC} -g -qfree=f90 -c -qmaxmem=-1 ${FFLAG}
FCOMP = ${FC} -c ${FFLAG}
#LINK = mpif90 -g -mcmodel=large -shared-intel 
LINK = ${FC} -g -mcmodel=medium 
LIBS = -lfftw3 -lm 

OBJ = var_inc.o main.o mpifft3.o para.o interface.o feqEXP3thE3Q27A7up.o evol.o gaussian.o initial1.o saveload.o

.SUFFIXES: .o .f90

.f90.o:
	${FCOMP} $*.f90 ${IDIR} ${LDIR} ${LIBS}  
main: ${OBJ}
	${LINK} ${LDIR} ${IDIR} ${FFLAG} ${OBJ} ${LIBS} -o main

clean:
	-rm -rf *.o *.mod main core* output *.err *.out
	-rm -rf coredir*

