#
# Makefile for SPME_NlogN
#
# Byoungseon Jeon
# Dept. Applied Science, UC.Davis
# July 31. 2005
#
# suffix definition for compiling
.SUFFIXES: .o .f90
#
#
# select compiler
F90 = /opt/openmpi.intel/bin/mpif90
FLAGS = -g #-O3 -fast -ipo -no-prec-div #-Wall -g -O3 # -march=native -ffast-math -funroll-loops -O3 # -g -Wall
# Object files

# Object files
OBJTS = datastr.o main.o force.o vverlet.o post.o potential.o pme.o cell.o\
	fft235.o kernel.o mfft235.o pzfft3dv.o
TARGET = cascade_morelon
#
# generation of executable
${TARGET}:${OBJTS}
	${F90} ${FLAGS} -o ${TARGET} ${OBJTS} ${LIB} ${INCLUDE}
#
# generation of object files
.f90.o:
	${F90} ${FLAGS} -c $< ${LIB}
.f.o:
	${F90} ${FLAGS} -c $< ${LIB}
#
# clean object files and executable
clean:
	rm -rf *.mod *.o *.f90~ core ${TARGET}
