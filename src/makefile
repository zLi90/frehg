# makefile for the implicit advdiff solver
HOME=laspack
#include $(HOME)

CC = mpicc
#CC = mpicc -O3 -fp-model precise -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512
# all:
# 	$(CC) bathymetry.c configuration.c fileio.c groundwater.c initialize.c map.c \
# 		mpifunctions.c nsfunctions.c nssolve.c scalar.c subgrid.c utilities.c \
# 	  $(HOME)/eigenval.c $(HOME)/errhandl.c $(HOME)/factor.c \
# 		$(HOME)/itersolv.c $(HOME)/matrix.c $(HOME)/mlsolv.c \
# 		$(HOME)/operats.c $(HOME)/precond.c $(HOME)/qmatrix.c $(HOME)/vector.c \
# 		$(HOME)/rtc.c FREHD.c -O3 -lm -o runthis.o

all:
	$(CC) configuration.c groundwater.c initialize.c map.c mpifunctions.c \
		  scalar.c shallowwater.c solve.c utility.c linsys.c\
		  $(HOME)/eigenval.c $(HOME)/errhandl.c $(HOME)/factor.c $(HOME)/itersolv.c \
		  $(HOME)/matrix.c $(HOME)/mlsolv.c $(HOME)/operats.c $(HOME)/precond.c \
		  $(HOME)/qmatrix.c $(HOME)/vector.c $(HOME)/rtc.c FREHG.c -lm -O3 -fopenmp -o frehg
