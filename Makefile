SRC = B0inhom.c OCroutines.c allocation.c auxmath.c averaging.c blockdiag.c \
  cm.c complx.c cryst.c crystdat.c fft.c fidcalc.c ftools.c ham.c iodata.c isotopes.c lbfgs.c \
  main.c matrix.c pthread_barrier_mac.c pulse.c readsys.c relax.c rfprof.c rfshapes.c \
  sim.c simpson.c spinach.c spinsys.c tclcode.c tclutil.c wigner.c distortions.c
OBJ = $(SRC:.c=.o)

# Linux Metacentrum tarkil:
INCLUDES = -I$(MKLROOT)/include/fftw -I/storage/praha1/home/tosner/my_libs/nfft_intelcdk/include -I$(MKLROOT)/include
LIBRARIES = -L/cvmfs/software.metacentrum.cz/spack1/software/tcl/linux-debian9-x86_64/8.6.8-intel-fijli3/lib -ltcl8.6 -L/storage/praha1/home/tosner/my_libs/nfft_intelcdk/lib -lnfft3 -lpthread  -lm 
FLAGS = -DINTEL_MKL -DNO_CONST -O3

CC = icc -mkl
RM = rm
TAR = tar

simpson: $(OBJ)
	$(CC) $(FLAGS) $(OBJ) $(LIBRARIES) -o simpson_intelcdk 
.c.o:
	$(CC) $(FLAGS) $(INCLUDES) -c $<
clean:
	$(RM) -f *.o simpson
dist:
	$(TAR) cvzf simpson.tgz *.c *.h simpson.xcodeproj Makefile
