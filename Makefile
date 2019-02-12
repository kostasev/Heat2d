MPICC = mpicc
FLAGS = -Wall
THREADS=4
SIZE=900
STEPS=10000
STEP=20
SOURCE = mpi_heat2D_opt.c
all: heat_mpi

heat_mpi: $(SOURCE)
	$(MPICC) $(FLAGS) -o heat_mpi $(SOURCE) -lm

clean:
	rm -f $(BINS)
