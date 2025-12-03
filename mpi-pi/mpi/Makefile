CC = mpicc

CFLAGS = -I. -std=c11 -O2
LIBS = -lm -lmpi

DEPS = proto.h
OBJ  = main.o calculate_pi.o
EXE = calcpi

.PHONY: all clean 

all: $(EXE)

$(EXE): $(OBJ) $(DEPS)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LIBS)

# make sure accuracy.o is rebuilt if proto.h changes


clean:
	-rm -f $(OBJ) $(EXE) 
