# ----------------------------------------------- #
# ** Makefile for FMR code **
# ----------------------------------------------- #


# ----------------------------------------------- #
# ** User defined variables **
# ----------------------------------------------- #
# MPI C++ compiler
#CPP = mpicxx
CPP = mpicc-openmpi-mp 
#CPP = CC

# Debugging flags
#DEBUG = -g -DFMR_DEBUG
DEBUG = -g 

# Compiler flags
CFLAGS = -O3 $(DEBUG) -DMPICH_IGNORE_CXX_SEEK

# Libraries, if any
LIBDIR = -L/usr/lib
LIB = -lm -llapack

# ----------------------------------------------- #
# Compilation instructions
# ----------------------------------------------- #

SRC = $(wildcard *.cpp)
INC = $(wildcard *.h)
OBJ = $(SRC:.cpp=.o)

install: $(OBJ)
	$(CPP) $(OBJ) $(LIB) -o fmr.exe

clean: 
	rm -f *.o fmr.exe

%.o:%.cpp
	$(CPP) $(CFLAGS) -c $<
