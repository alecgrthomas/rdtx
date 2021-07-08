
##################################################################################################
# RDTX Makefile
#
# To compile RDTX on your system, please edit the sections below setting your options for:
#   - C++ compiler
#
##################################################################################################
NPROCS=1 # note not yet parallelized
EXEC = rdtx_beta
SDIR = src
ODIR = ${SDIR}/obj
LDIR = lib
IDIR = include
BDIR = bin
LIBS = 

# CXX - C++ compiler command
# CXXFLAGS - compilation flags
CXX = g++
CXXFLAGS = -Wno-variadic-macros -Wno-long-long -Wall -O3 -funroll-loops -ffast-math -pedantic 
DBGFLAGS = -v -O0

##################################################################################################


##################################################################################################
# MPI 
##################################################################################################

# Not yet parallel

##################################################################################################
################################### Do not edit below this line! ################################# 
##################################################################################################

# Set the CPP object files
CF=-I$(IDIR) -I. ${CXXFLAGS} ${WXFLAGS}
DF=-I$(IDIR) -I. ${DBGFLAGS}

_DEPS = radiation_force.h fourtens.h fourvec.h constants.h input_params.h \
	fresnel.h global_constants.h headers.h namespaces.h main.h  \
	

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = derived_params.o field_funcs.o fields.o particles.o spectralfunctions.o \
		write_data.o functions.o spin.o rdtx.o main.o  \
		 run_params.o import.o\
		 

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all : $(EXEC) 

$(ODIR)/%.o: ${SDIR}/%.cpp $(DEPS) 
	$(CXX) -c -o $@ $< $(CF)

$(EXEC) : $(OBJ)
	${CXX} -o ${BDIR}/$@ $^ $(CF) $(LIBS)



debug : ${OBJ}
	${CXX} -o ${BDIR}/$@ $^ $(DF) $(LIBS)

.PHONY: clean run

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
	rm -rf $(EXEC).app

run:
	@echo Running on ${NPROCS} processors	
	mpiexec -n ${NPROCS} ${BDIR}/${EXEC}
