DIRS = ./
ODIR = ./

# compiler
CC_SERIAL     =  g++
CC_MPI        =  mpicxx
OMPI_CC       =  mpicxx
OMPI_CLINKER  =  mpicxx
MPI_LIB       =  ${MPI_HOME}
OPTFLAGS	  =  -Wall -O3
CPPFLAGS	  =  -std=c++11

# define any directories containing header files other than /usr/include
#   define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
#   define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname
#   option, something like (this will link in libmylib.so and libm.so:

LMP_LIB     =  ${SRC}/lammps/src/
MY_INC 		=  /home/chaomy/src/armadillo-8.300.3/include/
CINCLUDE 	=  -I../include  -I${HOME}/install/include  -I${MPI_HOME}/include  -I${LMP_LIB} -I${MY_INC} -I${BOOST_INCLUDE} 
LIBS   		=  -L${HOME}/install/lib  -L${MPI_HOME}/lib -L${LMP_LIB} -L${BOOST_LIB} 
CDLINK  	=  ${LIBS} -lm -lnlopt -lmpi -lpthread -llammps -DARMA_DONT_USE_WRAPPER -lblas -llapack -lboost_mpi-mt -lboost_serialization-mt -lboost_serialization-mt

MAKETARGET 	=  pfmpi 

# The source files
POTFITSRC 	=  ./pfConfig.cpp     \
			   ./pfForceEAM.cpp   \
			   ./pfForceADP.cpp   \
			   ./pfForceMEAM.cpp  \
			   ./pfStressMEAM.cpp \
			   ./pfInputForce.cpp \
			   ./pfInputPot.cpp   \
			   ./pfInputLMP.cpp   \
			   ./pfOptAnneal.cpp  \
			   ./pfOptEvo.cpp \
			   ./pfOptCMAES.cpp \
			   ./pfRand.cpp   \
			   ./pfSpline.cpp \
			   ./pfUtil.cpp   \
			   ./pfForce.cpp  \
			   ./pfHome.cpp   \
			   ./pfOutPot.cpp  \
			   ./pfOutConfig.cpp \
			   ./pfUpdateRho.cpp \
			   ./pfRescale.cpp \
			   ./pfOptimizer.cpp \
			   ./pfParam.cpp \
			   ./pfUpgrade.cpp \
			   ./pfNeighbors.cpp \
			   ./pfLmpDrv.cpp \
			   ./pfLmpBCC.cpp \
			   ./pfLmpFCC.cpp \
			   ./pfLmpHCP.cpp \
			   ./pfLmpEla.cpp \
			   ./pfLmpSuf.cpp \
			   ./pfLmpGSF.cpp \
			   ./pfLmpPV.cpp \
			   ./pfLmpVac.cpp \
			   ./pfLmpItenOpt.cpp \
			   ./pfLmpItenRun.cpp \
			   ./pfEle.cpp 	\
			   ./pfElePV.cpp \
			   ./pfEleGSF.cpp \
			   ./pfEleDFT.cpp \
			   ./pfAnaly.cpp \
			   ./pfPhyBcc.cpp \
			   ./pfPhyFcc.cpp \
			   ./pfPhyHcp.cpp \
			   ./pfPhyStrain.cpp \
			   ./pfPhyEla.cpp \
			   ./pfPhyTar.cpp \
			   ./pfPhyPV.cpp \
			   ./pfPhySur.cpp \
			   ./pfPhyD03.cpp \
			   ./pfMPI.cpp  \
			   ./pfMain.cpp

			   # ./pfForceMEAMtb.cpp  \

# parallel
PARALLEL = MPI
CC = ${CC_MPI}

# OPT_FLAGS   += ${${PARALLEL}_FLAGS} ${OPT_${PARALLEL}_FLAGS} -DNDEBUG
###########################################################################
# Rules
###########################################################################

# all objects depend on headers
OBJECTS := $(subst .cpp,.o,${POTFITSRC})

%.o: %.cpp
	${CC} -c  $<  ${CINCLUDE}  ${OPTFLAGS}  ${CPPFLAGS}  # ${LIBS}

${MAKETARGET}:$(OBJECTS)

all: ${MAKETARGET}
	${CC}  -o  ${MAKETARGET}  ${OBJECTS}  ${CDLINK}  ${OPTFLAGS}  ${CPPFLAGS}

clean:
	rm -f  *.o  *.u  *~ # \#* *.V *.T *.O *.il
