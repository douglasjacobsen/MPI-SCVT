CC = mpic++
BDIR = 
LIBS = -L${BDIR}/lib/ -lboost_mpi -lboost_serialization

ifeq ($(NETCDF),yes)
NCDFDIR = 
LIBS += -L$(NCDFDIR)/lib/ -lnetcdf -lnetcdf_c++
NCDF_FLAGS = -DUSE_NETCDF -I$(NCDFDIR)/include
endif

FLAGS = -O3 -m64 $(NCDF_FLAGS)
DFLAGS = -g -m64 -D_DEBUG
EXE=MpiScvt.x
TRISRC=Triangle/

PLATFORM=_MACOS
PLATFORM=_LINUX

ifeq ($(PLATFORM),_LINUX)
	FLAGS = -O3 -m64 -DLINUX $(NCDF_FLAGS) -I$(BDIR)/include/boost
	DFLAGS = -g -m64 -D_DEBUG -DLINUX
endif

ifeq ($(PLATFORM),_MACOS)
	FLAGS = -O3 -m64 $(NCDF_FLAGS) -I$(BDID)/include/boost
	DFLAGS = -g -m64 -D_DEBUG
endif

TRILIBDEFS= -DTRILIBRARY

all: trilibrary
	${CC} scvt-mpi.cpp ${TRISRC}triangle.o ${FLAGS} -o ${EXE} ${LIBS} 

debug: trilibrary-debug
	${CC} scvt-mpi.cpp ${TRISRC}triangle.o ${DFLAGS} -o ${EXE} ${LIBS} 

trilibrary:
	$(CC) $(CSWITCHES) $(TRILIBDEFS) ${FLAGS} -c -o ${TRISRC}triangle.o ${TRISRC}triangle.c

trilibrary-debug:
	$(CC) $(CSWITCHES) $(TRILIBDEFS) ${DFLAGS} -c -o ${TRISRC}triangle.o ${TRISRC}triangle.c

clean:
	rm -f *.dat ${EXE} ${TRISRC}triangle.o
