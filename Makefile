CC = mpic++
BDIR = /opt/local/boost_1_48_0-gcc
NCDFDIR = /opt/local/netcdf-4.1.3-gcc
LIBS = -I${BDIR}/include/ -L${BDIR}/lib/ -lboost_mpi -lboost_serialization

ifeq ($(NETCDF),yes)
NCDFDIR = /opt/local/netcdf-4.1.3-gcc
LIBS += -I$(NCDFDIR)/include/ -L$(NCDFDIR)/lib/ -lnetcdf -lnetcdf_c++
NCDF_FLAGS = -DUSE_NETCDF
endif

FLAGS = -O3 -m64 $(NCDF_FLAGS)
DFLAGS = -g -m64 -D_DEBUG
EXE=MpiScvt.x
TRISRC=Triangle/
XMLSRC=Pugixml/

PLATFORM=_MACOS
PLATFORM=_LINUX

ifeq ($(PLATFORM),_LINUX)
	FLAGS = -O3 -m64 -DLINUX $(NCDF_FLAGS)
	DFLAGS = -g -m64 -D_DEBUG -DLINUX
endif

ifeq ($(PLATFORM),_MACOS)
	FLAGS = -O3 -m64 $(NCDF_FLAGS)
	DFLAGS = -g -m64 -D_DEBUG
endif

TRILIBDEFS= -DTRILIBRARY

all: trilibrary pugixml-library
	${CC} scvt-mpi.cpp ${TRISRC}triangle.o ${XMLSRC}pugixml.o ${LIBS} ${FLAGS} -o ${EXE}

debug: trilibrary-debug pugixml-library-debug
	${CC} scvt-mpi.cpp ${TRISRC}triangle.o ${XMLSRC}pugixml.o ${LIBS} ${DFLAGS} -o ${EXE}

trilibrary:
	$(CC) $(CSWITCHES) $(TRILIBDEFS) ${FLAGS} -c -o ${TRISRC}triangle.o ${TRISRC}triangle.c

trilibrary-debug:
	$(CC) $(CSWITCHES) $(TRILIBDEFS) ${DFLAGS} -c -o ${TRISRC}triangle.o ${TRISRC}triangle.c

pugixml-library:
	$(CC) $(CSWITCHES) $(FLAGS) -c -o $(XMLSRC)pugixml.o $(XMLSRC)pugixml.cpp

pugixml-library-debug:
	$(CC) $(CSWITCHES) $(DFLAGS) -c -o $(XMLSRC)pugixml.o $(XMLSRC)pugixml.cpp

clean:
	rm -f *.dat ${EXE} ${TRISRC}triangle.o $(XMLSRC)/pugixml.o
