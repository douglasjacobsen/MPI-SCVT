CC = mpic++
LIBS = -lboost_mpi -lboost_serialization
FLAGS = -O3 -m64
DFLAGS = -g -m64 -D_DEBUG
EXE=MpiScvt.x
TRISRC=Triangle/

TRILIBDEFS= -DTRILIBRARY

all: trilibrary
	${CC} scvt-mpi.cpp ${TRISRC}triangle.o ${LIBS} ${FLAGS} -o ${EXE}

debug: trilibrary-debug
	${CC} scvt-mpi.cpp ${TRISRC}triangle.o ${LIBS} ${DFLAGS} -o ${EXE}

trilibrary:
	$(CC) $(CSWITCHES) $(TRILIBDEFS) ${FLAGS} -c -o ${TRISRC}triangle.o ${TRISRC}triangle.c

trilibrary-debug:
	$(CC) $(CSWITCHES) $(TRILIBDEFS) ${DFLAGS} -c -o ${TRISRC}triangle.o ${TRISRC}triangle.c

clean:
	rm -f *.dat ${EXE} ${TRISRC}triangle.o