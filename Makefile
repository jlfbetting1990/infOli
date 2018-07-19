# Example Makefile

#MACROS
CAT_HOME = ${MGC_HOME}
TARGET = InferiorOlive.x
OBJECTS = infoli_tb.o infoliNew.o infoliOld.o
DEPENDS = infoli.h infoliNew.h
INCLUDES = -I"${CAT_HOME}/shared/include"
DEFINES = 
CXX = /usr/bin/g++
CXXFLAGS = -g -o3 ${DEFINES} ${INCLUDES}

#my_tb target is dependent on main.o and hello.o

${TARGET} : ${OBJECTS}
	${CXX} ${CXXFLAGS} -o ${TARGET} ${OBJECTS}


#phony target to remove all objects and execs
.PHONY: clean
clean:
	rm -f *.o *.exe *.x
