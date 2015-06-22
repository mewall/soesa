SOESA_ROOT = <Insert SOESA root dir here>
VPATH=${SOESA_ROOT}/src:${SOESA_ROOT}/bin
S=${SOESA_ROOT}/src
B=${SOESA_ROOT}/bin
I=${SOESA_ROOT}/src
CCFLGS= -Ofast -c -I${I}
CCEXECFLAGS= -Ofast -I${I}
CCLIBS= -lm
CC=cc

all: soesa tellsoesa
soesa:soesa.c pdftool.o
	${CC} ${CCEXECFLAGS} -o ${B}/soesa ${S}/soesa.c ${S}/pdftool.o ${CCLIBS}
tellsoesa:tellsoesa.c pdf.h
	${CC} ${CCEXECFLAGS} -o ${B}/tellsoesa ${S}/tellsoesa.c ${CCLIBS}
pdftool.o:pdftool.c pdf.h
	${CC} ${CCFLGS} -o ${S}/pdftool.o ${S}/pdftool.c






