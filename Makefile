##
## Makefile for coseg project
##
VERSION=0.2.2
INSTALLDIR=/usr/local/coseg-${VERSION}

## Basic
CC = cc -O4 -lm
## A nice memory leak checker:
#CC = bgcc -O4 -fbounds-checking -lm

all: coseg 

version.c: Makefile
	echo "char const* Version = \"$(VERSION)\";" > version.c

coseg: version.o coseg.o
	${CC} version.o coseg.o -o coseg 

beautify:
	indent -bap -cdb -bl -bli0 -npcs -nut -lp coseg.c

install: coseg
	-mkdir ${INSTALLDIR}
	cp coseg ${INSTALLDIR}
	cp README ${INSTALLDIR}
	cp preprocessAlignments.pl ${INSTALLDIR}
	cp postprocess.pl ${INSTALLDIR}
	cp runcoseg.pl ${INSTALLDIR}
	cp extractSubSeqs.pl ${INSTALLDIR}
	cp refineConsSeqs.pl ${INSTALLDIR}
	cp ALU.cons ${INSTALLDIR}
	cp ALU.ins ${INSTALLDIR}
	cp ALU.seqs ${INSTALLDIR}
	cp LINE1 ${INSTALLDIR}


dist:
	-mkdir dist
	-mkdir dist/coseg
	cp ALU.cons dist/coseg
	cp ALU.ins dist/coseg
	cp ALU.seqs dist/coseg
	cp LINE1 dist/coseg
	cp Makefile dist/coseg
	cp README dist/coseg
	cp coseg.c dist/coseg
	cp coseg.h dist/coseg
	cp postprocess.pl dist/coseg
	cp preprocessAlignments.pl dist/coseg
	cp runcoseg.pl dist/coseg
	cp extractSubSeqs.pl dist/coseg
	cp refineConsSeqs.pl dist/coseg
	(cd dist; tar zcvf coseg-$(VERSION).tar.gz coseg)

clean:
	-rm *.o
	-rm coseg
	-rm version.c
	-rm ALU.seqs.subfamilies.seq
	-rm ALU.seqs.assign
	-rm ALU.seqs.log
	-rm ALU.seqs.subfamilies
	-rm ALU.seqs.tree.viz


ALU.seqs.subfamilies.seq: coseg
	./runcoseg.pl -d -filePrefix ALU

t/kothi.seqs.subfamilies.seq: coseg
	./runcoseg.pl -u1 -t -m 5 -filePrefix t/kothi

test: ALU.seqs.subfamilies.seq t/kothi.seqs.subfamilies.seq
	diff ALU.seqs.subfamilies.seq t
	diff t/kothi.seqs.subfamilies.seq t/kothi-baseline
