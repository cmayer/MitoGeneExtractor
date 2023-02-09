# CFLAGS  = -g3 -ggdb -Wall -std=c++11
CFLAGS = -O2 -std=c++11 # -Wall

INCL =    -I . -I tclap-1.2.5/include
#LIBS    = -lc -Wall

SRC    = MitoGeneExtractor.cpp \
	 global-types-and-parameters_MitoGeneExtractor.cpp

HEADER = CDnaString2.h CSequence_Mol2_1.h CSequences2.h CSplit2.h Ctriple.h \
         basic-DNA-RNA-AA-routines.h fast-realloc-vector.h faststring2.h \
         global-types-and-parameters_MitoGeneExtractor.h primefactors.h statistic_functions.h


all:    MitoGeneExtractor-v1.9.5


MitoGeneExtractor-v1.9.5: $(SRC) $(HEADER)
	g++ $(CFLAGS) $(INCL) $(SRC) -o MitoGeneExtractor-v1.9.5


clean:
	rm -f MitoGeneExtractor-v1.9.5


