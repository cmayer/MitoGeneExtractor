# CFLAGS  = -g3 -ggdb -Wall -std=c++11
CFLAGS = -O2 -std=c++11 # -Wall

INCL =    -I . -I tclap-1.2.5/include
#LIBS    = -lc -Wall

SRC    = MitoGeneExtractor.cpp \
	 global-types-and-parameters_MitoGeneExtractor.cpp \
         exonerate_wrapper_and_parser.cpp

HEADER = CDnaString3.h CSequence_Mol3.h CSequences3.h CSplit2.h Ctriple.h \
         basic-DNA-RNA-AA-routines.h fast-realloc-vector.h faststring3.h \
         global-types-and-parameters_MitoGeneExtractor.h primefactors.h statistic_functions.h \
         exonerate_wrapper_and_parser.hpp fastq.h


all:    MitoGeneExtractor-v1.9.6


MitoGeneExtractor-v1.9.6: $(SRC) $(HEADER)
	g++ $(CFLAGS) $(INCL) $(SRC) -o MitoGeneExtractor-v1.9.6


clean:
	rm -f MitoGeneExtractor-v1.9.6


