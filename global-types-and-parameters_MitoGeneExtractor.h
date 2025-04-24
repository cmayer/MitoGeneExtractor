#ifndef GLOBAL_TYPES_AND_PARAMETERS_MITOGENEEXTRACTOR_H
#define GLOBAL_TYPES_AND_PARAMETERS_MITOGENEEXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>
//#include "faststring3.h"

// #define DEBUG

#define PROGNAME "MitoGeneExtractor"
#define VERSION  "1.9.6beta3"

#define macromax(x,y) ((x)<(y) ? (y) : (x))
#define macromin(x,y) ((x)<(y) ? (x) : (y))

#ifdef  DEBUG
#define DEBUGOUT1(x)        std::cerr << x                << std::endl;
#define DEBUGOUT2(x,y)      std::cerr << x << y           << std::endl;
#define DEBUGOUT3(x,y,z)    std::cerr << x << y << z      << std::endl;
#define DEBUGOUT4(x,y,z,w)  std::cerr << x << y << z << w << std::endl;
#else
#define DEBUGOUT1(x)
#define DEBUGOUT2(x,y)
#define DEBUGOUT3(x,y,z)
#define DEBUGOUT4(x,y,z,w)

#endif


void good_bye_and_exit(int);
void init_param();
void read_and_init_parameters(int argc, char** argv);
void print_parameters(std::ostream &os, const char *s);

extern std::vector<std::string>    global_input_dna_fasta_filenames;
extern std::vector<std::string>    global_input_dna_fastq_filenames;
extern std::string                 global_input_prot_reference_sequence;
extern std::string                 global_tmp_directory;
extern unsigned                    global_verbosity;
extern unsigned                    global_num_bp_beyond_exonerate_alignment_if_at_start_or_end;
extern std::string                 global_exonerate_binary;
extern std::string                 global_vulgar_directory;
extern std::string                 global_alignment_output_file;
extern float                       global_consensus_threshold;
extern std::string                 global_consensus_sequence_output_filename;
extern char                        global_include_frameshift_alignments;
extern char                        global_notinclude_gap_alignments;
extern char                        global_include_only_gap_alignments;
extern char                        global_include_double_hits;
extern float                       global_relative_score_threshold;
extern char                        global_genetic_code_number;
extern int                         global_frameshift_penalty;
extern char                        global_run_mode; // Exonerate run mode. Supported modes: 'd': dry run.
extern unsigned                    global_minimum_seq_coverage_uppercase;
extern unsigned                    global_minimum_seq_coverage_total;
extern unsigned                    global_exonerate_score_threshold;
extern int                         global_gap_frameshift_mode;
extern int                         global_report_gap_mode;
extern bool                        global_keep_concatenated_input_file;

extern unsigned                    global_ends_width;
extern unsigned                    global_weight_fraction_in_ends;
extern bool                        global_treat_references_as_individual;
//extern bool                        global_report_internal_gaps_as_tilde;


#endif
