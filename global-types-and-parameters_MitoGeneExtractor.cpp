#include "global-types-and-parameters_MitoGeneExtractor.h"
#include "CmdLine.h"
#include "climits"

using namespace TCLAP;
using namespace std;

string                      global_input_dna_fasta_file;
string                      global_input_prot_reference_sequence;
unsigned                    global_verbosity;
unsigned                    global_num_bp_beyond_exonerate_alignment_if_at_start_or_end;
string                      global_exonerate_binary;
string                      global_vulgar_file;
string                      global_alignment_output_file;
float                       global_consensus_threshold;
string                      global_consensus_sequence_output_filename;
char                        global_include_frameshift_alignments;
char                        global_include_gap_alignments;
char                        global_include_double_hits;
float                       global_relative_score_threshold;
// Not implemented as parameters yet:
char                        global_genetic_code_number;
int                         global_frameshift_penalty;
char                        global_run_mode;  // default: normal run, 'd' exonerate dry-run,
unsigned                    global_minimum_seq_coverage_uppercase;
unsigned                    global_minimum_seq_coverage_total;
unsigned                    global_exonerate_score_threshold;

void good_bye_and_exit(int error)
{
  if (error == 0)
    cout << "## Finished successfully. ##" << endl;
  cerr << endl << PROGNAME << " says goodbye." << endl << endl;
  exit(error);
}


void init_param()
{
  global_num_bp_beyond_exonerate_alignment_if_at_start_or_end = 0;
  global_consensus_threshold = 0.5;
  global_verbosity = 1;
  global_include_gap_alignments = 0;
  global_include_frameshift_alignments = 0;
  global_include_double_hits  = 0;
  global_relative_score_threshold = 0;
  global_genetic_code_number = 2;
  global_frameshift_penalty  = -9;
  global_run_mode = 0;
  global_minimum_seq_coverage_uppercase = 1;
  global_minimum_seq_coverage_total = 1;
}


void read_and_init_parameters(int argc, char** argv)
{
  init_param();

  bool option_set_consensus_sequence_output_filename;
  bool option_set_consensus_threshold;

  try
  {
    CmdLine cmd("Description: This program aligns multiple nucleotide sequences against a protein reference "
		"sequence to obtain a multiple sequence alignment. The recommended use case is to extract mitochondrial genes from sequencing libraries, "
                "e.g. from hybrid enrichment libraries sequenced on the Illumina platform."
		"The individual alignments are computed by calling the exonerate program. "
		"Exonerate is a very efficient alignment program which allows to align protein and nucleotide sequences."
		"Nucleotide sequences which cannot be aligned to the protein reference will not be included in the output. "
                "Exonerate should be able to align 100k short reads in a few minutes using a single CPU core, so this approach can be used for projects of any size."
		,
		' ', VERSION);
    /*
    ValueArg<string> path_exonerate_binary_Arg("E", "exonerateBinary",
	"Path and filename or only the filename of the exonerate program. Alternatively, to letting "
	PROGNAME " call exonerate, the user can provide an existing vulgar output file. "
        "If " PROGNAME " has to call exonerate, which is the case if specified vulgar file does not exist, "
	"a Protein reference sequence has to be specified. If no executable is specified with this command, "
       PROGNAME " tries to run the \"exonerate\" command.", false,
       global_exonerate_binary.c_str(), "string");
    cmd.add(path_exonerate_binary_Arg);
    */

    ValueArg<unsigned> verbosity_Arg("", "verbosity",
	"Specifies how much run time information is printed to the console. Values: 0: minimal output, 1: important notices, 2: more notices, 3: basic progress, 4: detailed progress, 50-100: debug output, 1000: all output.",
	false, global_verbosity, "int");
    cmd.add(verbosity_Arg);

    ValueArg<unsigned> min_exonerate_score_threshold_Arg("s", "minExonerateScoreThreshold",
						       "The score threshold passed to exonerate to decide whether to include or not include the hit in the output.",
						       false, UINT_MAX, "int");
    cmd.add(min_exonerate_score_threshold_Arg);

    ValueArg<unsigned> min_seq_coverage_upper_case_Arg("", "minSeqCoverageInAlignment_uppercase",
						       "Specifies the absolute value of the minimum alignment coverage for computing the consensus sequence. As coverage, only upper case nucleotides are taken into account, i.e. no nucleotides are counted that have been added beyond the exonerate alignment region. Bases beyond the exonerate alignment are added with the -n or --numberOfBpBeyond option. If no bases are added beyond the exonerate alignment (default), the effect of this option is identical to the minSeqCoverageInAlignment_total option. Default: 1. Increasing this value increases the number of unknown nucleotides in the consensus sequence.",
						       false, global_minimum_seq_coverage_uppercase, "int");
    cmd.add(min_seq_coverage_upper_case_Arg);

    ValueArg<unsigned> min_seq_coverage_total_Arg("", "minSeqCoverageInAlignment_total",
						       "Specifies the absolute value of the minimum alignment coverage for computing the consensus sequence. For the coverage, all nucleotides count, also those lower case nucleotides that have been added beyond the exonerate alingmnet region. Default: 1. Increasing this value increases the number of unknown nucleotides in the consensus sequence.",
						       false, global_minimum_seq_coverage_total, "int");
    cmd.add(min_seq_coverage_total_Arg);
    
    
    ValueArg<float> relative_score_threshold_Arg("r", "relative_score_threshold",
	"Specified the relative alignment score threshold for exonerate hits to be considered. "
	"The relative score is the score reported by exonerate divided by the alignment length. Default 0. "
	"Reasonable thresholds are between 1 and 2.",
	false, global_relative_score_threshold, "float");
    cmd.add(relative_score_threshold_Arg);

    ValueArg<int> genetic_code_number_Arg("C", "genetic_code",
	 "The number of the genetic code to use in exonerate, if this step is required. Default: 2, i.e. "
	"vertebrate mitochondrial code.",
	false, global_genetic_code_number, "int");
    cmd.add(genetic_code_number_Arg);

    ValueArg<int> frameshift_penalty_Arg("f", "frameshift_penalty",
	 "The frameshift penalty passed to exonerate. The option value has to be a negative integer. Default: -9. Higher values lead to lower scores and by this can have the following effects: (i) hit regions are trimmed since trimming can lead to a better final alignment score, (ii) they can also lead to excluding a read as a whole if the final score is too low and trimming does lead to a higher score. The default of the exonerate program is -28. A value of -9 (or other values lower than -28) lead to more reads in which the best alignment has a frameshift. In order to remove reads that do not align well, one can use a smaller frameshift penalty and then exclude hits with a frameshift, see -F option).",
	false, global_frameshift_penalty, "int");
    cmd.add(frameshift_penalty_Arg);

    SwitchArg include_F_Arg("F", "includeFrameshift",
			 "Include reads which aligned with a frameshift. Default: false.",
			 false);
    cmd.add(include_F_Arg);

    SwitchArg include_G_Arg("G", "includeGap",
			 "Include reads which aligned with a gap.",
			 false);
    cmd.add(include_G_Arg);

    SwitchArg include_D_Arg("D", "includeDoubleHits",
			 "Include reads with two alignment results found by exonerate.",
			 false);
    cmd.add(include_D_Arg);

    ValueArg<float> consensus_threshold_Arg("t", "consensus_threshold",
	"This option modifies the consensus threshold. Default: 0.7 which corresponds to 70%.",
	false, global_consensus_threshold, "float");
    cmd.add(consensus_threshold_Arg);

    ValueArg<string> consensus_output_file_Arg("c", "consensus_file",
	"If this option is specified, a consensus sequence of all aligned reads is written to the file with the specified name. Normally, this is the intended output. "
	"Default: No consensus is written, since no good default output file is known.",
	false, global_consensus_sequence_output_filename, "string");
    cmd.add(consensus_output_file_Arg);

    ValueArg<unsigned> bp_beyond_Arg("n", "numberOfBpBeyond",
	"Specifies the number of base pairs that are shown beyond the exonerate alignment. A value of 0 means that the sequence is clipped at the point the exonerate alignment ends. Values of 1 and 2 make sense, since exonerate does not consider partial matches of the DNA to the amino acid sequence, so that partial codons would always be clipped, even if the additional base pairs would match with the expected amino acid. Values >0 lead to the inclusion of sequence segments that do not align well with the amino acid sequence and have to be treated with caution. They might belong to chimera, numts, or other problematic sequences. Larger values might be included e.g. if problematic sequences with a well matching seed alignment are of interest. CAUTION: Bases included with this option might not be aligned well or could even belong to stop codons! They should be considered as of lower quality compared to other bases. Bases that are added with this option are added as lower case characters to the output alignment file. A sequence coverage of bases not belonging to these extra bases can be requested with the --minSeqCoverageInAlignment_uppercase option. Default: 0.",
	false, global_num_bp_beyond_exonerate_alignment_if_at_start_or_end, "int");
    cmd.add(bp_beyond_Arg);

    ValueArg<string> exonerate_path_Arg("e", "exonerate_program",
     				        "Name of the exonerate program in system path OR the path to the exonerate program including the program name. Default: exonerate",
				        false, global_exonerate_binary, "string");
    cmd.add(exonerate_path_Arg);

    ValueArg<string> vulgar_file_Arg("V", "vulgar_file",
				     "Name of exonerate vulgar file. If the specified file exists, it will be used for the analysis. "
				     "If it does not exist " PROGNAME " will run exonerate in order to create the file. The created file will then be used to proceed. If no file is specified with this option, a temporary file called tmp-vulgar.txt will be created and removed after the program run. In this case a warning will be printed to the console.",
				     false, global_vulgar_file, "string");
    cmd.add(vulgar_file_Arg);

    ValueArg<string> alignment_output_file_name_Arg("o", "",
	"Name of alignment output file. ",
	true, global_alignment_output_file, "string");
    cmd.add(alignment_output_file_name_Arg);

    ValueArg<string> prot_fasta_input_file_name_Arg("p", "prot_reference_file",
	"Protein sequence file in the fasta format. This is the sequence used to align the reads against. File is expected to have exactly one reference sequence. ",
	false, global_input_prot_reference_sequence, "string");
    cmd.add(prot_fasta_input_file_name_Arg);

    ValueArg<string> dna_fasta_input_file_name_Arg("d", "dna_sequences_file",
	"Nucleotide sequence file in the fasta format. Sequences are expected to be unaligned. ",
	true, global_input_dna_fasta_file, "string");
    cmd.add(dna_fasta_input_file_name_Arg);

    cmd.parse( argc, argv );

    // Assigning parameters to variables:
    global_input_dna_fasta_file                        = dna_fasta_input_file_name_Arg.getValue();
    global_input_prot_reference_sequence               = prot_fasta_input_file_name_Arg.getValue();
    global_verbosity                                   = verbosity_Arg.getValue();
    global_num_bp_beyond_exonerate_alignment_if_at_start_or_end = bp_beyond_Arg.getValue();
    global_exonerate_binary                            = exonerate_path_Arg.getValue();
    global_vulgar_file                                 = vulgar_file_Arg.getValue();
    global_alignment_output_file                       = alignment_output_file_name_Arg.getValue();
    global_consensus_sequence_output_filename          = consensus_output_file_Arg.getValue();
    global_consensus_threshold                         = consensus_threshold_Arg.getValue();

    global_include_double_hits                         = include_D_Arg.getValue();
    global_include_frameshift_alignments               = include_F_Arg.getValue();
    global_include_gap_alignments                      = include_G_Arg.getValue();
    global_relative_score_threshold                    = relative_score_threshold_Arg.getValue();

    option_set_consensus_sequence_output_filename      = consensus_output_file_Arg.isSet();
    option_set_consensus_threshold                     = consensus_threshold_Arg.isSet();

    global_genetic_code_number                         = genetic_code_number_Arg.getValue();
    global_frameshift_penalty                          = frameshift_penalty_Arg.getValue();
    global_minimum_seq_coverage_total                  = min_seq_coverage_total_Arg.getValue();
    global_minimum_seq_coverage_uppercase              = min_seq_coverage_upper_case_Arg.getValue();

    global_exonerate_score_threshold                   = min_exonerate_score_threshold_Arg.getValue();
  }
  catch (ArgException &e)
  {
    cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
    good_bye_and_exit(-1);
  }

  if (global_vulgar_file.empty())
  {
    cerr << "WARNING: You did not specify a vulgar file, so a temporary vulgar file will be created for this run that will be removed at the end of this program run. Therefor the vulgar file cannot be reused in other runs.\n"; 
    global_vulgar_file = "tmp-vulgar.txt";
  }

  if (global_consensus_threshold > 1 || global_consensus_threshold < 0)
  {
    cerr << "ERROR: You specified a consensus threshold of " << global_consensus_threshold << ". Allowed values are in the range 0..1, "
      "where 1 corresponds to a strict consensus where 100% identity is required and 0 implies that the dominant nucleotide "
      "is chosen independent of its proportion." << endl;
    cerr << "Please choose a threshold in the range 0..1 and restart this program. Exiting." << endl;
    exit(-1);
  }
  
  if (option_set_consensus_threshold && !option_set_consensus_sequence_output_filename ) 
  {
    cerr << "You specified a consensus threshold but you do not request to write a consensus sequence file. In order to request to write a consensus sequence file, you have to specify the -c or --consensus_file option followed by a file name";
    exit(-1);
  }
}



void print_parameters(std::ostream &os, const char *s)
{
  os << s <<   "Parameter settings:"
     << std::endl
     << s <<   "===================" 
     << std::endl;
  
  os << s <<   "DNA input file name:                                    " << global_input_dna_fasta_file
     << std::endl;

  os << s <<   "Protein reference input file name:                      " << global_input_prot_reference_sequence
     << std::endl;

  os << s <<   "Output file name:                                       " << global_alignment_output_file
     << std::endl;

  if (!global_vulgar_file.empty())
    os << s << "Vulgar file name:                                       " << global_vulgar_file
       << std::endl;

  if (!global_exonerate_binary.empty()) 
    os << s << "Exonerate binary:                                       " << global_exonerate_binary
       << std::endl;

  os << s <<   "Print this number of bp beyond exonerate alignment:     " << global_num_bp_beyond_exonerate_alignment_if_at_start_or_end
     << std::endl;

  if (!global_consensus_sequence_output_filename.empty())
  {
    os << s << "Write consensus sequence to file :                      " << "Yes"
       << std::endl;
    os << s << "Filename for consensus sequence output:                 " << global_consensus_sequence_output_filename
       << std::endl;
    os << s << "Consensus sequence threshold value:                     " << global_consensus_threshold
       << std::endl;
  }
  else
  {
    os << s << "Write consensus sequences to file:                      " << "No"
       << std::endl;
  }

  os << s <<   "Frameshift penalty:                                     " << global_frameshift_penalty
     << std::endl;

  os << s <<   "Relative score threshold:                               " << global_relative_score_threshold
     << std::endl;

  os << s <<   "Minimum coverage in exonerate alignment:                " << global_minimum_seq_coverage_total
     << std::endl;

  os << s << "Minimum coverage in exonerate alignment (upper case):   " << global_minimum_seq_coverage_uppercase
     << std::endl;
  
  if (global_exonerate_score_threshold != UINT_MAX)
    os << s << "Exonerate score threshold:                              " << global_exonerate_score_threshold
       << std::endl;    

  os << s <<   "Verbosity:                                              " << global_verbosity
     << std::endl;
  os << s <<   "===================" 
     << std::endl;
  


}


void print_output_creation_message(std::ostream &os, char *s)
{
  os << s << "Output created by " << PROGNAME << ", version " << VERSION
     << std::endl;
  os << s << std::endl;
  print_parameters(os, s);
  os << s << std::endl;
}


