#include "global-types-and-parameters_MitoGeneExtractor.h"
#include "tclap/CmdLine.h"
#include "climits"

#include <sys/types.h>
#include <sys/stat.h>


using namespace TCLAP;
using namespace std;

vector<string>              global_input_dna_fasta_filenames;
vector<string>              global_input_dna_fastq_filenames;
string                      global_input_prot_reference_sequence;
string                      global_tmp_directory;
unsigned                    global_verbosity;
unsigned                    global_num_bp_beyond_exonerate_alignment_if_at_start_or_end;
string                      global_exonerate_binary;
string                      global_vulgar_directory;
string                      global_alignment_output_file;
float                       global_consensus_threshold;
string                      global_consensus_sequence_output_filename;
char                        global_include_frameshift_alignments;
char                        global_notinclude_gap_alignments;
char                        global_include_only_gap_alignments;
char                        global_include_double_hits;
float                       global_relative_score_threshold;
// Not implemented as parameters yet:
char                        global_genetic_code_number;
int                         global_frameshift_penalty;
char                        global_run_mode;  // default: normal run, 'd' exonerate dry-run,
unsigned                    global_minimum_seq_coverage_uppercase;
unsigned                    global_minimum_seq_coverage_total;
unsigned                    global_exonerate_score_threshold;
//bool                        global_report_internal_gaps_as_tilde;

int                         global_gap_frameshift_mode; // 0: No gaps or frameshift
                                                        // 1: Gaps allowed (default)
                                                        // 2: Only gaps allowed

int                         global_report_gap_mode;     // 1: all with '-'
                                                        // 2: leading, trailing with '-', internal with '~'
                                                        // 3: remove all gaps in output sequence
bool                        global_treat_references_as_individual;
bool                        global_keep_concatenated_input_file;

unsigned                    global_ends_width;
unsigned                    global_weight_fraction_in_ends;

void good_bye_and_exit(int error)
{
  if (error == 0)
    cerr << "## Finished successfully. ##" << endl;
  cerr << endl << PROGNAME << " says goodbye." << endl << endl;
  exit(error);
}

bool directory_exists(const char* dir_path)
{
  // Most backword compatible way to determine whether a directory exists.
  struct stat info;

  int res = stat( dir_path, &info );

  if (res != 0)
    return false;
  return (info.st_mode & S_IFDIR);
}


void init_param()
{
  global_tmp_directory = ".";
  global_num_bp_beyond_exonerate_alignment_if_at_start_or_end = 0;
  global_consensus_threshold = 0.5;
  global_verbosity = 1;
  global_notinclude_gap_alignments = 0;
  global_include_frameshift_alignments = 0;
  global_include_double_hits  = 0;
  global_relative_score_threshold = 1;
  global_genetic_code_number = 2;
  global_frameshift_penalty  = -9;
  global_run_mode = 0;
  global_minimum_seq_coverage_uppercase = 1;
  global_minimum_seq_coverage_total = 1;
  global_report_gap_mode = 1;
  global_exonerate_binary = "exonerate";

  global_keep_concatenated_input_file = false;
  
  global_ends_width = 9;
  global_weight_fraction_in_ends = 20;
}


void read_and_init_parameters(int argc, char** argv)
{
  init_param();

  bool option_set_consensus_sequence_output_filename;
  bool option_set_consensus_threshold;

  try
  {
    CmdLine cmd("Description: This program aligns multiple nucleotide sequences "
		"(typically sequencing reads) against a protein reference "
		"sequence to obtain a multiple sequence alignment. The "
		"recommended use case is to extract protein coding mitochondrial "
		"genes from sequencing libraries. "
		"The individual alignments are computed by calling the Exonerate program. "
		"Exonerate is a very efficient alignment program which allows to align protein and nucleotide sequences."
		"Nucleotide sequences which cannot be aligned to the protein reference will not be included in the output. "
                "Exonerate should be able to align 100k short reads in a few "
		"minutes using a single CPU core, so this approach can be used for projects of any size.\n"
		"Uses cases:\n"
		"- extract protein coding mitochondrial genes (PCMGs) form Illumina sequencing libraries, "
		"e.g. transcriptomic or hybrid enrichment libraries.\n"
		"- extract PCMGs from PacBio libraries, "
                "  attempts to use MinIon long reads have so far not been successful,\n"
		"- extract PCMGs from complete mitochondrial genomes,\n"
		"- extract PCMGs from assembled libraries. This typically has a lower success rate"
		"  compared to extracting genes from not-assembled reads.\n"
		,
		' ', VERSION);
    /*
    ValueArg<string> path_exonerate_binary_Arg("E", "exonerateBinary",
	"Path and filename or only the filename of the Exonerate program. Alternatively, to letting "
	PROGNAME " call Exonerate, the user can provide an existing vulgar output file. "
        "If " PROGNAME " has to call Exonerate, which is the case if specified vulgar file does not exist, "
	"a Protein reference sequence has to be specified. If no executable is specified with this command, "
       PROGNAME " tries to run the \"Exonerate\" command.", false,
       global_exonerate_binary.c_str(), "string");
    cmd.add(path_exonerate_binary_Arg);
    */

    ValueArg<unsigned> verbosity_Arg("", "verbosity",
                                     "Specifies how much run time information is printed to the console. Values: 0: minimal output, 1: important notices, 2: more notices, 3: and basic progress, 4: and detailed progress, 50-100: and debug output, 1000: all output.",
                                     false, global_verbosity, "int");
    cmd.add(verbosity_Arg);


    SwitchArg keep_concat_input_files_Arg("", "keep-concat-input-file",
					  "If multiple input files are "
					  "specified, MGE first creates a "
					  "concatenated file. By default this "
					  "file is removed. Use this option if "
					  "you want to keep this file.", global_keep_concatenated_input_file);
    cmd.add(keep_concat_input_files_Arg);

    SwitchArg treat_references_as_individua_Arg("", "treat-references-as-individual",
					  "Input sequences which can be aligned with "
					  "different reference sequences are by default assigned "
					  "only to the references for which the alignment score is "
					  "equal to the best score achieved by this input sequence. "
					  "This score competition is switched off if this option "
					  "is specified. This treats multiple references as if they "
					  "are specified in independent program runs. ",
				          global_treat_references_as_individual);
    cmd.add(treat_references_as_individua_Arg);
    
    ValueArg<string>   tmp_directory_Arg("", "temporaryDirectory",
					 "MGE has to create potentially large "
					 "temporary files, e.g. if multiple "
					 "input files are specified, or if fastq file "
					 "are specified. With this option these files "
					 "will not be created in the directory the program "
					 "was launched, but in the specified tmp directory. "
                     "Temporary files are created without file extensions.",
					 false, global_tmp_directory, "string");
    cmd.add(tmp_directory_Arg);

    ValueArg<unsigned> min_seq_coverage_upper_case_Arg("", "minSeqCoverageInAlignment_uppercase",
						       "Specifies the absolute value of the minimum alignment coverage for computing the consensus sequence. As coverage, only upper case nucleotides are taken into account, i.e. no nucleotides are counted that have been added beyond the Exonerate alignment region. Bases beyond the Exonerate alignment are added with the -n or --numberOfBpBeyond option. If no bases are added beyond the Exonerate alignment (default), the effect of this option is identical to the minSeqCoverageInAlignment_total option. Default: 1. Increasing this value increases the number of unknown nucleotides in the consensus sequence.",
						       false, global_minimum_seq_coverage_uppercase, "int");
    cmd.add(min_seq_coverage_upper_case_Arg);

    ValueArg<unsigned> min_seq_coverage_total_Arg("", "minSeqCoverageInAlignment_total",
						       "Specifies the absolute value of the minimum alignment coverage for computing the consensus sequence. For the coverage, all nucleotides count, including lower case nucleotides that have been added beyond the Exonerate alignment region. Default: 1. Increasing this value increases the number of unknown nucleotides in the consensus sequence.",
						       false, global_minimum_seq_coverage_total, "int");
    cmd.add(min_seq_coverage_total_Arg);
    
    
    ValueArg<float> relative_score_threshold_Arg("r", "relative_score_threshold",
	"Specified the relative alignment score threshold for Exonerate hits to be considered. "
	"The relative score is the score reported by Exonerate divided by the alignment length. Default 1. "
	"Reasonable thresholds are between 0.7 and 2.0.",
	false, global_relative_score_threshold, "float");
    cmd.add(relative_score_threshold_Arg);

    ValueArg<unsigned> min_exonerate_score_threshold_Arg("s", "minExonerateScoreThreshold",
							 "The score threshold passed to Exonerate to decide whether to include "
							 "or not include the hit in the output. Note that the minimum score implies "
                             "a minimum sequence length that will be reported by exonerate. "
                             "Default: Use exonerate default. For the version 2.4 this is 100.",
						       false, UINT_MAX, "int");
    cmd.add(min_exonerate_score_threshold_Arg);
    
    ValueArg<int> genetic_code_number_Arg("C", "genetic_code",
	 "The number of the genetic code to use in Exonerate, if this step is required. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details. Default: 2, i.e. "
	"vertebrate mitochondrial code. Exonerate only accepts one genetic code, so for multiple references only one genetic code can be used.",
	false, global_genetic_code_number, "int");
    cmd.add(genetic_code_number_Arg);

    ValueArg<int> frameshift_penalty_Arg("f", "frameshift_penalty",
	 "The frameshift penalty passed to Exonerate. Default: -9. Higher values lead to lower scores and by this can have the following effects: (i) hit regions are trimmed since trimming can lead to a better final alignment score, (ii) they can also lead to excluding a read as a whole if the final score is too low and trimming does lead to a higher score. The default of the Exonerate program is -28. A value of -9 (or other values lower than -28) lead to more reads in which the best alignment has a frameshift. In order to remove reads that do not align well, one can use a smaller frameshift penalty and then exclude hits with a frameshift, see -F option).",
	false, global_frameshift_penalty, "int");
    cmd.add(frameshift_penalty_Arg);

    ValueArg<int> report_gaps_mode_Arg("", "report_gaps_mode",
				       "Gaps can be reported in different ways. With this option the reporting mode can be specified. 1: report leading and trailing gaps with \'-\' character. Report internal gaps (introduced with options -G or -g) with \'~\' character. 2: report leading and trailing gaps with \'-\' character. Report internal gaps (introduced with options -G or -g) with \'-\' characters. 3: Remove all gap characters in output. In this case sequences are extracted but are reported with respect to the reference. Default: 1.",
	false, global_report_gap_mode, "int");
    cmd.add(report_gaps_mode_Arg);

    
    //    SwitchArg report_internal_gaps_as_tilde_Arg("", "report_internal_gaps_as_tilde",
    //			 "Report gaps included with -G or -g option as tilde characters. Default: false.",
    //			 false);
    //    cmd.add(report_internal_gaps_as_tilde_Arg);

    
    
    /*
    SwitchArg include_F_Arg("F", "includeFrameshift",
			 "Include reads which aligned with a frameshift. Default: false.",
			 false);
    cmd.add(include_F_Arg);
    */

    SwitchArg include_gaps_in_Reference_Arg("", "reportSequencesWithGapInAAReference_removingNucleotidesInAlignment",
                                            "Include sequences in alignment and consensus sequence for which a gap is suggested in the reference. "
                                            "Note that this removes the nucleotides in sequences that have no "
                                            "corresponding AAs in the reference. Default: false",
                                            false);
    cmd.add(include_gaps_in_Reference_Arg);

    SwitchArg include_gaps_in_NucSeqe_Arg("", "reportSequencesWithGapInNucleotideSequence",
                                          "Include sequences in alignment and conensus sequence that have a gap in the nucleotide sequence. "
                                          "This is recommended. Default: true.",
                                          true);
    cmd.add(include_gaps_in_NucSeqe_Arg);

    SwitchArg only_gaps_in_Reference_Arg("", "reportOnlySequencesWithGapInAAReference_removingNucleotidesInAlignment",
                                         "Include only sequences in alignment and consensus sequence for which a gap is suggested in the reference. "
                                         "Note that this removes the nucleotides in sequences that have no "
                                         "corresponding AAs in the reference. Default: false",
                                         false);
     cmd.add(only_gaps_in_Reference_Arg);

     SwitchArg only_gaps_in_NucSeqe_Arg("", "reportOnlySequencesWithGapInNucleotideSequence",
                                           "Include only sequences in alignment and consensus that have a gap in the nucleotide sequence. "
                                           "This is recommended. Default: true.",
                                           true);
     cmd.add(only_gaps_in_NucSeqe_Arg);


    SwitchArg include_g_only_Arg("g", "onlyGap",
			 "Include only reads which aligned with a gap.",
			 false);
    cmd.add(include_g_only_Arg);
    
    SwitchArg notinclude_Gaps_Arg("", "noGaps",
				  "Do not include reads for which the alignment with the reference contains gaps.",
				  false);
    cmd.add(notinclude_Gaps_Arg);

    //    SwitchArg include_G_Arg("G", "includeGap",
    //			 "Include reads which aligned with a gap.",
    //			 false);
    //    cmd.add(include_G_Arg);


    
    SwitchArg include_D_Arg("D", "includeDoubleHits",
			    "Include input sequences with two alignment results against the same reference.", false);
    cmd.add(include_D_Arg);

    ValueArg<float> consensus_threshold_Arg("t", "consensus_threshold",
	"This option modifies the consensus threshold. Default: 0.5 which corresponds to 50%.",
	false, global_consensus_threshold, "float");
    cmd.add(consensus_threshold_Arg);

    ValueArg<string> consensus_output_file_Arg("c", "consensus_file",
	"Specifies the base name of the consensus sequence output file(s). A "
	"consensus sequence with the name "
	"baseName + reference-sequence-name + .fas is written for each "
	"reference sequence.",
	true, global_consensus_sequence_output_filename, "string");
    cmd.add(consensus_output_file_Arg);

    ValueArg<unsigned> bp_beyond_Arg("n", "numberOfBpBeyond",
	"Specifies the number of base pairs that are shown beyond the Exonerate alignment. A value of 0 means that the sequence is clipped at the point the Exonerate alignment ends. Values >0 can lead to the inclusion of sequence segments that do not align well with the amino acid sequence and have to be treated with caution. They might belong to chimera, Numts, or other problematic sequences. Larger values might be included e.g. if problematic sequences with a well matching seed alignment are of interest. CAUTION: Bases included with this option might not be aligned well or could even belong to stop codons! They should be considered as of lower quality compared to other bases. Bases that are added with this option are added as lower case characters to the output alignment file. A sequence coverage of bases not belonging to these extra bases can be requested with the --minSeqCoverageInAlignment_uppercase option. Default: 0.",
	false, global_num_bp_beyond_exonerate_alignment_if_at_start_or_end, "int");
    cmd.add(bp_beyond_Arg);

    ValueArg<string> exonerate_path_Arg("e", "exonerate_program",
     				        "Specifies the name of the Exonerate program in system path OR the path "
					"to the Exonerate program including the program name. Default: Exonerate",
				        false, global_exonerate_binary, "string");
    cmd.add(exonerate_path_Arg);

    ValueArg<string> vulgar_directory_Arg("V", "vulgar_directory",
				     "Specifies the name of the directory in which vulgar files, generated as exonerate output, are stored. "
                     "For each input file a vular file with the same base name as the input file, but with the extension "
                     ".vulgar is read as input or generated if it does not exist. "
				     "If the vulgar files do not exist MitoGeneExtractor will run Exonerate in "
				     "order to create the files. The created file will then be used to proceed. "
				     "If no directory is specified with this option, the current working directory is used to store the vulgar files. ",
				     false, global_vulgar_directory, "string");
    cmd.add(vulgar_directory_Arg);

    ValueArg<string> alignment_output_file_name_Arg("o", "",
       "Specifies the base name of alignment output file(s). Aligned input sequences "
       "are written to a file with the name: BaseName + sequenceNameOfRefernce + .fas "
       "for each reference sequence.",
       true, global_alignment_output_file, "string");
    cmd.add(alignment_output_file_name_Arg);

    ValueArg<string> prot_fasta_input_file_name_Arg("p", "prot_reference_file",
       "Specifies the fasta file containing the amino acid reference sequences. This file can contain one "
       "or multiple reference sequences. All input nucleotide sequences are aligned against all references. "
       "Hits with a score higher than the minimum are considered. If a sequence matches multiple reference genes/variants, "
       "the sequence will be assigned to the reference for which the alignment score is higher or to both if the scores are equal. ",
       false, global_input_prot_reference_sequence, "string");
    cmd.add(prot_fasta_input_file_name_Arg);

    MultiArg<string> dna_fastq_input_file_names_Arg("q", "dna_fastq_file",
       "Specifies the input query nucleotide sequence files in the fastq format. This option can be "
       "specified multiple times if multiple input files shall be analysed "
       "in one run. All input files will be converted to a fasta file without "
       "taking into account the quality scores. "
       "Sequence files should be quality filtered before being used as input for this program. This "
       "option can be combined with multiple input files in the fasta format (see -d option).",
	false, "string");
    cmd.add(dna_fastq_input_file_names_Arg);

    MultiArg<string> dna_fasta_input_file_names_Arg("d", "dna_fasta_file",
	"Specifies the input query nucleotide sequence files in the fasta format. Sequences are expected not to include gap characters. "
	"This option can be specified multiple times if multiple input files shall be analysed in one run. "
        "If sequence files contain reads, they should have been quality filtered before being used as input for this program. "
	"This option can be combined with multiple input files in the fastq format (see -q option).",
	false, "string");
    cmd.add(dna_fasta_input_file_names_Arg);

    cmd.parse( argc, argv );

    // Assigning parameters to variables:
    global_input_dna_fasta_filenames                   = dna_fasta_input_file_names_Arg.getValue();
    global_input_dna_fastq_filenames                   = dna_fastq_input_file_names_Arg.getValue();

    global_input_prot_reference_sequence               = prot_fasta_input_file_name_Arg.getValue();
    global_tmp_directory                               = tmp_directory_Arg.getValue();

    global_verbosity                                   = verbosity_Arg.getValue();
    global_num_bp_beyond_exonerate_alignment_if_at_start_or_end = bp_beyond_Arg.getValue();
    global_exonerate_binary                            = exonerate_path_Arg.getValue();
    global_vulgar_directory                          = vulgar_directory_Arg.getValue();
    global_alignment_output_file                       = alignment_output_file_name_Arg.getValue();
    global_consensus_sequence_output_filename          = consensus_output_file_Arg.getValue();
    global_consensus_threshold                         = consensus_threshold_Arg.getValue();

    global_include_double_hits                         = include_D_Arg.getValue();
    //    global_include_frameshift_alignments               = include_F_Arg.getValue();
    global_notinclude_gap_alignments                   = notinclude_Gaps_Arg.getValue();
    global_include_only_gap_alignments                 = include_g_only_Arg.getValue();
    global_relative_score_threshold                    = relative_score_threshold_Arg.getValue();

    option_set_consensus_sequence_output_filename      = consensus_output_file_Arg.isSet();
    option_set_consensus_threshold                     = consensus_threshold_Arg.isSet();

    global_genetic_code_number                         = genetic_code_number_Arg.getValue();
    global_frameshift_penalty                          = frameshift_penalty_Arg.getValue();
    global_minimum_seq_coverage_total                  = min_seq_coverage_total_Arg.getValue();
    global_minimum_seq_coverage_uppercase              = min_seq_coverage_upper_case_Arg.getValue();

    global_exonerate_score_threshold                   = min_exonerate_score_threshold_Arg.getValue();
    //    global_report_internal_gaps_as_tilde               = report_internal_gaps_as_tilde_Arg.getValue();
    global_report_gap_mode                             = report_gaps_mode_Arg.getValue();
    global_treat_references_as_individual              = treat_references_as_individua_Arg.getValue();
    global_keep_concatenated_input_file                = keep_concat_input_files_Arg.getValue();
  }

  catch (ArgException &e)
  {
    cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
    good_bye_and_exit(-1);
  }

  //  if (global_)

  if (global_input_dna_fasta_filenames.size() + global_input_dna_fastq_filenames.size() == 0)
  {
    cerr << "ERROR: At least one input filename needs to be specified with the "
            "-d and/or -q option.\nExiting\n";
    good_bye_and_exit(-2);
  }

  if (global_tmp_directory.find('~') != string::npos)
  {
    cerr << "The character \'~\' (tilde) is not allowed in the string "
      "specifying the path to a temporary directory.\nExiting.\n" << endl;
    exit(-13);
  }

  if (global_include_only_gap_alignments && global_notinclude_gap_alignments)
  {
    cerr << "Error: Conflicting command line parameters. It is not possible to "
            "specify to include only sequences that align with a gap and at the "
            "same time no reads that align with a gap in the output.\n";
    good_bye_and_exit(-2);
  }

  global_gap_frameshift_mode = 1;
  if (global_notinclude_gap_alignments)
    global_gap_frameshift_mode = 0;
  else if (global_include_only_gap_alignments)
    global_gap_frameshift_mode = 2;

  if (global_vulgar_directory.empty())
  {
    global_vulgar_directory = "./";
  }
  else if (global_vulgar_directory.back() != '/')
    global_vulgar_directory.push_back('/');

  if (global_consensus_threshold > 1 || global_consensus_threshold < 0)
  {
    cerr << "ERROR: You specified a consensus threshold of "
	 << global_consensus_threshold
	 << ". Allowed values are in the range 0..1, "
            "where 1 corresponds to a strict consensus where 100% identity is "
            "required and 0 implies that the dominant nucleotide "
            "is chosen independent of its proportion."
	 << endl;
    cerr << "Please choose a threshold in the range 0..1 and restart this program. "
            "Exiting."
	 << endl;
    good_bye_and_exit(-4);
  }

  if (option_set_consensus_threshold && !option_set_consensus_sequence_output_filename ) 
  {
    cerr << "ERROR: You specified a consensus threshold but you do not request "
            "to write a consensus sequence file. In order to request to write a "
            "consensus sequence file, you have to specify the -c or "
            "--consensus_file option followed by a file name";
    good_bye_and_exit(-5);
  }

  global_tmp_directory.erase(global_tmp_directory.find_last_not_of(" \n\r\t/")+1);
  if (!directory_exists(global_tmp_directory.c_str()))
  {
    cerr << "ERROR: The specified temporary directory:\n"
         << global_tmp_directory
         << "\ndoes not exist. Exiting\n";
    good_bye_and_exit(-6);
  }

  if (!directory_exists(global_vulgar_directory.c_str()))
  {
    cerr << "ERROR: The specified directory for vulgar files:\n"
    << global_vulgar_directory
    << "\ndoes not exist. Exiting\n";
    good_bye_and_exit(-6);
  }


}



void print_parameters(std::ostream &os, const char *s)
{
  os << s <<   "Parameter settings:"
     << std::endl
     << s <<   "===================" 
     << std::endl;

  if (global_input_dna_fasta_filenames.size() > 0)
  {
    os << s <<   "DNA fasta input file names:                             ";
    size_t i;
    for (i=0; i<global_input_dna_fasta_filenames.size()-1; ++i)
      os << global_input_dna_fasta_filenames[i] << ",";
    os << global_input_dna_fasta_filenames[i] << '\n';
  }

  if (global_input_dna_fastq_filenames.size() > 0)
  {
    os << s <<   "DNA fastq input file names:                             ";
    size_t i;
    for (i=0; i<global_input_dna_fastq_filenames.size()-1; ++i)
      os << global_input_dna_fastq_filenames[i] << ",";
    os << global_input_dna_fastq_filenames[i] << '\n';
  }

  if (global_input_dna_fastq_filenames.size() + global_input_dna_fasta_filenames.size() > 1)
  {
    os << s <<   "Keep concatenated input file:                           " << (global_keep_concatenated_input_file ? "yes":"no")
       << std::endl;
  }

  os << s << "Protein reference input file name:                      " << global_input_prot_reference_sequence
     << std::endl;

  os << s << "Directory for temporary files:                          " << global_tmp_directory
     << std::endl;

  os << s << "Alignment output file:                                  " << global_alignment_output_file
     << std::endl;

  os << s << "Vulgar directory:                                       " << global_vulgar_directory
     << std::endl;

  if (!global_exonerate_binary.empty()) 
    os << s << "Exonerate binary:                                       " << global_exonerate_binary
       << std::endl;

  os << s <<   "Genetic code (NCBI genetic code number):                " << (short) global_genetic_code_number
     << std::endl;
  
  os << s <<   "Print this number of bp beyond Exonerate alignment:     " << global_num_bp_beyond_exonerate_alignment_if_at_start_or_end
     << std::endl;

  if (!global_consensus_sequence_output_filename.empty())
  {
    os << s << "Write consensus sequence to file :                      " << "yes"
       << std::endl;
    os << s << "Filename for consensus sequence output:                 " << global_consensus_sequence_output_filename
       << std::endl;
    os << s << "Consensus sequence threshold value:                     " << global_consensus_threshold
       << std::endl;
  }
  else
  {
    os << s << "Write consensus sequences to file:                      " << "no"
       << std::endl;
  }

  os << s <<   "Frameshift penalty:                                     " << global_frameshift_penalty
     << std::endl;

  os << s <<   "Relative score threshold:                               " << global_relative_score_threshold
     << std::endl;

  os << s <<   "Minimum coverage in Exonerate alignment:                " << global_minimum_seq_coverage_total
     << std::endl;

  os << s << "Minimum coverage in Exonerate alignment (upper case):   " << global_minimum_seq_coverage_uppercase
     << std::endl;
  
  if (global_exonerate_score_threshold != UINT_MAX)
    os << s << "Exonerate score threshold:                              " << global_exonerate_score_threshold
       << std::endl;    

  if (global_gap_frameshift_mode == 0)
  {
    os << s <<   "Gappy reads used:                                       no\n";
    os << s <<   "Frameshift reads used:                                  no\n";
  }
  else if (global_gap_frameshift_mode == 1)
  {
    os << s <<   "Gappy reads used:                                       yes\n";
    //    os << s <<   "Report internal gaps as tilde:                          " << (report_internal_gaps_as_tilde_Arg ? "yes\n" : "no\n");
    os << s <<   "Frameshift reads used:                                  no\n";
  }
  else if (global_gap_frameshift_mode == 2)
  {
    os << s <<   "Only gappy reads used:                                  yes\n";
    //    os << s <<   "Report internal gaps as tilde:                          " << (report_internal_gaps_as_tilde_Arg ? "yes\n" : "no\n");
    os << s <<   "Frameshift reads used:                                  no\n";
  }

  os << s <<   "Treat all references as independent:                    " << (global_treat_references_as_individual ? "yes\n" : "no\n");
  
  if (global_report_gap_mode == 1)
    os << s <<   "Report gaps mode:                                       Report all (leading, trailing, internal) gaps with \'-\' character.\n";
  else if (global_report_gap_mode == 2)
    os << s <<   "Report gaps mode:                                       Report leading and trailing gaps with \'-\' character and internal gaps with \'~\' character.\n";
  else if (global_report_gap_mode == 3)
    os << s <<   "Report gaps mode:                                       Remove all gaps in output.\n";

  //  os << s << "\n";

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


