#include <iostream>
#include <fstream>
#include "faststring2.h"
#include <cstdlib>
#include "CSequences2.h"
#include "CSequence_Mol2_1.h"
#include "Ctriple.h"
#include <vector>
#include "CFile/CFile2_1.h"
#include <map>
#include "CDnaString2.h"
#include "global-types-and-parameters_MitoGeneExtractor.h"
#include <iomanip>
#include <cstdio>
#include "statistic_functions.h"
/// #include <utility>

using namespace std;

typedef pair<unsigned, unsigned> Key;
typedef map< Key , unsigned> Mymap;

// Verbosity rules:
//      verbosity 0 -> only error messages that lead to exit() (goes to cerr). Welcome is printed to cout. Parameters are printed to cout.
//      verbosity 1 -> warnings, ALERT (goes to cerr)
//      verbosity 2 -> more warnings, e.g. skipping features
//      verbosity 3 -> Basic progress (goes to cout)
//      verbosity 4 -> more progress (goes to cout)
//      verbosity >=20: DEGUB
// Default: 1


// What the exonerate manual says about coordinates:
// The coordinates are the coordinates in between the bases:
//      A C G T
//     0 1 2 3 4
// This is equivalent to a 0 based system in which the end is indicated by the
// next index after the range.

//bool 	  add_partial_codons_at_ends      = true;
//unsigned  length_constaint                = 2000;

unsigned num_WARNINGS = 0;

inline void add_or_count(Mymap &m, Key k)
{
  if (m.find(k) == m.end())
    m.insert(make_pair(k,1));
  else
    ++m[k];
}

void print_Mymap(ostream &os, Mymap &m)
{
  Mymap::iterator it, it_end;

  it     = m.begin();
  it_end = m.end();

  while (it != it_end)
  {
    os << "(" << it->first.first << "," << it->first.second << "):" << it->second << endl;
    ++it;
  }
}

inline void add_or_count(std::map<faststring, unsigned> &m, faststring &x)
{
  std::map<faststring,unsigned>::iterator it;
  
  it = m.find(x);
  if (it == m.end() )
  {
    m[x] = 1;
  }
  else
  {
    ++it->second;
  }
}

void print_DNA_profile(ostream &os, unsigned **profile, unsigned N_sites)
{
  os << "pos\tA\tC\tG\tT\t-\tambig\ttilde\n";
  for (unsigned i=0; i<N_sites; ++i)
  {
    os << i+1 << '\t' << profile[i][0] << '\t' << profile[i][1] << '\t' << profile[i][2]
       << '\t' << profile[i][3] << '\t' << profile[i][4] << '\t' << profile[i][5] << '\t' << profile[i][6] << '\n';
  }
}


class vulgar
{
public:
  
  faststring query;          // In this program, this is the AA-COI consensus seq. exonerate aligned against.
  unsigned   query_start;    // 0 based numbers.
  unsigned   query_end;      // The first position after the specified range.
  char       query_strand;
  faststring target;         // In this program, they are the reads.
  unsigned   target_start;   // 0 based numbers.
  unsigned   target_end;     // The first position after the specified range.
  char       target_strand;
  unsigned   score;
  vector<Ctriple<char, unsigned, unsigned> > attributes;
  
  // Frequent labels:
  bool       bool_has_M;
  bool       bool_has_G;
  bool       bool_has_F;
  // Rare labels:
  bool       bool_has_5;
  bool       bool_has_3;
  bool       bool_has_C;
  bool       bool_has_N;
  bool       bool_has_I;
  bool       bool_has_S;
  bool       bool_has_non_M;
  
  void parse(faststring &str)
  {
    bool_has_M = false;
    bool_has_G = false;
    bool_has_F = false;
    bool_has_5 = false;
    bool_has_3 = false;
    bool_has_C = false;
    bool_has_N = false;
    bool_has_I = false;
    bool_has_S = false;
    bool_has_non_M = false;
    
    attributes.clear();
    
    vector<faststring> l(19);
    split(l, str);
    
    if (l.size() < 13) // 1+9+3
    {
      
      cerr << "ERROR: Bad vulgar string encountered with a wrong number of elements: " << str << endl;
      exit(-4);
    }

    if (l[0] != "vulgar:")
    {
      cerr << "ERROR: Bad vulgar string: The vulgar string is expected to start with \"vulgar:\" but found: " << str << endl;
      exit(-4);
    }

    // It is a convention in exonerate that the protein sequence is the query and the DNA sequence the target:
    query         = l[1];
    query_start   = l[2].ToUnsigned();
    query_end     = l[3].ToUnsigned();
    query_strand  = l[4][0];

    target        = l[5];
    target_start  = l[6].ToUnsigned();
    target_end    = l[7].ToUnsigned();
    target_strand = l[8][0];

    score         = l[9].ToUnsigned();

    int m = l.size()-10;
    if (m%3 != 0)
    {
      cerr << "ERROR: Bad vulgar string encountered. The number of fields must by 10+3*n where n is an integer. But the number is "
	   << m << "." << endl;
      cerr << "       Vulgar string: " << str << endl;
      exit(-5);
    }

    unsigned k = 10;
    Ctriple<char, unsigned, unsigned> t('\0', 0, 0);
    while (k < l.size())
    {
      char c        = l[k][0]; ++k;
      unsigned num1 = l[k].ToUnsigned(); ++k;
      unsigned num2 = l[k].ToUnsigned(); ++k;
      
      if (c == 'M')
        bool_has_M = true;
      else if (c == 'G')
      {
        bool_has_G     = true;
        bool_has_non_M = true;
      }
      else if (c == 'F')
      {
        bool_has_F = true;
        bool_has_non_M = true;
      }
      else if (c == '5')
      {
        bool_has_5 = true;
        bool_has_non_M = true;
      }
      else if (c == '3')
      {
        bool_has_3 = true;
        bool_has_non_M = true;
      }
      else if (c == 'C')
      {
        bool_has_C = true;
        bool_has_non_M = true;
      }
      else if (c == 'N')
      {
        bool_has_N = true;
        bool_has_non_M = true;
      }
      else if (c == 'I')
      {
        bool_has_I = true;
        bool_has_non_M = true;
      }
      else if (c == 'S')
      {
        bool_has_S = true;
        bool_has_non_M = true;
      }
      
      t = Ctriple<char, unsigned, unsigned>(c, num1, num2);
      attributes.push_back(t);
    }
  }
  
  /* Does not seem to make sense to combine F and G. There is a simple and consistent interpretation of the two independent of the other.
   void combine_FG_attributes()
   {
   unsigned i, N=attributes.size();
   
   for (i=1; i<N; ++i)
   {
   if (attributes[i-1].first() == 'G' && attributes[i].first() == 'F')
   {
   if (attributes[i-1].second()>0 && attributes[i-1].third() == 0 && (attributes[i].second()==0 && attributes[i-1].third() > 0) )
   {
   Ctriple<char, unsigned, unsigned> t('\0', 0, 0);
   t = Ctriple<char, unsigned, unsigned>('f', attributes[i-1].second(), attributes[i].third());
   attributes[i-1] = t;
   attributes.erase(attributes.begin()+i);
   // Do we have to adapt i after removing an element. In this case no, since we will not combine anything with the newly inserted feature 'f'.
   }
   else
   {
   cout << "WARNING: Unexpected combination of lengths in successive G and F attribute in vulgar line:\n";
   print(cout);
   }
   }
   else if (attributes[i-1].first() == 'F' && attributes[i].first() == 'G')
   {
   if (attributes[i-1].second()==0 && attributes[i-1].third() >0 && (attributes[i].second()>0 && attributes[i].third() == 0) )
   {
   Ctriple<char, unsigned, unsigned> t('\0', 0, 0);
   t = Ctriple<char, unsigned, unsigned>('f', attributes[i].second(), attributes[i-1].third());
   attributes[i-1] = t;
   attributes.erase(attributes.begin()+i);
   // Do we have to adapt i after removing an element. In this case no, since we will not combine anything with the newly inserted feature 'f'.
   }
   else
   {
   cout << "WARNING: Unexpected combination of lengths in successive F and G attributes in vulgar line:\n";
   print(cout);
   }
   }
   }
   }
   */
  
public:
  void print(ostream &os)
  {
    os << "vulgar: "
    << query << " "
    << query_start   << " "
    << query_end     << " "
    << query_strand  << " "
    << target        << " "
    << target_start   << " "
    << target_end     << " "
    << target_strand  << " "
    << score << " ";
    for (unsigned i=0; i<attributes.size(); ++i)
    {
      os << attributes[i].first()  << " ";
      os << attributes[i].second() << " ";
      os << attributes[i].third() << " ";
    }
  }
  
  vulgar(faststring &str)
  {
    parse(str);
    //    combine_FG_attributes();
  }
  
  const faststring & get_target()
  {
    return target;
  }
  
  bool has_F()
  {
    return bool_has_F;
  }
  
  bool has_G()
  {
    return bool_has_G;
  }
  
  bool has_M()
  {
    return bool_has_M;
  }
  
  bool has_non_M()
  {
    return bool_has_non_M;
  }
  
  unsigned num_attributes()
  {
    return attributes.size();
  }
  
  unsigned get_query_start()
  {
    return query_start;
  }

  unsigned get_target_start()
  {
    return target_start;
  }
  
  unsigned get_target_end()
  {
    return target_end;
  }
  
  bool is_revcomp()
  {
    return (target_strand == '-');
  }

  double relative_score()
  {
    double dist;
    if (target_end > target_start)
      dist = target_end - target_start;
    else
      dist = target_start - target_end;
    if (dist == 0)
      return 0;
    return score/dist;
  }
};

bool fileExists(const char *filename)
{
  ifstream is(filename);
  if (is.fail())
    return false;
  else
    return true;
}


int run_exonerate(const faststring &exonerate_binary, const faststring &protein_reference, const faststring &fasta_reads, const faststring &vulgar_base_file, char genetic_code_number, int frameshift_penalty, int max_trials, char mode=0)
{
  faststring code_string;

  if (genetic_code_number == 2) // The Vertebrate Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
  else if (genetic_code_number == 3)  // The Yeast Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 4)  // The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 5)  // The Invertebrate Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 9)  // The Echinoderm and Flatworm Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 13) // The Ascidian Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 14) // The Alternative Flatworm Mitochondrial Code
    code_string = "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 16) // Chlorophycean Mitochondrial Code
    code_string = "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else
  {
    cerr << "ERROR: Unknown genetic code number passed to function run_exonerate. Code nuber: " << genetic_code_number << endl;
    exit(-20);
  }
  
  //  system "exonerate --geneticcode FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG --frameshift -9 --query ${dir2}/${gene}.fasta 
  // -Q protein --target $dir1/${name}.fas -T dna --model protein2dna --showalignment 0 --showvulgar 1 > $dir3/${gene}_${name}_vulgar.txt 2> $dir3/${gene}_${name}_vulgar.log";

  int system_return_value;
  faststring score_threshold;
  if (global_exonerate_score_threshold != UINT_MAX)
    score_threshold = " -s " + faststring(global_exonerate_score_threshold) + " ";
  faststring command =   exonerate_binary + " --geneticcode " + code_string + score_threshold
                                          + " --frameshift " + faststring(frameshift_penalty)
                       + " --query " + protein_reference + " -Q protein --target " + fasta_reads
                       + " -T dna --model protein2dna --showalignment 0 --showvulgar 1 1> "
                       + vulgar_base_file + " 2> " + vulgar_base_file + ".log";

  // Special debug mode to examine command string:
  if (mode == 'd')
  {
    cout << "DEBUG: Exonerate command: " << command << endl;
    return 0;
  }

  int i;
  for (i=1; i <= max_trials; ++i)
  {
    system_return_value = system(command.c_str());
    if (system_return_value == 0)
      return i;
    // Else
    // Clean up the unsuccessful run:
    //    sleep(1); // Sleep for 1 second to allow system to close all files.
    faststring remove_file = vulgar_base_file;

    remove(remove_file.c_str());
    faststring log_file_name = vulgar_base_file + ".log";
    remove(log_file_name.c_str());
  }
  return max_trials +1;
}


int main(int argc, char **argv)
{
  cout << "Welcome to the " PROGNAME " program, version " VERSION << endl;
  
  // In exonerate the protein sequence is by convention the query sequence.
  // This makes things a bit confusing here:
  
  read_and_init_parameters(argc, argv);
  print_parameters(cout, "");
  
  CSequences2 seqs_DNA_reads_target(CSequence_Mol::dna);
  CSequences2 seqs_prot_query(CSequence_Mol::protein);
  CSequences2 seqs_DNA_result(CSequence_Mol::dna);

  unsigned skipped_double_vulgar=0, skipped_G=0, skipped_F=0, skipped_relative_score=0, skipped_no_G=0;
  
  //**************************
  // Read query sequences
  //**************************
  CFile query_file;
  
  if (global_verbosity >= 3)
  {
    cout << "Reading query sequences from fasta file: " << global_input_dna_fasta_file << endl;
  }
  
  query_file.ffopen(global_input_dna_fasta_file.c_str());
  if (query_file.fail())
  {
    cerr << "ERROR: Could not open specified file: " << global_input_dna_fasta_file << endl;
    exit(-1);
  }
  // cerr << "Line:   " <<  file.line() << endl;
  // cerr << "Status: " <<  (int)file.status() << endl;
  
  CSequence_Mol::processing_flag pflag(CSequence_Mol::convert_toupper);
  int error = seqs_DNA_reads_target.read_from_Fasta_File(query_file, pflag, 0, -1, false);
  query_file.ffclose();
  
  if (error < 0)
  {
    cerr << "ERROR: Exiting: An error occurred while reading the input file: " << global_input_dna_fasta_file
    << " Double check the input type or the file you specified."    << endl;
    exit(-3);
  }

  // TODO:
  // Instead of printing this error message and exiting, we could fix this problem at this point, at least
  // if we compute the vulgar file ourselves.

  if (!seqs_DNA_reads_target.are_short_names_unique())
  {
    cerr << "ERROR: The sequence names in the DNA read file are not unique when trimming them at the first space, "
            "which is done by many programs, including exonerate. Please verify your input data. Often it "
            "helps to rename sequences by replacing spaces with e.g. underscores.\n";
    exit(-3);
  }
  
  if (global_verbosity >= 3)
  {
    cout << "Finished reading target read sequences. Found this number of target sequences: " << seqs_DNA_reads_target.GetTaxaNum()  << endl;
  }
  
  // Print memory usage details:
  {
    unsigned long mem1, mem2, mem3;
    float mean2;
    seqs_DNA_reads_target.memory_usage(mem1, mem2, mem3, mean2);
    
    if (global_verbosity > 2000)
    {
      cout << "Estimated memory usage of sequences in memory: " << endl
      << "class CSequences2:           " << mem1 << endl
      << "seqData vector content:      " << mem2 << endl
      << "seqData vector content: (GB) " << (float)mem2/1024/1024/1024 << endl
      << "sn_map content:              " << mem3 << endl
      << "mean per seq object:         " << mean2 << endl;
    }
  }
  
  
  //**************************
  // Read query sequence
  //**************************
  CFile query_prot_file;
  
  if (global_verbosity >= 3)
    cout << "Reading fasta file with protein reference sequence: " << global_input_prot_reference_sequence << endl;
  
  query_prot_file.ffopen(global_input_prot_reference_sequence.c_str());
  if (query_prot_file.fail())
  {
    cerr << "ERROR: Could not open specified file: " << global_input_prot_reference_sequence << endl;
    exit(-1);
  }
  // cerr << "Line:   " <<  file.line() << endl;
  // cerr << "Status: " <<  (int)file.status() << endl;
  
  error = seqs_prot_query.read_from_Fasta_File(query_prot_file, pflag, 0, -1, false);
  query_prot_file.ffclose();
  
  if (error < 0)
  {
    cerr << "ERROR: Exiting: An error occurred while reading the input file. "
            "Double check the input type or the file you specified."
	 << endl;
    exit(-3);
  }
  
  if (seqs_prot_query.GetTaxaNum() != 1)
  {
    cerr << "ERROR: Found multiple sequences in the query protein file. Only one sequence is allowed in this program." << endl;
    exit(-23);
  }
  unsigned query_prot_length = seqs_prot_query.GetPosNum();
  
  if (global_verbosity >= 3)
  {
    cout << "Finished reading fasta file with protein reference sequence." << endl;
  }
  
  
  //**************************
  // Export query if sequence names had to be modified.
  //**************************
  
  //**************************
  // Run exonerate
  //**************************
  bool vulgar_file_read_successfully = false;
  unsigned count_appempts_to_read_vulgar_file = 0;

  list<pair<faststring, vulgar> > list_exonerate_results;
  map<faststring, unsigned>       map_exonerate_count_results;

  while (!vulgar_file_read_successfully)
  {
    list_exonerate_results.clear();
    map_exonerate_count_results.clear();

    ++count_appempts_to_read_vulgar_file;
    faststring exonerate_output_file_name = global_vulgar_file.c_str();

    // Do we need to compute the vulgar file?
    if ( !fileExists(global_vulgar_file.c_str()))
    {
      if(global_exonerate_binary.empty() )
	global_exonerate_binary = "exonerate";

      if (global_verbosity >= 1 && count_appempts_to_read_vulgar_file == 1)
	cout << "No vulgar file with the specified name has been found. So the binary \"" << global_exonerate_binary << "\" will be used to create the vulgar file." << endl;

      int num_trials = run_exonerate(global_exonerate_binary.c_str(), global_input_prot_reference_sequence.c_str(),
				     global_input_dna_fasta_file.c_str(), global_vulgar_file.c_str(), global_genetic_code_number,
				     global_frameshift_penalty, 10, global_run_mode);

      if (num_trials == 10+1)
      {
	cerr << "ERROR: Running exonerate failed. The generated vulgar file is incomplete and should be removed manually. Exiting.\n";
	exit(-1);
      }
    
      if (global_verbosity >= 1)
	cout << "Exonerate result file was created successfully after N = " << num_trials << " trials." << endl;
    }
    else
    {
      if (global_verbosity >= 1)
	cout << "The specified vulgar file exists, so it does not have to be recomputed with exonerate." << endl;
    }

    //**************************
    // Parse exonerate output:
    //**************************

    CFile exonerate_output_file;
    //  map<faststring, vulgar>  m_exonerate_results;

    exonerate_output_file.ffopen(exonerate_output_file_name.c_str());

    if (exonerate_output_file.fail())
    {
      cerr << "ERROR: Could not open file " << exonerate_output_file_name << endl;
      cerr << "Exiting.\n\n";
      exit(-6);
    }

    faststring line;
    while(!exonerate_output_file.fail())
    {
      //    char c = exonerate_output_file.getchar();
      //    cout << "char: " << c << endl;

      exonerate_output_file.getline(line);
      if (line.size() > 7 && line.substr(0,7) == "vulgar:")
      {
	//       cout << "Adding line: " << line << endl;
	vulgar v(line);
	//       m_exonerate_results.insert(make_pair(v.get_target(), v) );
	list_exonerate_results.push_back(make_pair(v.get_target(),v));
	add_or_count(map_exonerate_count_results, v.get_target());
      }
      else
      {
	if (line.size() >= 13 && line.substr(0,13) == "Command line:")
	  continue;
	if (line.size() >= 9 && line.substr(0,9) == "Hostname:")
	  continue;
	if (line.size() >= 31 && line.substr(0,31) ==  "-- completed exonerate analysis")
	{
	  vulgar_file_read_successfully = true;
	  break;
	}
	line.removeSpacesBack();
	if (line.empty()) // Skip empty lines even though they normally do not occur in vulgar files.
	  continue;
	cerr << "WARNING: Ignoring non standard line in vulgar file: \"" << line << "\""<< endl;
      }
    } // END  while(!exonerate_output_file.fail())
    if (!vulgar_file_read_successfully)
    {
      cerr << "WARNING: Vulgar file is incomplete. It does not end with \"-- completed exonerate analysis\"\n";
      cerr << "This is attempt number " << count_appempts_to_read_vulgar_file
	   << " to read the vulgar file. I will try to recompute the vulgar file.\n";
      exonerate_output_file.ffclose();
      remove(exonerate_output_file_name.c_str());
    }
    else
    {
      exonerate_output_file.ffclose();
    }
  } // END  while (!vulgar_file_read_successfully)
  ////////////////////////

  if (global_vulgar_file == "tmp-vulgar.txt")
    remove("tmp-vulgar.txt");

  // Debug output:
  if (0)
  {
    std::cout << "Finished reading exonerate results. Please hit any key to continue.\n" << std::endl;
    std::cin.get();
  }
  
  if (global_verbosity >= 3)
  {
    unsigned s1 = list_exonerate_results.size();
    cout << "Number of successful exonerate alignments (could include multiple hits) to the amino acid sequence found in input file: " << s1 << endl;
    unsigned s2 = map_exonerate_count_results.size();
    cout << "Number of successful exonerate alignments without multiple hits. Normally the number of aligned reads:                  " << s2 << endl;
  }
  
  // Debug output:
  if (global_verbosity >= 100)
  {
    cout << "DEBUG OUTPUT: Exonerate results.\n";
    
    list<pair<faststring, vulgar> >::iterator it;
    it = list_exonerate_results.begin();
    while (it !=  list_exonerate_results.end())
    {
      faststring key = it->first;
      cout << " => item " << key << " ";
      it->second.print(cout);
      cout << endl;
      ++it;
    }
  }
  
  // Debug output:
  if (global_verbosity >= 100)
  {
    vector<faststring> names;
    seqs_DNA_reads_target.get_sequence_names(names);
    
    for (unsigned i=0; i<names.size(); ++i)
    {
      cout << "Name " << i << " " << names[i] << endl;
    }
  }
  
  ofstream os;
  os.open(global_alignment_output_file.c_str());
  
  // After reading the data, this is the main block.
  {
    list<pair<faststring, vulgar> >::iterator l_it;
    //    map<faststring, vulgar>::iterator m_it;
    
    unsigned seq_count = 0;
    
    //**************************
    // For each sequence, use exonerate output to align the sequences:
    //**************************
    // Loop over all vulgar strings and align the corresponding read:
    // Sequences without vulgar string are ignored. - Diagnostic output could be written for these sequences.
    
    bool consider_this_read;
    unsigned count_not_considered = 0;

    map<pair<unsigned, unsigned>, unsigned> map_of_insertion_sites_suggested_in_reference, map_of_gap_sites_in_queries;
    
    l_it = list_exonerate_results.begin();
    while (l_it !=  list_exonerate_results.end()) // For all hits:
    {
      consider_this_read = true;

      ++seq_count;
      faststring &key = l_it->first;
      vulgar     &vul = l_it->second;
      
      if (global_verbosity >= 4)
      {
        cout << "Working on key: " << key << endl;
      }
      
      //      m_it = m_exonerate_results.find(key);
      
      //       if (m_it == m_exonerate_results.end() )
      //       {
      // 	cerr << "INTERNAL ERROR: No sequence name found for map key: " << key << endl;
      // 	exit(-8);
      //       }
      
      //      vulgar vul = m_it->second;
      
      //       cout << "XX-Working on vulgar : ";
      //       vul.print(cout);
      //       cout << endl;
      
      map<faststring, unsigned>::iterator find_it = map_exonerate_count_results.find(key);
      bool double_vulgar;
      
      if (find_it != map_exonerate_count_results.end() && find_it->second > 1)
      {
        double_vulgar = true;
      }
      else
      {
        double_vulgar = false;
      }
      
      //      if ( vul.has_non_M() || double_vulgar )
      if ( double_vulgar )
      {
        ++num_WARNINGS;
        if (global_verbosity >= 1)
          cerr << "WARNING: Found two or more exonerate hits for a single target read." << endl;
        
        if (!global_include_double_hits)
        {
          consider_this_read = false;
	  ++skipped_double_vulgar;
        }
        
      }
      if (vul.has_F() && !global_include_frameshift_alignments)
      {
        consider_this_read = false;
	++skipped_F;
      }
      
      if ((vul.has_G() && global_gap_frameshift_mode == 0) )
      {
        consider_this_read = false;
	if (!vul.has_F() )
	  ++skipped_G;
      }

      if ((!vul.has_G() && global_gap_frameshift_mode == 2) )
      {
        consider_this_read = false;
	if (!vul.has_F() )
	  ++skipped_no_G;
      }

      if (global_relative_score_threshold && vul.relative_score() < global_relative_score_threshold)
      {
	cerr << "NOTE: Exonerate hit skipped due to low relative alignment score: " << vul.relative_score()
	     << " " << vul.get_target() << endl;
	consider_this_read = false;
	if (!vul.has_F() && !vul.has_G())
	  ++skipped_relative_score;
      }

      if (!consider_this_read)
      {
	++count_not_considered;
      }
      else //      if (consider_this_read)
      {
        CDnaString seqDNA, newseqDNA;
        unsigned start_in_seqDNA, end_in_seqDNA;

        // The query prot sequence start coordinate sets the number of gaps before:
        unsigned gaps_before = vul.query_start*3;

        CSequence_Mol* theSeq = seqs_DNA_reads_target.get_seq_by_name(key);
        if (theSeq == NULL)
        {
          cerr << "ERROR: No sequence for key " << key << " found. This can only happen if the vulgar file does not correspond to the fasta file. Exiting.\n";
          exit(-13);
        }

        unsigned length_this_nuc_seq = theSeq->length();
        // Compute the rev comp once, so we do not have to handle so many different cases.
        if (vul.is_revcomp())
        {
          seqDNA.setToReverseComplementOf( theSeq->getSeq_faststring());
/*
          unsigned   len = vul.get_target_start() - vul.get_target_end();
          CDnaString tmp = theSeq->getPartialSeq(vul.get_target_end(),len);
          seqDNA.setToReverseComplementOf( tmp.c_str() );
*/
          start_in_seqDNA = length_this_nuc_seq - vul.get_target_start();
          end_in_seqDNA   = length_this_nuc_seq - vul.get_target_end();
        }
        else
        {
          // See coordinate system above. The coordinates are 0 based. The end is the index after the range.
          seqDNA = theSeq->getSeq_faststring();
          /*
            unsigned len = vul.get_target_end() - vul.get_target_start();
            //->getPartialSeq(vul.get_target_start(),len);
          */
          start_in_seqDNA = vul.get_target_start();
          end_in_seqDNA   = vul.get_target_end();
        }

        /* // TODO: What do we do with this code. Add options!!!
         if (vul.num_attributes() != 1) // if we have no non-M labels we expect only one attribute in the list of attributes.
         {
         cerr << "ERROR: Unexpected number of attributes for exonerate output: ";
         vul.print(cerr);
         cerr << "\n";
         exit(-19);
         }
         */

        // Problem: exonerate only aligns complete codons to amino acids

        // Test for consistent bases in partial codons:
        
        // Exonerate aligns the reads to the amino acid sequence. Only complete codons are considered as alignment matches.
        // So using the coordinates obtained from exonerate, a partial codon is removed in 2/3 of the cases at the beginning as well as at the end.
        //

        // We might want to add unmatched bp before and after, e.g. for partial codons or in order to see how much they differ from other sequences.
        // In order not to intervene with with multiple attributes, we treat that case
        // that we add unmatched residues before and after as special cases.
        // If would only expect to find a single M attribute, the logic would be much simpler.

        faststring add_seq_before;
        faststring add_seq_after;

	//        global_num_bp_beyond_exonerate_alignment_if_at_start_or_end = 6;

        // Code can still be simplified !!!
	// Handling recomp is done by computing the reverse complement already for the template sequence.
	// This makes the following code simpler! The final step, i.e. not to distinguishing the two cases any more could now be taken.
        if (global_num_bp_beyond_exonerate_alignment_if_at_start_or_end > 0)
        {
          if (gaps_before > 0)
          {
            if (global_verbosity >= 1000)
            {
              cout << endl;
              cout << "Sequence:    " << seq_count << endl;
              cout << "gaps_before: " << gaps_before << endl;
            }

            if (vul.is_revcomp())
            {
              unsigned add_length_beginning = 0;
              unsigned new_start_in_seqDNA;
              
              if (global_verbosity >= 1000)
                cout << "Number add_length_beginning: " << add_length_beginning << endl;

              if (start_in_seqDNA > 0) // // Could be 0, and nothing should be done
              {
		add_length_beginning = start_in_seqDNA;
		if (global_num_bp_beyond_exonerate_alignment_if_at_start_or_end < add_length_beginning)
		  add_length_beginning = global_num_bp_beyond_exonerate_alignment_if_at_start_or_end;
		if (gaps_before < add_length_beginning)
		  add_length_beginning = gaps_before;

		new_start_in_seqDNA = start_in_seqDNA - add_length_beginning;
		if (add_length_beginning > 0)
		  add_seq_before = seqDNA.substr(new_start_in_seqDNA, add_length_beginning);
              }
            }
            else // Not revcomp:
            {
              unsigned add_length_beginning = 0;
              unsigned new_start_in_seqDNA;

              // We could also handle the case:
              // vul.get_target_start() > global_num_bp_beyond_exonerate_alignment_if_at_end
              // not a condition any more.
              if ( start_in_seqDNA > 0) // Could be == 0.
              {
		add_length_beginning = start_in_seqDNA;
		if (global_num_bp_beyond_exonerate_alignment_if_at_start_or_end < add_length_beginning)
		  add_length_beginning = global_num_bp_beyond_exonerate_alignment_if_at_start_or_end;
		if (gaps_before < add_length_beginning)
		  add_length_beginning = gaps_before;

		new_start_in_seqDNA = start_in_seqDNA - add_length_beginning;
		if (add_length_beginning > 0)
		  add_seq_before = seqDNA.substr(new_start_in_seqDNA, add_length_beginning);
              }

	      if (0)
	      {
		cerr << "DEBUG: Information on sequence data printed in lower case before the matching sequence:\n";
		cerr << "Sequence name: " << key << endl;
		cerr << "Starting coordinate of + strand: " << new_start_in_seqDNA  << endl;
		cerr << "Length:                          " << add_length_beginning << endl;
		cerr << "String:                          " << add_seq_before       << endl;
		
		//            faststring seq_before_start = add_partial_condons_before;
                
                // 		cout << endl;
                // 		cout << "key: " << key << endl;
                
                
                // 		seq_before_start += "|";
                // 		seq_before_start += theSeq->getPartialSeq(vul.get_target_start(), 12);
                // 		faststring aa_at_seq_before_start = seqs_target.get_seq_by_index(0)->getPartialSeq(vul.query_start-1, 5);
                // 		cout << "gaps_before: " << gaps_before << " vul.get_target_start(): " << vul.get_target_start() << endl;
                // 		cout << "Sequence before start: " << seq_before_start << endl;
                // 		cout << "AA at this posiiton:   " << aa_at_seq_before_start << endl;
	      }
	    } // END Not revcomp
	    
	    if (global_verbosity >= 1000)
	      cout << "Add partial codon front \"" << add_seq_before << "\" might have length 0." << endl;
            
	    if (global_verbosity >= 50)
	      cout << "Add front length: " << add_seq_before.length() << " " << key << endl;
            
	  }  // END if (gaps_before > 0)
	} // END if (global_num_bp_beyond_exonerate_alignment_if_at_end > 0)

        //	cout << "Add gaps before: " << gaps_before << endl;
	gaps_before -= add_seq_before.length();
        newseqDNA = faststring('-', gaps_before);

        // Debug code:
        add_seq_before.ToLower();

        if (add_seq_before.size() > 0)
          newseqDNA += add_seq_before;

        unsigned DNA_pos   = start_in_seqDNA;
	unsigned amino_pos = vul.query_start; // Only needed to store insertions with respect to reference.
        unsigned count     = gaps_before + add_seq_before.length();

        // For all vulgar features: Usually only one: The M feature.
        for (unsigned i=0; i < vul.attributes.size(); ++i)
        {
	  if (global_relative_score_threshold && vul.relative_score() < global_relative_score_threshold)
	  {
	    cerr << "NOTE: Exonerate hit skipped due to low relative alignment score: " << vul.relative_score()
		 << " " << vul.get_target() << endl;
	  }
          //	  cout << "Working on vulgar with target_strand: " << vul.target_strand << endl;
          else if (vul.attributes[i].first() == 'M')
          {
            //	    cout << "Working on attribute i: " << i << " with " << endl;
            if (!vul.is_revcomp())
            {
              //	      cout << "non-revcomp" << endl;
              newseqDNA += seqDNA.substr(DNA_pos, vul.attributes[i].third());
              DNA_pos   += vul.attributes[i].third();
              count     += vul.attributes[i].third();
	      amino_pos += vul.attributes[i].second();
            }
            else
            {
              //	      cout << "revcomp" << endl;
              //	      CDnaString revcompDNA = setToReverseComplementOf(seqDNA);
              newseqDNA += seqDNA.substr(DNA_pos, vul.attributes[i].third());
              DNA_pos   += vul.attributes[i].third();
              count += vul.attributes[i].third();
	      amino_pos += vul.attributes[i].second();
            }
          }
          else if (vul.attributes[i].first() == 'F')
          {
	    // Handling frameshift alignments is difficult. Often a 'F' and 'G' attribute appear together and can only
	    // be interpreted together.
	    
            if (global_include_frameshift_alignments)
            {
              // Frameshift rules:
              // (I)   If there is no G left or right, then the F nucleotides have to be removed with respect to the AA sequence.
              // (II)  If there is a G left or right, then the F nucleotides have to be added but the codon is not complete.
              // ***** There is a remaining problem: How do we know whether the nucleotides are left or right
              // ***** justified in the frameshift codon?
              // ***** In the case (II) we have a G to the left or to the right. This might indicate the justification of the nucleotides.
              
              // if (vul.attributes[i].third() <= 3)
              {
                cout << "Hello. F" << endl;
                //	    cout << "Working on attribute i: " << i << " with " << endl;
                
                if (!vul.is_revcomp())
                {
                  if (vul.attributes[i].second() > 0)
                  {
		    amino_pos += vul.attributes[i].second();
                    ++num_WARNINGS;
                    if (global_verbosity >= 1)
                      cerr << "WARNING: Unhandled rare case: F x y with x > 0." << endl;
                  }
                  else // second() == 0, => third() > 0
                  {
                    ++num_WARNINGS;
                    if (global_verbosity >= 1)
                      cerr << "NOTE: Partial codons removed. F 0 x." << endl;
                    DNA_pos   += vul.attributes[i].third();
                  }
                }
                else
                {
                  if (vul.attributes[i].second() > 0)
                  {
		    amino_pos += vul.attributes[i].second();
                    ++num_WARNINGS;
                    if (global_verbosity >= 1)
                      cerr << "WARNING: Unhandled rare case: F x y with x > 0." << endl;
                  }
                  else // second() == 0, => third() > 0
                  {
                    ++num_WARNINGS;
                    if (global_verbosity >= 1)
                      cerr << "NOTE: Partial codons removed. F 0 x." << endl;
                    DNA_pos   += vul.attributes[i].third();
                  }
                }
              }
            }
          }
          else if (vul.attributes[i].first() == 'G')
          {
	    // Hits with gaps and without global_include_gap_alignments are not conidered at all, so we do not need this check here.
	    //            if (global_include_gap_alignments)
            {
	      // Gap in nucleotide sequence:
              if (vul.attributes[i].second() > 0)
              {
		unsigned gap_chars = vul.attributes[i].second()*3;
                newseqDNA += faststring('~', gap_chars);
                count += gap_chars;
		add_or_count(map_of_gap_sites_in_queries, make_pair(3*amino_pos, 3*vul.attributes[i].second() ));
		amino_pos += vul.attributes[i].second();
              }
              else if (vul.attributes[i].third() > 0)
              {
                // This case is difficult to handle. It requires to add gaps in the reference sequence.
                // Currently, this case is handled by removing the additional bases in the sequence added to the alignment.
                
                ++num_WARNINGS;
                if (global_verbosity >= 1)
                  cout << "WARNING: Rare gap case: G 0 x. Additional bases are skipped." << endl;
                DNA_pos   += vul.attributes[i].third();
		add_or_count(map_of_insertion_sites_suggested_in_reference, make_pair(3*amino_pos, 3*vul.attributes[i].third()));
              }
            }
          }
          else
          {
            cerr << "ERROR: Unexpected attribute in vulgar string: " << vul.attributes[i].first() << endl;
            cerr << "Normally, a \'grep \" " << vul.attributes[i].first()
		 << " \" vulgar-file-name \' should find the offending attribute in the vulgar file. " << endl;
            exit(-21);
          }
        } // END for (unsigned i=0; i < vul.attributes.size(); ++i)
        
        if (query_prot_length*3 < count)
        {
          cerr << "ERROR: count is larger than length of target in nucleotides. "
	          "This internal error should be reported.\n";
          exit(-131);
        }
        
        // Now work on gaps after and add - after:
        unsigned gaps_after = query_prot_length*3 - count;
        
        if (global_verbosity >= 40)
          cout << "DEBUG: " << key << " gaps_after: " << gaps_after << endl;
        
        if (global_num_bp_beyond_exonerate_alignment_if_at_start_or_end > 0)
        {
          if (gaps_after > 0)
          {
            if (global_verbosity >= 1000)
            {
              cout << endl;
              cout << "Sequence:   " << seq_count << endl;
              cout << "gaps_after: " << gaps_after << endl;
            }
            if (vul.is_revcomp())
            {
              // Default values for start and length of extracted region:
              unsigned add_length_end = length_this_nuc_seq - end_in_seqDNA;
		//		vul.get_target_end();
              // Where in seqDNA do we start to extract the extra bases
              unsigned new_start_in_seqDNA = end_in_seqDNA;
              
              if (global_verbosity >= 1000)
                cout << "Number potential_add_at_end_of_inserted_revcomp_seq: " << add_length_end << endl;

              if (add_length_end > 0) // Could be 0, and nothing should be done
              {
                if (add_length_end > gaps_after)
                {
                  add_length_end = gaps_after;
                }
                if (add_length_end > global_num_bp_beyond_exonerate_alignment_if_at_start_or_end)
                {
                  add_length_end= global_num_bp_beyond_exonerate_alignment_if_at_start_or_end;
                }

		if (add_length_end > 0)
		  add_seq_after = seqDNA.substr(new_start_in_seqDNA, add_length_end);

                //		CDnaString tmp1 = theSeq->getPartialSeq(0, add_at_end);
                CDnaString tmp2;
                //      tmp2.setToReverseComplementOf(tmp1);
                //		add_partial_condons_after = tmp2.c_str();
                //		faststring seq_before_start = add_partial_condons_before = theSeq->getPartialSeq(0, vul.get_target_start());
                
                //		cout << "Add after: " << add_partial_condons_after << endl;
                
              } // END if (len_at_start_of_seq > 0 && (len_at_start_of_seq <= global_num_bp_beyond_exonerate_alignment_if_at_end) )
            }
            else  // Not vul.is_revcomp()
            {
              unsigned add_length_end = length_this_nuc_seq - end_in_seqDNA;
	      // length_this_nuc_seq - vul.get_target_end();
              unsigned new_start_in_seqDNA = end_in_seqDNA;
              
              if (global_verbosity >= 1000)
                cout << "Number add_length_end: " << add_length_end << endl;

              if (add_length_end > 0) // // Could be 0, and nothing should be done
              {
                if (add_length_end > gaps_after)
                {
                  add_length_end = gaps_after;
                }
                if (add_length_end > global_num_bp_beyond_exonerate_alignment_if_at_start_or_end)
                {
                  add_length_end = global_num_bp_beyond_exonerate_alignment_if_at_start_or_end;
                } // else nothing has to be done.
                
                add_seq_after = seqDNA.substr(new_start_in_seqDNA, add_length_end);
              }
            } // END Not vul.is_revcomp()

            if (global_verbosity >= 1000)
              cout << "Add at back \"" << add_seq_after << "\" might have length 0." << endl;
            
            if (global_verbosity >= 50)
              cout << "Add back length: " << add_seq_after.length() << " " << key << endl;
            
            gaps_after -= add_seq_after.length();
          } // END  if (gaps_after > 0)
        } // END if (global_num_bp_beyond_exonerate_alignment_if_at_start_or_end > 0)

	// Debug code:
	add_seq_after.ToLower();

        //	cout << "Add gaps after: " << gaps_after << endl;
        if (add_seq_after.length() > 0)
          newseqDNA += add_seq_after;
        
        newseqDNA += faststring('-', gaps_after);
        
        //	cout << "gaps_before count : " << gaps_before << " " << count << endl;
        
        seqs_DNA_result.add_seq_to_alignment(CSequence_Mol::dna, key, newseqDNA, 'N');
        os << ">" << key << "\n" << newseqDNA << "\n";
        
      } // else if (!consider_this_read), i.e. if (consider_this_read)
      
      ++l_it;
    } // END  while (l_it !=  l_keys.end()) // Loop over all "vulgar strings", i.e. "reads with hit" and align them.
    os.close();
    
    if (global_verbosity >= 1)
    {
      cout << "\n\n";
      cout << "Some stats for the exonerate results:\n";
      cout.setf(ios::left);
      unsigned s1 = list_exonerate_results.size();
      cout << "Number of successful exonerate alignments (includes skipped hits/reads)\n"
      << setw(50)  << "to the amino acid sequence found in input file: " << s1 << endl;
      unsigned s2 = map_exonerate_count_results.size();
      cout << "Number of successful exonerate alignments without skipped hits.\n"
      << setw(50) << "Normally this is the number of aligned reads: " << s2-count_not_considered << endl;
      cout << setw(50) << "Number of sequences in result file: " << seqs_DNA_result.GetTaxaNum() << endl;

      cout << "-------------------------------------------\n";
      if (global_gap_frameshift_mode == 0)
	cout << setw(50) << "# skipped reads due to gappy alignment: " << skipped_G << endl;
      if (global_gap_frameshift_mode == 2)
      cout << setw(50) << "# skipped reads due to no gaps in aln.  : " << skipped_no_G << endl;
      cout << setw(50) << "# skipped reads due to multiple hits:   " << skipped_double_vulgar  << endl;
      cout << setw(50) << "# skipped reads due to low rel. score:  " << skipped_relative_score << endl;

      cout << "Gap insertion sites and lengths suggested from reads for the reference:\n";
      print_Mymap(cout, map_of_insertion_sites_suggested_in_reference);

      cout << "Gap sites in query sequences:\n";
      print_Mymap(cout, map_of_gap_sites_in_queries);
      cout << '\n';

      unsigned **coverage_profile = seqs_DNA_result.get_DNA_site_profile();
      print_DNA_profile(cout, coverage_profile, seqs_DNA_result.GetPosNum());
    }

    if (!global_consensus_sequence_output_filename.empty())
    {
      faststring consensus_seq;
      unsigned   replaced_tilde_gaps = 0;
      vector<unsigned> coverage_vec;

      seqs_DNA_result.WeightedConsensusSequence_DNA(consensus_seq, global_consensus_threshold,
						    global_minimum_seq_coverage_uppercase, global_minimum_seq_coverage_total,
						    &coverage_vec, global_ends_width, global_weight_fraction_in_ends,
						    'N');
      // Default: Internal gaps are encoded as tilde.
      if (global_report_gap_mode == 0)
      {
	replaced_tilde_gaps = consensus_seq.replace_char_count('~', '-');      
      }

      if (global_report_gap_mode == 3)
      {
	consensus_seq.removeSpacesBack("-");
	consensus_seq.removeSpaces("-~");
      }

      /*
      if (global_verbosity >= 50)
      {
	cerr << "Coverage values obtained while determining the consensus sequence:\n";
	for (unsigned i=0; i<coverage_vec.size(); ++i)
	  cerr << i << "\t" << coverage_vec[i] << endl;
	cerr << "=====================\n";
      }
      */

      ofstream os_con(global_consensus_sequence_output_filename.c_str());
      if (os_con.fail() )
      {
        cerr << "ERROR when trying to write the consensus file. The output file could not be opened. "
	        "Please make sure the specified path exists!" << endl;
        exit(-23);
      }
      os_con << ">Consensus" << global_input_dna_fasta_file << endl;
      os_con << consensus_seq << endl;
      os_con.close();

      if (global_verbosity >= 1)
      {
	unsigned mi, ma;
	double mean, sd, Q1, Q2, Q3;
	vec_min_max(coverage_vec, mi, ma);
	vec_mean_sd(coverage_vec, mean, sd);
	// Be careful: coverage_vec will be sorted after calling this function:
	vec_median_quartile_sort_method3(coverage_vec, Q1, Q2, Q3);

	cout.precision(2);
	cout.setf(ios::fixed);
	cout << "-------------------------------------------\n";
	cout << setw(50) << "Coverage minimum: " << mi   << "\n";
	cout << setw(50) << "Coverage maximum: " << ma   << "\n";
	cout << setw(50) << "Coverage mean:    " << mean << "\n";
	cout << setw(50) << "Coverage median:  " << Q2   << "\n";
	cout << "-------------------------------------------\n";
      }
    }
    
  } // END: After reading the data, this is the main block.
  
  if (num_WARNINGS > 0)
  {
    if (global_verbosity >= 1)
      cout << "There have been " << num_WARNINGS << " warnings or notes. See above for details." << endl;
    else
      cout << "Due to the verbosity parameter set to 0, you missed " << num_WARNINGS << " warnings or notes. Increase the verbosity if you want to see the details and run again." << endl;
  }
  
  cout << PROGNAME " finished successfully." << endl << endl;
} // END: main(...)
