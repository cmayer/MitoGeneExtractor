#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include <ctime>
#include <iomanip>
#include "faststring2.h"
#include "CSequences2.h"
#include "CSequence_Mol2_1.h"
#include "Ctriple.h"
#include "CFile/CFile2_1.h"
#include "CDnaString2.h"
#include "global-types-and-parameters_MitoGeneExtractor.h"
#include "statistic_functions.h"
#include "fastq.h"
#include "Cfastq-sequences.h"


/// #include <utility>

using namespace std;

typedef pair<unsigned, unsigned> Key;
typedef map< Key , unsigned> Mymap;


// Verbosity rules:
//      verbosity 0 -> only error messages that lead to exit() (goes to cerr).
//                     Welcome is printed to cout. Parameters are printed to cout.
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
//unsigned  length_constaint              = 2000;

unsigned num_WARNINGS = 0;

void append_fasta_files_to_exonerate_input_file(FILE *ofp, vector<string> global_input_dna_fasta_filenames)
{
  for (vector<string>::size_type i=0; i < global_input_dna_fasta_filenames.size(); ++i)
  {
    FILE* ifp; // input file pointer
    unsigned nread;

    ifp = fopen(global_input_dna_fasta_filenames[i].c_str(), "r");

    if (!ifp)
    {
      cerr << "The fasta input file " << global_input_dna_fasta_filenames[i] << "\ncould not be opened. "
	      "Maybe it does not exist. Exiting.\n";
      good_bye_and_exit(-27);
    }
    
    char const_buf [1000000];
    while ((nread = fread(const_buf, sizeof(char), sizeof(const_buf), ifp)) > 0)
    {
      fwrite(const_buf, sizeof(char), nread, ofp);
    }
    fclose(ifp);
  }
}

void append_fastq_files_to_exonerate_input_file(FILE *ofp, vector<string> global_input_dna_fastq_filenames)
{
  for (vector<string>::size_type i=0; i < global_input_dna_fastq_filenames.size(); ++i)
  {
    CFile      infile;
    infile.ffopen(global_input_dna_fastq_filenames[i]);

    if (infile.fail())
    {
      cerr << "The fastq input file " << global_input_dna_fastq_filenames[i] << "\ncould not be opened. "
	      "Maybe it does not exist. Exiting.\n";
	good_bye_and_exit(-27);
    }

    fastq_sequences fq_in_collection;
    fq_in_collection.read_fastq(infile, global_input_dna_fastq_filenames[i].c_str(), 1);
    infile.close();

    /*
      unsigned N_seq_in_fq_file = fq_in_collection.size();

      if (global_verbosity >= 1)
      {
         cout << "Found " << N_seq_in_fq_file << " sequences in input file " << global_input_dna_fastq_filenames[i]
              << endl;
      }
    */

    //    fq_in.print(cout);

    CSequences2 *pseqs = new CSequences2(CSequence_Mol::dna);
    fq_in_collection.add_sequences_to_CSequences_object(pseqs);
    pseqs->ExportSequences(ofp, 'f', UINT_MAX);
    delete pseqs;
  }
}


struct stats_for_given_target
{
  // All numbers will be for this specific target
  unsigned skipped_double_vulgar;
  unsigned skipped_G;
  unsigned skipped_F;
  unsigned skipped_relative_score;
  unsigned skipped_no_G;
  unsigned not_considered;
  unsigned not_best_score_for_query;
  unsigned good_hits_added_to_alignment;
  unsigned hits_in_vulgar_file;

  stats_for_given_target():skipped_double_vulgar(0), skipped_G(0),
			   skipped_F(0), skipped_relative_score(0),
			   skipped_no_G(0), not_considered(0),
			   not_best_score_for_query(0),
			   good_hits_added_to_alignment(0),
			   hits_in_vulgar_file(0)
  {}

  void increment_skipped_double_vulgar()
  {
    ++skipped_double_vulgar;
  }

  void increment_skipped_G()
  {
    ++skipped_G;
  }

  void increment_skipped_F()
  {
    ++skipped_F;
  }

  void increment_skipped_relative_score()
  {
    ++skipped_relative_score;
  }

  void increment_skipped_no_G()
  {
    ++skipped_no_G;
  }

  void increment_not_considered()
  {
    ++not_considered;
  }

  void increment_not_best_score()
  {
    ++not_best_score_for_query;
  }

  void increment_good_hits_added_to_alignment()
  {
    ++good_hits_added_to_alignment;
  }

  void increment_hits_in_vulgar_file()
  {
    ++hits_in_vulgar_file;
  }
};

/*
void copy_file(const char * sourcename, const char * destname, const char * out_mode)
{
  FILE* ifp; // input file pointer                                                                                   
  FILE* ofp; // output file pointer                                                                                  

  unsigned nread;
  ifp = fopen(sourcename, "r");
  ofp = fopen(destname, out_mode);

  char const_buf [100000];
  while ((nread = fread(const_buf, sizeof(char), sizeof(const_buf), ifp)) > 0)
  {
    fwrite(const_buf, sizeof(char), nread, ofp);
  }

  fclose(ifp);
  fclose(ofp);
}
*/


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

// IMPORTANT: Reads/sequences are targets.
//            Reference sequences are queries.

class vulgar
{
public:
  static map<faststring, short>  queryID2Index;

  faststring queryID;          // In this program, this is the AA-COI consensus seq. exonerate aligned against.
  short      queryIndex;       // Every query sequence (reference sequence) has a query number.
                               // This avoids many lookups for index.
  unsigned   query_start;      // 0 based numbers.
  unsigned   query_end;        // The first position after the specified range.
  char       query_strand;
  faststring targetID;         // In this program, they are the reads.
  unsigned   target_start;     // 0 based numbers.
  unsigned   target_end;       // The first position after the specified range.
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

  public:
  static short add_query_to_query2short_map_or_get_index(faststring q)
  {
    map<faststring, short>::iterator find_it = queryID2Index.find(q);
    if (find_it != queryID2Index.end())
      return find_it->second;

    unsigned s = queryID2Index.size();
    queryID2Index[q] = s;
    return s;
  }

  short get_index_in_query2short_map(faststring q)
  {
    map<faststring, short>::iterator find_it = queryID2Index.find(q);
    if (find_it != queryID2Index.end())
      return find_it->second;
    else
    {
      cerr << "Reference sequence appears in Exonerate output, but not in "
    	      "reference sequences. This should not be possible.\n";
      good_bye_and_exit(-127);
      return 0; // Only introduced to silence warnings about no return for non-void function.
    }
  }

  private:
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
      good_bye_and_exit(-4);
    }

    if (l[0] != "vulgar:")
    {
      cerr << "ERROR: Bad vulgar string: The vulgar string is expected to start with \"vulgar:\" but found: " << str << endl;
      good_bye_and_exit(-4);
    }

    // It is a convention in exonerate that the protein sequence is the query and the DNA sequence the target:
    queryID       = l[1];
    query_start   = l[2].ToUnsigned();
    query_end     = l[3].ToUnsigned();
    query_strand  = l[4][0];

    targetID      = l[5];
    target_start  = l[6].ToUnsigned();
    target_end    = l[7].ToUnsigned();
    target_strand = l[8][0];

    score         = l[9].ToUnsigned();

    queryIndex    = get_index_in_query2short_map(queryID);
    
    int m = l.size()-10;
    if (m%3 != 0)
    {
      cerr << "ERROR: Bad vulgar string encountered. The number of fields must "
	      "by 10+3*n where n is an integer. But the number is "
           << m << "." << endl;
      cerr << "       Vulgar string: " << str << endl;
      good_bye_and_exit(-5);
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
  void print(ostream &os) const
  {
    os << "vulgar: "
    << queryID       << " "
      //    << queryIndex    << " "
    << query_start   << " "
    << query_end     << " "
    << query_strand  << " "
    << targetID      << " "
    << target_start  << " "
    << target_end    << " "
    << target_strand << " "
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

  const faststring & get_queryID() const
  {
    return queryID;
  }

  short get_queryIndex() const
  {
    return queryIndex;
  }
  
  const faststring & get_targetID() const
  {
    return targetID;
  }

  faststring get_hitkey() const
  {
    return queryID + "##" + targetID;
  }
  
  bool has_F() const
  {
    return bool_has_F;
  }
  
  bool has_G() const
  {
    return bool_has_G;
  }
  
  bool has_M() const
  {
    return bool_has_M;
  }
  
  bool has_non_M() const
  {
    return bool_has_non_M;
  }
  
  unsigned num_attributes() const
  {
    return attributes.size();
  }
  
  unsigned get_query_start() const
  {
    return query_start;
  }

  unsigned get_target_start() const
  {
    return target_start;
  }
  
  unsigned get_target_end() const
  {
    return target_end;
  }
  
  bool is_revcomp() const
  {
    return (target_strand == '-');
  }

  // Divides the exonerate score by the length of the hit and returns the result.
  double relative_score() const
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

  unsigned get_score()
  {
    return score;
  }
};


// Global variable definition for static member variable of vulgar class:
map<faststring, short>  vulgar::queryID2Index;


bool fileExists(const char *filename)
{
  ifstream is(filename);
  if (is.fail())
    return false;
  else
    return true;
}

inline unsigned size_keyset_multimap(multimap<faststring, vulgar *> &mm)
{
  set<pair<faststring, faststring> > s;
  multimap<faststring, vulgar *>::iterator it = mm.begin();
  while (it != mm.end())
  {
    vulgar &vul = *(it->second);
    s.insert(make_pair(vul.queryID, vul.targetID));
    ++it;
  }
  return s.size();
}

inline unsigned number_of_entries_targetID(multimap<faststring, vulgar *> &mm, faststring targetID)
{
  return mm.count(targetID);
}

inline unsigned number_of_entries_targetID_and_queryID(multimap<faststring, vulgar *> &mm,
						       faststring targetID,
						       faststring &quereyID)
{
  pair<multimap<faststring, vulgar *>::iterator, multimap<faststring, vulgar *>::iterator> it_bounds_pair;
  multimap<faststring, vulgar *>::iterator it, it_end;

  it_bounds_pair = mm.equal_range(targetID);
  it     = it_bounds_pair.first;
  it_end = it_bounds_pair.second;
  unsigned count = 0;
  
  while (it != it_end)
  {
    if ( (*(it->second)).get_queryID() == quereyID )
      ++count;
    ++it;
  }
  return count;
}


  
// Determines the aligned sequence newseqDNA
void determine_alignment_string(CSequence_Mol* theSeq, const vulgar &vul, unsigned query_prot_length,
				CDnaString &newseqDNA, unsigned seq_count,
                                map<pair<unsigned, unsigned>, unsigned> &map_of_gap_sites_in_queries,
				map<pair<unsigned, unsigned>, unsigned> &map_of_insertion_sites_suggested_in_reference)
{
  // The query prot sequence start coordinate sets the number of gaps before:
  unsigned    gaps_before = vul.query_start*3;
  unsigned    start_in_seqDNA, end_in_seqDNA;
  CDnaString  seqDNA;


  if (global_verbosity >= 100)
  {
    cerr << "determine_alignment_string was called with the following parameters:\n"
	 << "query_prot_length:  " << query_prot_length << endl;
    cerr << "Reference (query):  " << vul.get_queryID() << endl;
    cerr << "Reference (target): " << vul.get_targetID() << endl;
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
   good_bye_and_exit(-19);
   }
   */
  
  // Problem: exonerate only aligns complete codons to amino acids
  
  // Test for consistent bases in partial codons:
  
  // Exonerate aligns the reads to the amino acid sequence. Only complete codons are considered as alignment matches.
  // So using the coordinates obtained from exonerate, a partial codon is removed in 2/3 of the cases at the beginning as well as at the end.
  //
  
  // We might want to add unmatched bp before and after, e.g. for partial codons or in order to see how much they differ from other sequences.
  // In order not to intervene with multiple attributes, we treat that case
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
        
        if (start_in_seqDNA > 0) // Could be 0, and nothing should be done
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
          cerr << "Sequence name: " << theSeq->getName() << endl;
          cerr << "Starting coordinate of + strand: " << new_start_in_seqDNA  << endl;
          cerr << "Length:                          " << add_length_beginning << endl;
          cerr << "String:                          " << add_seq_before       << endl;
          
          //            faststring seq_before_start = add_partial_condons_before;
          
          //       cout << endl;
          //       cout << "key: " << key << endl;
          //       seq_before_start += "|";
          //       seq_before_start += theSeq->getPartialSeq(vul.get_target_start(), 12);
          //       faststring aa_at_seq_before_start =     seqs_target.get_seq_by_index(0)->getPartialSeq(vul.query_start-1, 5);
          //       cout << "gaps_before: " << gaps_before << " vul.get_target_start(): "   << vul.get_target_start() << endl;
          //       cout << "Sequence before start: " << seq_before_start << endl;
          //       cout << "AA at this posiiton:   " << aa_at_seq_before_start << endl;
        }
      } // END Not revcomp
      
      if (global_verbosity >= 1000)
        cout << "Add partial codon front \"" << add_seq_before << "\" might have length 0." << endl;
      
      if (global_verbosity >= 50)
        cout << "Add front length: " << add_seq_before.length() << " " << theSeq->getName() << endl;
      
    }  // END if (gaps_before > 0)
  } // END if (global_num_bp_beyond_exonerate_alignment_if_at_end > 0)
  
  //   cout << "Add gaps before: " << gaps_before << endl;
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
    /* These hits should never be handled here:
    if (global_relative_score_threshold && vul.relative_score() < global_relative_score_threshold)
    {
      cerr << "NOTE: Exonerate hit skipped due to low relative alignment score: "
	   << vul.relative_score()
	   << " " << vul.get_targetID() << endl;
    }
    //     cout << "Working on vulgar with target_strand: " << vul.target_strand << endl;
    else */
    if (vul.attributes[i].first() == 'M')
    {
      //       cout << "Working on attribute i: " << i << " with " << endl;
      if (!vul.is_revcomp())
      {
        //         cout << "non-revcomp" << endl;
        newseqDNA += seqDNA.substr(DNA_pos, vul.attributes[i].third());
        DNA_pos   += vul.attributes[i].third();
        count     += vul.attributes[i].third();
        amino_pos += vul.attributes[i].second();
      }
      else
      {
        //         cout << "revcomp" << endl;
        //         CDnaString revcompDNA = setToReverseComplementOf(seqDNA);
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
          //       cout << "Working on attribute i: " << i << " with " << endl;
          
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
      // Hits with gaps and without global_include_gap_alignments are not considered at all, so we do not need this check here.
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
      good_bye_and_exit(-21);
    }
  } // END for (unsigned i=0; i < vul.attributes.size(); ++i)
  
  if (query_prot_length*3 < count)
  {
    cerr << "ERROR: count is larger than length of target in nucleotides. "
            "This internal error should be reported.\n"
            "Vulgar hit for this problematic sequence:\n";
    vul.print(cerr);
    cerr << endl;
    cerr << "query_prot_length   " << query_prot_length   << endl;
    cerr << "query_prot_length*3 " << query_prot_length*3 << endl;
    cerr << "count (nuc-length)  " << count  << endl;
    good_bye_and_exit(-131);
  }
  
  // Now work on gaps after and add - after:
  unsigned gaps_after = query_prot_length*3 - count;
  
  if (global_verbosity >= 40)
    cout << "DEBUG: " << theSeq->getName() << " gaps_after: " << gaps_after << endl;
  
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
        //      vul.get_target_end();
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
          
          //      CDnaString tmp1 = theSeq->getPartialSeq(0, add_at_end);
          CDnaString tmp2;
          //      tmp2.setToReverseComplementOf(tmp1);
          //      add_partial_condons_after = tmp2.c_str();
          //      faststring seq_before_start = add_partial_condons_before = theSeq->getPartialSeq(0, vul.get_target_start());
          
          //      cout << "Add after: " << add_partial_condons_after << endl;
          
        } // END if (len_at_start_of_seq > 0 && (len_at_start_of_seq <= global_num_bp_beyond_exonerate_alignment_if_at_end) )
      }
      else  // Not vul.is_revcomp()
      {
        unsigned add_length_end = length_this_nuc_seq - end_in_seqDNA;
        // length_this_nuc_seq - vul.get_target_end();
        unsigned new_start_in_seqDNA = end_in_seqDNA;
        
        if (global_verbosity >= 1000)
          cout << "Number add_length_end: " << add_length_end << endl;
        
        if (add_length_end > 0) // Could be 0, and nothing should be done
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
        cout << "Add back length: " << add_seq_after.length() << " " << theSeq->getName() << endl;
      
      gaps_after -= add_seq_after.length();
    } // END  if (gaps_after > 0)
  } // END if (global_num_bp_beyond_exonerate_alignment_if_at_start_or_end > 0)
  
  // Debug code:
  add_seq_after.ToLower();
  
  //   cout << "Add gaps after: " << gaps_after << endl;
  if (add_seq_after.length() > 0)
    newseqDNA += add_seq_after;
  
  newseqDNA += faststring('-', gaps_after);
  
  //   cout << "gaps_before count : " << gaps_before << " " << count << endl;
}


int run_exonerate(const faststring &exonerate_binary, const faststring &protein_reference, const faststring &fasta_reads, const faststring &vulgar_base_file, char genetic_code_number, int frameshift_penalty, int max_trials, char mode=0)
{
  faststring code_string;

  if (genetic_code_number == 1) // The Standard Code (transl_table=1)
    code_string = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 2) // The Vertebrate Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
  else if (genetic_code_number == 3)  // The Yeast Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 4)  // The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 5)  // The Invertebrate Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 6)  // The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
    code_string = "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";


  else if (genetic_code_number == 9)  // The Echinoderm and Flatworm Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 10) // The Euplotid Nuclear Code (transl_table=10)
    code_string = "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 11) // The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
    code_string = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 12) // The Alternative Yeast Nuclear Code (transl_table=12)
    code_string = "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 13) // The Ascidian Mitochondrial Code
    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 14) // The Alternative Flatworm Mitochondrial Code
    code_string = "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
  else if (genetic_code_number == 16) // Chlorophycean Mitochondrial Code
    code_string = "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  else
  {
    cerr << "ERROR: Unknown genetic code number passed to function run_exonerate. Code number: " << genetic_code_number << endl;
    good_bye_and_exit(-20);
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

  bool combined_input_sequence_file_created = false;
  CSequences2 seqs_DNA_reads_target(CSequence_Mol::dna);
  CSequences2 seqs_prot_query(CSequence_Mol::protein);
  //  CSequences2 seqs_DNA_result(CSequence_Mol::dna);

  // unsigned skipped_double_vulgar=0, skipped_G=0, skipped_F=0, skipped_relative_score=0, skipped_no_G=0;

  vector<faststring>   seqnames_of_references;
  vector<CSequences2*> result_alignments_for_references;

  string global_input_dna_fasta_file;

  // Single fasta input file:
  if (global_input_dna_fasta_filenames.size() == 1 && global_input_dna_fastq_filenames.size() == 0)
  {
    global_input_dna_fasta_file = global_input_dna_fasta_filenames[0];
  }
  else
  {
    if (global_tmp_directory.find('~') != string::npos)
    {
      cerr << "The character \'~\' (tilde) is not allowed in the string "
	      "specifying the path to a temporary directory.\nExiting.\n" << endl;
      exit(-13);
    }
    
    string tmpstring = global_tmp_directory + "/Concatenated_exonerate_input_XXXXXX";
    cerr << "Filename:" << tmpstring << endl;

      
    char  *tmpstr = new char [global_tmp_directory.length() + 40];
    strcpy (tmpstr, tmpstring.c_str());
    cerr << "Filename:" << tmpstr << endl;
    int   fd = mkstemp(tmpstr);
    cerr << "Filename:" << tmpstr << endl;

    FILE  *ofp = fdopen(fd, "w");
    combined_input_sequence_file_created = true;
    if (fd == -1 || ofp == 0)
    {
      cerr << "ERROR: Could not open the temporary file for the fasta file "
	      "passed to the Exonerate program." << endl;
      good_bye_and_exit(-13);
    }

    global_input_dna_fasta_file = string(tmpstr);

    if (global_input_dna_fasta_filenames.size() > 0)
    {
      append_fasta_files_to_exonerate_input_file(ofp, global_input_dna_fasta_filenames);
    }
    
    if (global_input_dna_fastq_filenames.size() > 0)
    {
      append_fastq_files_to_exonerate_input_file(ofp, global_input_dna_fastq_filenames);
    }

    //    global_input_dna_fasta_filenames = faststring(tmpstr);
    //    cout << "Temporary name: " << global_input_dna_fasta_file << endl;
    fclose(ofp);
  }
  
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
    good_bye_and_exit(-1);
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
    good_bye_and_exit(-3);
  }

  // TODO:
  // Instead of printing this error message and exiting, we could fix this
  // problem at this point, at least if we compute the vulgar file ourselves.
  // But this should be an extremely  rare case. NGS reads should always be
  // unique. I do not see the scenario in which this would occur.

  if (!seqs_DNA_reads_target.are_short_names_unique())
  {
    cerr << "ERROR: The sequence names in the DNA read file are not unique when "
            "trimming them at the first space, "
            "which is done by many programs, including exonerate. Please verify your "
            "input data. Often it "
            "helps to rename sequences by replacing spaces with e.g. underscores.\n";
    good_bye_and_exit(-3);
  }
  
  if (global_verbosity >= 3)
  {
    cout << "Finished reading target read sequences. Found this number of "
            "target sequences: "
	 << seqs_DNA_reads_target.GetTaxaNum()  << endl;
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
    good_bye_and_exit(-1);
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
    good_bye_and_exit(-3);
  }

  if (!seqs_prot_query.are_short_names_unique())
  {
    cerr << "ERROR: The sequence in the reference file are not unique."
            "Please rename the reference sequences to ensure they are unique.\n";
    good_bye_and_exit(-3);
  }

  seqs_prot_query.get_short_sequence_names(seqnames_of_references);

  // Create vector of output objects:
  for (std::vector<faststring>::size_type i=0; i < seqnames_of_references.size(); ++i)
  {
    result_alignments_for_references.push_back(new CSequences2(CSequence_Mol::dna) );
    vulgar::add_query_to_query2short_map_or_get_index(seqnames_of_references[i]);
  }
    
  if (global_verbosity >= 100)
  {
    cerr << "Sequence names found in reference file:\n";
    for (std::vector<faststring>::size_type i=0; i<seqnames_of_references.size(); ++i)
      cerr << seqnames_of_references[i] << endl;
  }
  
  // The program now accepts multiple reference sequences.
  /*
  if (seqs_prot_query.GetTaxaNum() != 1)
  {
    cerr << "ERROR: Found multiple sequences in the query protein file. Only one sequence is allowed in this program." << endl;
    good_bye_and_exit(-23);
  }
  */

  //  unsigned query_prot_length = seqs_prot_query.GetPosNum();

  map<faststring, unsigned> map_query_prot_lengths;
  for (std::vector<faststring>::size_type i=0; i < seqnames_of_references.size(); ++i)
    map_query_prot_lengths[seqnames_of_references[i]] = seqs_prot_query.get_seq(i)->length();

  if (global_verbosity >= 3)
  {
    cout << "Finished reading fasta file with protein reference sequence." << endl;
  }

  //**************************
  // Export query if sequence names had to be modified.
  //**************************
  
  bool     vulgar_file_read_successfully = false;
  unsigned count_appempts_to_read_vulgar_file = 0;

  // OLD
  //  vector<pair<faststring, vulgar *> > vec_exonerate_results;
  //  map<faststring, unsigned>       map_exonerate_count_results;
  // END OLD
  
  //*********
  vector<vulgar *>                   vec_of_hits_as_in_file;
  multimap<faststring, vulgar * >    map_of_vulgar_hits_targets_as_keys;

  vector<stats_for_given_target> vector_of_hit_stats_for_targets(seqnames_of_references.size() );

  //*********

  while (!vulgar_file_read_successfully)
  {
    // Reset vec and map and counters. They might have data from an incomplete vulgar file that we tried to read before.
    vec_of_hits_as_in_file.clear();
    map_of_vulgar_hits_targets_as_keys.clear();

    ++count_appempts_to_read_vulgar_file;
    faststring exonerate_output_file_name = global_vulgar_file.c_str();

    //*****************************
    // Run exonerate if necessary:
    //*****************************

    // Do we need to compute the vulgar file?
    if ( !fileExists(global_vulgar_file.c_str()))
    {
      if(global_exonerate_binary.empty() )
        global_exonerate_binary = "exonerate";
 
      if (global_verbosity >= 1 && count_appempts_to_read_vulgar_file == 1)
        cout << "No vulgar file with the specified name has been found. So the binary \""
	     << global_exonerate_binary << "\" will be used to create the vulgar file." << endl;

      int num_trials = run_exonerate(global_exonerate_binary.c_str(), global_input_prot_reference_sequence.c_str(),
                                     global_input_dna_fasta_file.c_str(), global_vulgar_file.c_str(),
				     global_genetic_code_number, global_frameshift_penalty, 10, global_run_mode);

      if (num_trials == 10+1)
      {
        cerr << "ERROR: Running exonerate failed. "
	        "The generated vulgar file is incomplete and should be removed manually. Exiting.\n";
	good_bye_and_exit(-1);
      }

      if (global_verbosity >= 1)
        cout << "Exonerate result file was created successfully after N = " << num_trials << " trials." << endl;
    }
    else
    {
      if (global_verbosity >= 1)
        cout << "The specified vulgar file exists, so it does not have to be recomputed with exonerate." << endl;
    }

    //*****************************
    // Parse exonerate output:
    //*****************************

    CFile exonerate_output_file;
    //  map<faststring, vulgar>  m_exonerate_results;

    exonerate_output_file.ffopen(exonerate_output_file_name.c_str());

    if (exonerate_output_file.fail())
    {
      cerr << "ERROR: Could not open file " << exonerate_output_file_name << endl;
      cerr << "Exiting.\n\n";
      good_bye_and_exit(-6);
    }

    faststring line;
    while(true) // Loop over all lines in file and create vulgar objects for all vulgar lines:
    {
      //    char c = exonerate_output_file.getchar();
      //    cout << "char: " << c << endl;

      exonerate_output_file.getline(line);

      if (exonerate_output_file.fail())
	break;

      line.removeSpacesBack();

      // These are valid hits:
      if (line.size() > 7 && line.substr(0,7) == "vulgar:")
      {
        //       cout << "Adding line: " << line << endl;
	// Create vulgar object:
        vulgar *pv = new vulgar(line);

	// Preselection of hits:
	// Depeding on the program parameters, some hits will not be considered. This is a good
	// point to remove them.
	{
	  vulgar &vul = *pv;
	  short refnumber = pv->get_queryIndex();
	  vector_of_hit_stats_for_targets[refnumber].increment_hits_in_vulgar_file();

	  // Currently, considering frameshift hits is not supported. Anyway we filter them first.
	  if (vul.has_F() && !global_include_frameshift_alignments)
	  {
	    vector_of_hit_stats_for_targets[refnumber].increment_skipped_F();
	    vector_of_hit_stats_for_targets[refnumber].increment_not_considered();
	    //	    ++skipped_F;
	    //	    ++count_not_considered;
	    if (!vul.has_G() )
	      vector_of_hit_stats_for_targets[refnumber].increment_skipped_G();
	      //	      ++skipped_G;
	    delete pv;
	    continue;
	  }

	  if ((!vul.has_G() && global_gap_frameshift_mode == 2) )
	  {
	    vector_of_hit_stats_for_targets[refnumber].increment_skipped_no_G();
	    //	    ++skipped_no_G;
	    vector_of_hit_stats_for_targets[refnumber].increment_not_considered();
	    //	    ++count_not_considered;
	    delete pv;
	    continue;
	  }

	  if ((vul.has_G() && global_gap_frameshift_mode == 0) )
	  {
	    vector_of_hit_stats_for_targets[refnumber].increment_not_considered();
	    //	    ++count_not_considered;
	    vector_of_hit_stats_for_targets[refnumber].increment_skipped_G();
	    //	    ++skipped_G;
	    delete pv;
	    continue;
	  }

	  if (global_relative_score_threshold && vul.relative_score() < global_relative_score_threshold)
	  {
	    cerr << "NOTE: Exonerate hit skipped due to low relative alignment score: " << vul.relative_score()
		 << " " << vul.get_targetID() << endl;
	    vector_of_hit_stats_for_targets[refnumber].increment_skipped_relative_score();
	    vector_of_hit_stats_for_targets[refnumber].increment_not_considered();
	    //	    ++skipped_relative_score;
	    //	    ++count_not_considered;
	    delete pv;
	    continue;
	  }
	} // END Preselection of hits

	vec_of_hits_as_in_file.push_back(pv);
	map_of_vulgar_hits_targets_as_keys.emplace(pv->get_targetID(),pv);

	////        vec_exonerate_results.push_back(make_pair(pv->get_targetID(),pv));
	////        faststring hit_key = pv->get_hitkey();
	////        add_or_count(map_exonerate_count_results, hit_key);
      } // END if (line.size() > 7 && line.substr(0,7) == "vulgar:")
      else  // These lines do not contain hits:
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
        if (line.empty()) // Skip empty lines even though they normally do not occur in vulgar files.
          continue;
        cerr << "WARNING: Ignoring non standard line in vulgar file: \"" << line << "\""<< endl;
      }
    } // END  while(true) // Reading all lines of the vulgar file.

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


  // Remove exonerate input file if this was created in this program
  if (combined_input_sequence_file_created)
  {
    if (!global_keep_concatenated_input_file)
    {
      remove(global_input_dna_fasta_file.c_str());
    }
    else
    {
      string final_concat_name = global_input_dna_fasta_file + ".fas";
      rename(global_input_dna_fasta_file.c_str(), final_concat_name.c_str());
      if (global_verbosity >= 1)
      {
	cout << "Concatenated sequence file written to: "
	     << global_input_dna_fasta_file.c_str() << ".fas" << endl;
      }
    }
  }

  // Remove the vulgar file, if the user has not specified a file name to save the file permanently.
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
    unsigned s1 = vec_of_hits_as_in_file.size();
    cout << "Number of exonerate hits (including multiple hits) to the amino "
            "acid sequence\nfound in input file:, without skipped hits: " << s1 << endl;
    //    cout << "Number of hits that have been filtered (gaps, frameshift, relative score): " << count_not_considered << endl;
    unsigned s2 = map_of_vulgar_hits_targets_as_keys.size();
    // TODO: Check what is the difference to s1. The text has been changed. I think s1 used to be without double hits?
    //       If s1 == s2, we can remove one line or indicate.
    cout << "Number of successful exonerate alignments with multiple hits: " << s2 << endl;
    //    cout << "Number of queryID, targetID hit pairs: " << size_query_target_set(map_of_vulgar_hits_targets_as_keys);
  }

  // Debug output:
  if (global_verbosity >= 100)
  {
    cout << "DEBUG OUTPUT: Exonerate results.\n";

    vector<vulgar *>::iterator it;
    it = vec_of_hits_as_in_file.begin();
    while (it !=  vec_of_hits_as_in_file.end())
    {
      faststring target_key = (**it).get_targetID();
      cout << " => item " << target_key << " ";
      (**it).print(cout);
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

  /*
  ofstream os;
  os.open(global_alignment_output_file.c_str());
  */  

  // After reading the data, this is the main analysis and alignment block.
  {
    multimap<faststring, vulgar *>::iterator seqs_it, equal_range_start, equal_range_end, tmp_it;
    //    pair<multimap<faststring, vulgar *>::iterator,multimap<faststring, vulgar *>::iterator> equal_range_its;

    
    //    vector<vulgar*>::iterator vec_it;
    //    map<faststring, vulgar>::iterator m_it;
    
    unsigned seq_count = 0;
    
    //**************************
    // For each sequence, use exonerate output to align the sequences:
    //**************************
    // Loop over all vulgar strings and align the corresponding read:
    // Sequences without vulgar string are ignored. - Diagnostic output
    // could be written for these sequences.
    
    /*
    bool consider_this_read;
    unsigned count_not_considered = 0;
    */

    vector<map<pair<unsigned, unsigned>, unsigned> > all_maps_of_insertion_sites_suggested_in_reference(seqnames_of_references.size());
    vector<map<pair<unsigned, unsigned>, unsigned> > all_maps_of_gap_sites_in_queries(seqnames_of_references.size());

    //    map<pair<unsigned, unsigned>, unsigned> map_of_insertion_sites_suggested_in_reference, map_of_gap_sites_in_queries;

    // The main analysis loop over all targets (reads or sequences).
    // They are called targets since exonerate treats them as targets.

    // For all target, i.e. read sequences do:
    seqs_it = map_of_vulgar_hits_targets_as_keys.begin();
    while (seqs_it !=  map_of_vulgar_hits_targets_as_keys.end())
    {
      const faststring &target_key = seqs_it->first;

      if (global_verbosity >= 4)
      {
	cout << "Creating alignments for (read) sequence ID: " << target_key << endl;
      }

      // Determine equal range. We could call equal range, but we already
      // have the beginning, so it is more efficient to only look for the end.
      tmp_it = equal_range_start = seqs_it;
      unsigned count_hits_this_target_before_final_filtering = 0;

      // A target sequence might have hits with multiple queries (references).
      // Furthermore, the hits have to be added to the right alignment.
      // This will be taken care of in the following:

      map<faststring, unsigned> count_hits_for_different_queries_of_this_target;
      vector<vulgar *>          vec_of_hits_for_this_target_key;
      vector<vulgar *>          hits_to_use_for_target(seqnames_of_references.size(), 0);

      // Loop over all hits for this target.
      while (tmp_it != map_of_vulgar_hits_targets_as_keys.end() && tmp_it->first == target_key)
      {
	vec_of_hits_for_this_target_key.push_back(tmp_it->second);
	++count_hits_this_target_before_final_filtering;
	add_or_count(count_hits_for_different_queries_of_this_target, (*(tmp_it->second)).get_queryID() );
	++tmp_it;
      }
      equal_range_end = tmp_it;
      
      // This should be the case for maybe 99.x% of the target sequences.
      if (count_hits_this_target_before_final_filtering == 1)
      {
	if (vec_of_hits_for_this_target_key.size() != 1)
	{
	  cerr << "vec_of_hits_for_this_target_key has unexpected size. "
	          "Should be 1. Has : "
	       << vec_of_hits_for_this_target_key.size() << endl;
	}
	vulgar *p = vec_of_hits_for_this_target_key[0];
	hits_to_use_for_target[p->get_queryIndex()] = p;
      }
      else // We have multiple hits for target (read) sequences:
      {
	// We have more than one hit. This should be an exception if we
	// look for different genes.
	// This situation is expected if we have competing references for the same gene.

	for (std::vector<vulgar*>::size_type i=0; i<vec_of_hits_for_this_target_key.size(); ++i)
	{
	  vulgar *p       = vec_of_hits_for_this_target_key[i];
	  vulgar &vul     = *p;
	  short ref_index = vul.get_queryIndex();

	  // Address 0          indicates: No hit found so far.
	  // Address (void *)-1 indicates: Double hit found for this reference already.  

	  // Is this the first hit against this reference?
	  if (hits_to_use_for_target[ref_index] == NULL )
	  {
	    hits_to_use_for_target[ref_index] = p;
	  }
	  else // Found more than one hit against this target. This is to be handled.
	  {
	    if (global_include_double_hits && hits_to_use_for_target[ref_index] !=  (void *)-1 )
	    {
	      // We include double hits, but only the best one of multiple hits:
	      // Is this score better than the one of the previous hit?
	      if (vul.get_score() > hits_to_use_for_target[ref_index]->get_score() )
		hits_to_use_for_target[ref_index] = p;
	      else
		vector_of_hit_stats_for_targets[ref_index].increment_skipped_double_vulgar();
	    }
	    else // We have a double hit, but do not want to include them.
	         // So we block this reference for this target. 
	    {
	      hits_to_use_for_target[ref_index] = (vulgar *)-1; // We only add a hit if pointer is 0.
	      //	      ++skipped_double_vulgar;                        // -> So we block this reference.
	      // We skip two hits, one that has already been entered and the one we just tried to enter.
	      vector_of_hit_stats_for_targets[i].increment_skipped_double_vulgar();
	      vector_of_hit_stats_for_targets[i].increment_skipped_double_vulgar();
	    }
	  }
	}

	// Hits can be assigned to multiple sequences.
	// Typically this is the case if we have multiple variants of the same reference gene.
	// In this case we do the following:
	// If one hit has a higher score than the hits to the other references, we only use this hit.
	// If multiple hits have the same score we use all.
	// This ensures that we cover all variants but prevent hits to be included in inferior
	// references 

	// Targets with pointer (void *)-1)
	for (std::vector<vulgar*>::size_type i=0; i<hits_to_use_for_target.size(); ++i)
	{
	  if (hits_to_use_for_target[i] == (vulgar *)-1 )
	    hits_to_use_for_target[i] = NULL;
	}

	if (!global_treat_references_as_individual)
	{
	  unsigned best_score = 0;

	  // Find best score:
	  for (std::vector<vulgar*>::size_type i=0; i<hits_to_use_for_target.size(); ++i)
	  {
	    if (hits_to_use_for_target[i] != NULL &&
		hits_to_use_for_target[i]->get_score() > best_score )
	    {
	      best_score = hits_to_use_for_target[i]->get_score();
	    }
	  }

	  // Remove all hits that do not have the best score:
	  for (std::vector<vulgar*>::size_type i=0; i < hits_to_use_for_target.size(); ++i)
	  {
	    if (hits_to_use_for_target[i] != NULL &&
		hits_to_use_for_target[i]->get_score() < best_score)
	    {
	      hits_to_use_for_target[i] = NULL;
	      vector_of_hit_stats_for_targets[i].increment_not_best_score();
	    }
	  }
	}

      }

      // Hits have been assigned to references. Now they should be aligned and
      // added to the alignment files.

      // Hier weiter:
      // Maybe determine if we have a hit at all??
      ////// XXXXXXXXXXXXXXXXXXXXXX

      CSequence_Mol* theSeq = seqs_DNA_reads_target.get_seq_by_name(target_key);
      for (std::vector<vulgar*>::size_type i=0; i<hits_to_use_for_target.size(); ++i)
      {
	vulgar     *pv         = hits_to_use_for_target[i];
	if (pv == NULL)
	  continue;

	// if pv != NULL:
	++seq_count;
	vulgar     &vul        = *pv;
	const faststring &target_key = vul.get_targetID(); // These are the reads.

	// Consider this hit:
	///   Create function for following artificial block:
        {
          CDnaString     newseqDNA;

          if (theSeq == NULL)
          {
            cerr << "ERROR: No sequence for sequence ID " << target_key
		 << " found. This can only happen if the vulgar file does not "
	            "correspond to the fasta file. Exiting.\n";
	    good_bye_and_exit(-13);
          }
          determine_alignment_string(theSeq, vul, map_query_prot_lengths[vul.get_queryID()],
				     newseqDNA, seq_count,
				     all_maps_of_gap_sites_in_queries[i],
				     all_maps_of_insertion_sites_suggested_in_reference[i]);

	  if (global_verbosity >= 5)
	    cerr << "Determining alignment for target key (read name) : " << target_key << endl;
	  
	  result_alignments_for_references[i]->add_seq_to_alignment(CSequence_Mol::dna, target_key, newseqDNA, 'N');
	  vector_of_hit_stats_for_targets[i].increment_good_hits_added_to_alignment();

	  //          seqs_DNA_result.add_seq_to_alignment(CSequence_Mol::dna, target_key, newseqDNA, 'N');
	  //          os << ">" << target_key << "\n" << newseqDNA << "\n";  
        }
      }

      //      m_it = m_exonerate_results.find(target_key);
      
      //       if (m_it == m_exonerate_results.end() )
      //       {
      // 	cerr << "INTERNAL ERROR: No sequence name found for map key (sequence ID): " << target_key << endl;
      //        good_bye_and_exit(-8);
      //       }
      
      //      vulgar vul = m_it->second;
      
      //       cout << "XX-Working on vulgar : ";
      //       vul.print(cout);
      //       cout << endl;
      ++seqs_it;
    } // END  while (vec_it !=  l_keys.end())
      // Loop over all "vulgar strings", i.e. "reads with hit" and align them.

    //os.close();

    // Final output loop:
    for (std::vector<faststring>::size_type i=0; i<seqnames_of_references.size(); ++i)
    {
      faststring &ref_name = seqnames_of_references[i];
      unsigned query_prot_length = map_query_prot_lengths[ref_name];
      
      // Write some stats about the hits:
      if (global_verbosity >= 1)
      {
	cout << "\n\n---------------------------------\n";
	cout << "Some stats for the exonerate results for reference: " << seqnames_of_references[i] << '\n';
	cout << "---------------------------------\n";
	cout.setf(ios::left);
	//	unsigned s1 = vec_of_hits_as_in_file.size; // vec_exonerate_results.size();
	cout << setw(50) << "Number of input sequences considered:" << seqs_DNA_reads_target.GetTaxaNum() << '\n';
	cout << "Number of input sequences successful aligned with exonerate (all)\n"
	     << setw(50) << "to the amino acid sequence found in vulgar file: "
	     <<  vector_of_hit_stats_for_targets[i].hits_in_vulgar_file << endl;
	//	unsigned s2 = map_exonerate_count_results.size();
	cout << "Number of successful exonerate alignments without skipped hits.\n"
	     << setw(50) << "This is the number of aligned reads: "
	     << vector_of_hit_stats_for_targets[i].good_hits_added_to_alignment
	     << endl;
	//	cout << setw(50) << "Number of sequences in result file: " << seqs_DNA_result.GetTaxaNum() << endl;
	
	cout << "-------------------------------------------\n";
	if (global_gap_frameshift_mode == 0)
	  cout << setw(50) << "# skipped reads due to gappy alignment:            "
	       << vector_of_hit_stats_for_targets[i].skipped_G << endl;
	if (global_gap_frameshift_mode == 2)
	  cout << setw(50) << "# skipped reads due to no gaps in alignment:       "
	       << vector_of_hit_stats_for_targets[i].skipped_no_G << endl;
	
	cout << setw(50) << "# skipped reads due to frameshifts in alignment:     "
	     << vector_of_hit_stats_for_targets[i].skipped_F << endl;

	cout << setw(50) << "# skipped reads due to multiple hits:                "
	     << vector_of_hit_stats_for_targets[i].skipped_double_vulgar  << endl;

	cout << setw(50) << "# skipped reads due to low rel. score:               "
	     << vector_of_hit_stats_for_targets[i].skipped_relative_score << endl;

	cout << setw(50) << "# Added to different reference due to better score:  "
	     << vector_of_hit_stats_for_targets[i].not_best_score_for_query << endl;
	
	cout << setw(50) << "# number of hits not considered:                     "
	     << vector_of_hit_stats_for_targets[i].not_considered << endl;
	
	cout << "Gap insertion sites and lengths suggested from reads for the reference:\n";
	if (all_maps_of_insertion_sites_suggested_in_reference[i].size() > 0)
	  print_Mymap(cout, all_maps_of_insertion_sites_suggested_in_reference[i]);
	else
	  cout << "None\n";

	cout << "Gap sites in query sequences:\n";
	if ( all_maps_of_gap_sites_in_queries[i].size() > 0)
	  print_Mymap(cout,  all_maps_of_gap_sites_in_queries[i]);
	else
	  cout << "None\n";
	cout << '\n';

	if (global_verbosity >= 2)
	{
	  unsigned **coverage_profile = result_alignments_for_references[i]->get_DNA_site_profile();
	  cout << "Alignment coverage for this gene:\n";
	  print_DNA_profile(cout, coverage_profile, result_alignments_for_references[i]->GetPosNum());
	}
      }

      // Export alignmentfile:
      {
	string this_output_filename = global_alignment_output_file + ref_name.c_str() + ".fas";
	ofstream os;	
	os.open(this_output_filename.c_str());
	if (os.fail() )
	{
	  cerr << "ERROR when trying to write the alignment file. "
	          "The output file could not be opened. "
	          "Please make sure the specified path exists! Path: "
	       << this_output_filename << endl;
	  exit(-23);
	}
	result_alignments_for_references[i]->ExportSequences(os, 'f', UINT_MAX, faststring(), UINT_MAX, false);
	os.close();
      }
      
      if (!global_consensus_sequence_output_filename.empty())
      {
	faststring       consensus_seq;
	unsigned         replaced_tilde_gaps;
	vector<unsigned> coverage_vec;

	result_alignments_for_references[i]->WeightedConsensusSequence_DNA(consensus_seq,
									   global_consensus_threshold,
									   global_minimum_seq_coverage_uppercase,
									   global_minimum_seq_coverage_total,
									   &coverage_vec, global_ends_width,
									   global_weight_fraction_in_ends,
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

	string outfilename = global_consensus_sequence_output_filename + ref_name.c_str() + ".fas";

	ofstream os_con(outfilename.c_str());
	if (os_con.fail() )
	{
	  cerr << "ERROR when trying to write the consensus file. The output file could not be opened. "
	          "Please make sure the specified path exists! Path: "
	       << outfilename << endl;
	  exit(-23);
	}
	os_con << ">Consensus__" << ref_name << endl;
	os_con << consensus_seq << endl;
	os_con.close();

	// Write some final coverage stats:
	if (global_verbosity >= 1)
	{
	  unsigned mi, ma;
	  double mean, sd, Q1, Q2, Q3;
	  vec_min_max(coverage_vec, mi, ma);
	  vec_mean_sd(coverage_vec, mean, sd);
	  // Be careful: coverage_vec will be sorted after calling this function:
	  vec_median_quartile_sort_method3(coverage_vec, Q1, Q2, Q3);

	  cout << "-------------------------------------------\n";
	  cout << "Alignment coverage stats:\n";
	  cout.precision(2);
	  cout.setf(ios::fixed);
	  cout << "-------------------------------------------\n";
	  cout << setw(50) << "Length of alignment: " << query_prot_length*3 << '\n';
	  cout << setw(50) << "Coverage minimum:    " << mi   << '\n';
	  cout << setw(50) << "Coverage maximum:    " << ma   << '\n';
	  cout << setw(50) << "Coverage mean:       " << mean << '\n';
	  cout << setw(50) << "Coverage median:     " << Q2   << '\n';
	  cout << "-------------------------------------------\n";
	} // END  if (global_verbosity >= 1) // Write some final coverage stats:
      } // END if (!global_consensus_sequence_output_filename.empty())
    } // for (int i=0; i<seqnames_of_references.size(); ++i) // Final output loop:
  } // END: After reading the data, this is the main analysis and alignment block.
  
  if (num_WARNINGS > 0)
  {
    if (global_verbosity >= 1)
      cout << "There have been " << num_WARNINGS << " warnings or notes. See above for details." << endl;
    else
      cout << "Due to the verbosity parameter set to 0, you missed "
      << num_WARNINGS
      << " warnings or notes. Increase the verbosity if you want to see the details and run again." << endl;
  }
  
  cout << PROGNAME " finished successfully." << endl << endl;
} // END: main(...)
