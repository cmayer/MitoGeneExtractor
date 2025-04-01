#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <ctime>
#include <iomanip>
#include "faststring3.h"
#include "CSequences3.1.h"
#include "CSequence_Mol3.1.h"
#include "Ctriple.h"
#include "CFile/CFile2_3.h"
#include "CDnaString3.h"
#include "global-types-and-parameters_MitoGeneExtractor.h"
#include "statistic_functions.h"
#include "fastq.h"
#include "Cfastq-sequences3.1.h"
#include <unordered_set>

#include "exonerate_wrapper_and_parser.hpp"



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

//bool     add_partial_codons_at_ends      = true;
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

    CSequences3 *pseqs = new CSequences3(SeqType_dna);
    fq_in_collection.add_sequences_to_CSequences_object(pseqs);
    pseqs->ExportSequences(ofp, 'f', UINT_MAX);
    delete pseqs;
  }
}


//struct stats_for_given_target
//{
//  // All numbers will be for this specific target
//  unsigned skipped_double_vulgar;
//  unsigned skipped_G;
//  unsigned skipped_F;
//  unsigned skipped_relative_score;
//  unsigned skipped_no_G;
//  unsigned not_considered;
//  unsigned not_best_score_for_query;
//  unsigned good_hits_added_to_alignment;
//  unsigned hits_in_vulgar_file;
//
//  stats_for_given_target():skipped_double_vulgar(0), skipped_G(0),
//  skipped_F(0), skipped_relative_score(0),
//  skipped_no_G(0), not_considered(0),
//  not_best_score_for_query(0),
//  good_hits_added_to_alignment(0),
//  hits_in_vulgar_file(0)
//  {}
//
//  void increment_skipped_double_vulgar()
//  {
//    ++skipped_double_vulgar;
//  }
//
//  void increment_skipped_G()
//  {
//    ++skipped_G;
//  }
//
//  void increment_skipped_F()
//  {
//    ++skipped_F;
//  }
//
//  void increment_skipped_relative_score()
//  {
//    ++skipped_relative_score;
//  }
//
//  void increment_skipped_no_G()
//  {
//    ++skipped_no_G;
//  }
//
//  void increment_not_considered()
//  {
//    ++not_considered;
//  }
//
//  void increment_not_best_score()
//  {
//    ++not_best_score_for_query;
//  }
//
//  void increment_good_hits_added_to_alignment()
//  {
//    ++good_hits_added_to_alignment;
//  }
//
//  void increment_hits_in_vulgar_file()
//  {
//    ++hits_in_vulgar_file;
//  }
//};

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

//class vulgar
//{
//public:
//  static map<faststring, short>  queryID2Index;
//
//  faststring queryID;          // Exonerate query sequences are the protein reference sequences, e.g. AA-COI consensus seq. exonerate aligns against.
//  short      queryIndex;       // Every query sequence (reference sequence) has a query number.
//                               // This avoids many lookups for index.
//  unsigned   query_start;      // 0 based numbers.
//  unsigned   query_end;        // The first position after the specified range.
//  char       query_strand;
//  faststring targetID;         // Exonerate targets are the DNA input sequences, typically reads.
//  unsigned   target_start;     // 0 based numbers.
//  unsigned   target_end;       // The first position after the specified range.
//  char       target_strand;
//  unsigned   score;
//  vector<Ctriple<char, unsigned, unsigned> > attributes;
//
//
//  // Frequent labels:
//  bool       bool_has_M;
//  bool       bool_has_G;
//  bool       bool_has_F;
//  // Rare labels:
//  bool       bool_has_5;
//  bool       bool_has_3;
//  bool       bool_has_C;
//  bool       bool_has_N;
//  bool       bool_has_I;
//  bool       bool_has_S;
//  bool       bool_has_non_M;
//
//public:
//  static short add_query_to_query2short_map_or_get_index(faststring q)
//  {
//    map<faststring, short>::iterator find_it = queryID2Index.find(q);
//    if (find_it != queryID2Index.end())
//      return find_it->second;
//
//    unsigned s = queryID2Index.size();
//    queryID2Index[q] = s;
//    return s;
//  }
//
//  short get_index_in_query2short_map(faststring q)
//  {
//    map<faststring, short>::iterator find_it = queryID2Index.find(q);
//    if (find_it != queryID2Index.end())
//      return find_it->second;
//    else
//    {
//      cerr << "Reference sequence appears in Exonerate output, but not in "
//      "reference sequences. This should not be possible.\n";
//      good_bye_and_exit(-127);
//      return 0; // Only introduced to silence warnings about no return for non-void function.
//    }
//  }
//
//private:
//  void parse(faststring &str)
//  {
//    bool_has_M = false;
//    bool_has_G = false;
//    bool_has_F = false;
//    bool_has_5 = false;
//    bool_has_3 = false;
//    bool_has_C = false;
//    bool_has_N = false;
//    bool_has_I = false;
//    bool_has_S = false;
//    bool_has_non_M = false;
//
//    attributes.clear();
//
//    vector<faststring> l(19);
//    split(l, str);
//
//    if (l.size() < 13) // 1+9+3
//    {
//
//      cerr << "ERROR: Bad vulgar string encountered with a wrong number of elements: " << str << endl;
//      good_bye_and_exit(-4);
//    }
//
//    if (l[0] != "vulgar:")
//    {
//      cerr << "ERROR: Bad vulgar string: The vulgar string is expected to start with \"vulgar:\" but found: " << str << endl;
//      good_bye_and_exit(-4);
//    }
//
//    // It is a convention in exonerate that the protein sequence is the query and the DNA sequence the target:
//    queryID       = l[1];
//    query_start   = l[2].ToUnsigned();
//    query_end     = l[3].ToUnsigned();
//    query_strand  = l[4][0];
//
//    targetID      = l[5];
//    target_start  = l[6].ToUnsigned();
//    target_end    = l[7].ToUnsigned();
//    target_strand = l[8][0];
//
//    score         = l[9].ToUnsigned();
//
//    queryIndex    = get_index_in_query2short_map(queryID);
//
//    int m = l.size()-10;
//    if (m%3 != 0)
//    {
//      cerr << "ERROR: Bad vulgar string encountered. The number of fields must "
//      "by 10+3*n where n is an integer. But the number is "
//      << m << "." << endl;
//      cerr << "       Vulgar string: " << str << endl;
//      good_bye_and_exit(-5);
//    }
//
//    unsigned k = 10;
//    Ctriple<char, unsigned, unsigned> t('\0', 0, 0);
//    while (k < l.size())
//    {
//      char c        = l[k][0]; ++k;
//      unsigned num1 = l[k].ToUnsigned(); ++k;
//      unsigned num2 = l[k].ToUnsigned(); ++k;
//
//      if (c == 'M')
//        bool_has_M = true;
//      else if (c == 'G')
//      {
//        bool_has_G     = true;
//        bool_has_non_M = true;
//      }
//      else if (c == 'F')
//      {
//        bool_has_F = true;
//        bool_has_non_M = true;
//      }
//      else if (c == '5')
//      {
//        bool_has_5 = true;
//        bool_has_non_M = true;
//      }
//      else if (c == '3')
//      {
//        bool_has_3 = true;
//        bool_has_non_M = true;
//      }
//      else if (c == 'C')
//      {
//        bool_has_C = true;
//        bool_has_non_M = true;
//      }
//      else if (c == 'N')
//      {
//        bool_has_N = true;
//        bool_has_non_M = true;
//      }
//      else if (c == 'I')
//      {
//        bool_has_I = true;
//        bool_has_non_M = true;
//      }
//      else if (c == 'S')
//      {
//        bool_has_S = true;
//        bool_has_non_M = true;
//      }
//
//      t = Ctriple<char, unsigned, unsigned>(c, num1, num2);
//      attributes.push_back(t);
//    }
//  }
//
//  /* Does not seem to make sense to combine F and G. There is a simple and consistent interpretation of the two independent of the other.
//   void combine_FG_attributes()
//   {
//   unsigned i, N=attributes.size();
//
//   for (i=1; i<N; ++i)
//   {
//   if (attributes[i-1].first() == 'G' && attributes[i].first() == 'F')
//   {
//   if (attributes[i-1].second()>0 && attributes[i-1].third() == 0 && (attributes[i].second()==0 && attributes[i-1].third() > 0) )
//   {
//   Ctriple<char, unsigned, unsigned> t('\0', 0, 0);
//   t = Ctriple<char, unsigned, unsigned>('f', attributes[i-1].second(), attributes[i].third());
//   attributes[i-1] = t;
//   attributes.erase(attributes.begin()+i);
//   // Do we have to adapt i after removing an element. In this case no, since we will not combine anything with the newly inserted feature 'f'.
//   }
//   else
//   {
//   cout << "WARNING: Unexpected combination of lengths in successive G and F attribute in vulgar line:\n";
//   print(cout);
//   }
//   }
//   else if (attributes[i-1].first() == 'F' && attributes[i].first() == 'G')
//   {
//   if (attributes[i-1].second()==0 && attributes[i-1].third() >0 && (attributes[i].second()>0 && attributes[i].third() == 0) )
//   {
//   Ctriple<char, unsigned, unsigned> t('\0', 0, 0);
//   t = Ctriple<char, unsigned, unsigned>('f', attributes[i].second(), attributes[i-1].third());
//   attributes[i-1] = t;
//   attributes.erase(attributes.begin()+i);
//   // Do we have to adapt i after removing an element. In this case no, since we will not combine anything with the newly inserted feature 'f'.
//   }
//   else
//   {
//   cout << "WARNING: Unexpected combination of lengths in successive F and G attributes in vulgar line:\n";
//   print(cout);
//   }
//   }
//   }
//   }
//   */
//
//public:
//  void print(ostream &os) const
//  {
//    os << "vulgar: "
//    << queryID       << " "
//    //    << queryIndex    << " "
//    << query_start   << " "
//    << query_end     << " "
//    << query_strand  << " "
//    << targetID      << " "
//    << target_start  << " "
//    << target_end    << " "
//    << target_strand << " "
//    << score << " ";
//    for (unsigned i=0; i<attributes.size(); ++i)
//    {
//      os << attributes[i].first()  << " ";
//      os << attributes[i].second() << " ";
//      os << attributes[i].third() << " ";
//    }
//  }
//
//  vulgar(faststring &str)
//  {
//    parse(str);
//    //    combine_FG_attributes();
//  }
//
//  const faststring & get_queryID() const
//  {
//    return queryID;
//  }
//
//  short get_queryIndex() const
//  {
//    return queryIndex;
//  }
//
//  const faststring & get_targetID() const
//  {
//    return targetID;
//  }
//
//  faststring get_hitkey() const
//  {
//    return queryID + "##" + targetID;
//  }
//
//  bool has_F() const
//  {
//    return bool_has_F;
//  }
//
//  bool has_G() const
//  {
//    return bool_has_G;
//  }
//
//  bool has_M() const
//  {
//    return bool_has_M;
//  }
//
//  bool has_non_M() const
//  {
//    return bool_has_non_M;
//  }
//
//  size_t num_attributes() const
//  {
//    return attributes.size();
//  }
//
//  unsigned get_query_start() const
//  {
//    return query_start;
//  }
//
//  unsigned get_target_start() const
//  {
//    return target_start;
//  }
//
//  unsigned get_target_end() const
//  {
//    return target_end;
//  }
//
//  bool is_revcomp() const
//  {
//    return (target_strand == '-');
//  }
//
//  // Divides the exonerate score by the length of the hit and returns the result.
//  double relative_score() const
//  {
//    double dist;
//    if (target_end > target_start)
//      dist = target_end - target_start;
//    else
//      dist = target_start - target_end;
//    if (dist == 0)
//      return 0;
//    return score/dist;
//  }
//
//  unsigned get_score()
//  {
//    return score;
//  }
//};


// Global variable definition for static member variable of vulgar class:
map<faststring, short>  vulgar::queryID2Index;


/* Currently not used:
inline size_t size_keyset_multimap(multimap<faststring, vulgar *> &mm)
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
*/

inline size_t number_of_entries_targetID(multimap<faststring, vulgar *> &mm, faststring targetID)
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



// For given vulgar data and corresponding input sequence, this function
// determines the aligned sequence newseqDNA.
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
    cerr << "DEBUG(100): determine_alignment_string was called with the following parameters:\n";
    cerr << "DEBUG(100): query_prot_length:  " << query_prot_length << endl;
    cerr << "DEBUG(100): Reference (query):  " << vul.get_queryID() << endl;
    cerr << "DEBUG(100): Reference (target): " << vul.get_targetID() << endl;
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
        cout << "DEBUG(1000): Sequence:    " << seq_count << endl;
        cout << "DEBUG(1000): gaps_before: " << gaps_before << endl;
      }

      if (vul.is_revcomp())
      {
        unsigned add_length_beginning = 0;
        unsigned new_start_in_seqDNA;

        if (global_verbosity >= 1000)
          cout << "DEBUG(1000): Number add_length_beginning: " << add_length_beginning << endl;

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

        if (/* DISABLES CODE */ (0))
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
        cout << "DEBUG(1000): Add partial codon front \"" << add_seq_before << "\" might have length 0." << endl;

      if (global_verbosity >= 50)
        cout << "DEBUG(50): Add front length: " << add_seq_before.length() << " " << theSeq->getName() << endl;

    }  // END if (gaps_before > 0)
  } // END if (global_num_bp_beyond_exonerate_alignment_if_at_end > 0)

  //   cout << "Add gaps before: " << gaps_before << endl;
  gaps_before -= add_seq_before.length();
  newseqDNA = faststring((faststring::size_type)gaps_before, '-');

  // Debug code:
  add_seq_before.ToLower();

  if (add_seq_before.size() > 0)
    newseqDNA += add_seq_before;

  unsigned DNA_pos   = start_in_seqDNA;
  unsigned amino_pos = vul.query_start; // Only needed to store insertions with respect to reference.
  unsigned count     = (unsigned)(gaps_before + add_seq_before.length());

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
          newseqDNA += faststring((faststring::size_type)gap_chars, '~');
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
            cout << "WARNING: Rare gap case: G 0 x. Additional bases are skipped for sequence ID " << vul.get_targetID() << endl;
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
    cout << "DEBUG(40): " << theSeq->getName() << " gaps_after: " << gaps_after << endl;

  if (global_num_bp_beyond_exonerate_alignment_if_at_start_or_end > 0)
  {
    if (gaps_after > 0)
    {
      if (global_verbosity >= 1000)
      {
        cout << endl;
        cout << "DEBUG(1000): Sequence:   " << seq_count << endl;
        cout << "DEBUG(1000): gaps_after: " << gaps_after << endl;
      }
      if (vul.is_revcomp())
      {
        // Default values for start and length of extracted region:
        unsigned add_length_end = length_this_nuc_seq - end_in_seqDNA;
        //      vul.get_target_end();
        // Where in seqDNA do we start to extract the extra bases
        unsigned new_start_in_seqDNA = end_in_seqDNA;

        if (global_verbosity >= 1000)
          cout << "DEBUG(1000): Number potential_add_at_end_of_inserted_revcomp_seq: " << add_length_end << endl;

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
          cout << "DEBUG(1000): Number add_length_end: " << add_length_end << endl;

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
        cout << "DEBUG(1000): Add at back \"" << add_seq_after << "\" might have length 0." << endl;

      if (global_verbosity >= 50)
        cout << "DEBUG(50): Add back length: " << add_seq_after.length() << " " << theSeq->getName() << endl;

      gaps_after -= add_seq_after.length();
    } // END  if (gaps_after > 0)
  } // END if (global_num_bp_beyond_exonerate_alignment_if_at_start_or_end > 0)

  // Debug code:
  add_seq_after.ToLower();

  //   cout << "Add gaps after: " << gaps_after << endl;
  if (add_seq_after.length() > 0)
    newseqDNA += add_seq_after;

  newseqDNA += faststring((faststring::size_type)gaps_after, '-');

  //   cout << "gaps_before count : " << gaps_before << " " << count << endl;
}


//int run_exonerate(const faststring &exonerate_binary, const faststring &protein_reference, const faststring &fasta_reads, const faststring &vulgar_base_file, char genetic_code_number, int frameshift_penalty, int max_trials, char mode=0)
//{
//  faststring code_string;
//
//  if (genetic_code_number == 1) // The Standard Code (transl_table=1)
//    code_string = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
//  else if (genetic_code_number == 2) // The Vertebrate Mitochondrial Code
//    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
//  else if (genetic_code_number == 3)  // The Yeast Mitochondrial Code
//    code_string = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
//  else if (genetic_code_number == 4)  // The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
//    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
//  else if (genetic_code_number == 5)  // The Invertebrate Mitochondrial Code
//    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
//  else if (genetic_code_number == 6)  // The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
//    code_string = "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
//
//
//  else if (genetic_code_number == 9)  // The Echinoderm and Flatworm Mitochondrial Code
//    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
//  else if (genetic_code_number == 10) // The Euplotid Nuclear Code (transl_table=10)
//    code_string = "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
//  else if (genetic_code_number == 11) // The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
//    code_string = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
//  else if (genetic_code_number == 12) // The Alternative Yeast Nuclear Code (transl_table=12)
//    code_string = "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
//  else if (genetic_code_number == 13) // The Ascidian Mitochondrial Code
//    code_string = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG";
//  else if (genetic_code_number == 14) // The Alternative Flatworm Mitochondrial Code
//    code_string = "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
//  else if (genetic_code_number == 16) // Chlorophycean Mitochondrial Code
//    code_string = "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
//  else
//  {
//    cerr << "ERROR: Unknown genetic code number passed to function run_exonerate. Code number: " << genetic_code_number << endl;
//    good_bye_and_exit(-20);
//  }
//
//  //  system "exonerate --geneticcode FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG --frameshift -9 --query ${dir2}/${gene}.fasta
//  // -Q protein --target $dir1/${name}.fas -T dna --model protein2dna --showalignment 0 --showvulgar 1 > $dir3/${gene}_${name}_vulgar.txt 2> $dir3/${gene}_${name}_vulgar.log";
//
//  int system_return_value;
//  faststring score_threshold;
//  if (global_exonerate_score_threshold != UINT_MAX)
//    score_threshold = " -s " + faststring(global_exonerate_score_threshold) + " ";
//  faststring command =   exonerate_binary + " --geneticcode " + code_string + score_threshold
//  + " --frameshift " + faststring(frameshift_penalty)
//  + " --query " + protein_reference + " -Q protein --target " + fasta_reads
//  + " -T dna --model protein2dna --showalignment 0 --showvulgar 1 1> "
//  + vulgar_base_file + " 2> " + vulgar_base_file + ".log";
//
//  // Special debug mode to examine command string:
//  if (mode == 'd')
//  {
//    cout << "DEBUG: Exonerate command: " << command << endl;
//    return 0;
//  }
//
//  int i;
//  for (i=1; i <= max_trials; ++i)
//  {
//    system_return_value = system(command.c_str());
//    if (system_return_value == 0)
//      return i;
//    // Else
//    // Clean up the unsuccessful run:
//    //    sleep(1); // Sleep for 1 second to allow system to close all files.
//    faststring remove_file = vulgar_base_file;
//
//    remove(remove_file.c_str());
//    faststring log_file_name = vulgar_base_file + ".log";
//    remove(log_file_name.c_str());
//  }
//  return max_trials +1;
//}


//************************************************
// seqnames_of_references
// map_of_vulgar_hits_targets_as_keys
// vector_of_hit_stats_for_targets

// Remember: Reads/sequences are targets.
//           Reference sequences are queries.

int align_sequences_using_vulgar(const vector<faststring>  &seqnames_of_references,
                                 const multimap<faststring, vulgar *> &map_of_vulgar_hits_targets_as_keys,
                                 vector<stats_for_given_target> &vector_of_hit_stats_for_query_reference,
                                 CSequences3 &seqs_DNA_reads_target,
                                 map<faststring, size_t> &map_query_prot_lengths,
                                 vector<CSequences3*> &result_alignments_for_references,
                                 vector<map<pair<unsigned, unsigned>, unsigned> > &all_maps_of_insertion_sites_suggested_in_reference,
                                 vector<map<pair<unsigned, unsigned>, unsigned> > &all_maps_of_gap_sites_in_queries
                                 )
{
  multimap<faststring, vulgar *>::const_iterator seqs_it, equal_range_start, equal_range_end, tmp_it;
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
      cout << "PROGRESS: Creating alignments for (read) sequence ID: " << target_key << endl;
    }

    // Determine equal range. We could call equal range, but we already
    // have the beginning, so it is more efficient to only look for the end.
    tmp_it = equal_range_start = seqs_it;
    unsigned count_hits_this_read_target_before_final_filtering = 0;

    // A DNA input sequence might have hits with multiple references.
    // Furthermore, the hits have to be added to the right alignment.
    // This will be taken care of in the following:
    // Exonerate targets: Input DNA sequences
    // Econerate queries: Protein references.

    map<faststring, unsigned> count_hits_for_different_queries_of_this_target; // Count hits associated to the reference sequences for this read.
    vector<vulgar *>          vec_of_vulgar_hits_for_this_read_target_key;
    vector<vulgar *>          hits_assigned_to_references_for_this_target(seqnames_of_references.size(), 0); // Stores for all target hits the references the hits are assigned to.

    // Loop over all hits for this target.
    while (tmp_it != map_of_vulgar_hits_targets_as_keys.end() && tmp_it->first == target_key)
    {
      vec_of_vulgar_hits_for_this_read_target_key.push_back(tmp_it->second);
      ++count_hits_this_read_target_before_final_filtering;
      add_or_count(count_hits_for_different_queries_of_this_target, (*(tmp_it->second)).get_queryID() ); // Query is reference
      ++tmp_it;
    }
    equal_range_end = tmp_it;

    // This is typically not true if we have mutliple queries/references, or if queries/references appear multiple times.
    if (count_hits_this_read_target_before_final_filtering == 1)
    {
      if (vec_of_vulgar_hits_for_this_read_target_key.size() != 1)
      {
        cerr << "vec_of_hits_for_this_target_key has unexpected size. "
        "Should be 1. Has : "
        << vec_of_vulgar_hits_for_this_read_target_key.size() << endl;
      }
      vulgar *p = vec_of_vulgar_hits_for_this_read_target_key[0];
      hits_assigned_to_references_for_this_target[p->get_queryIndex()] = p;   // The query index is the index of the AA reference sequence.
    }
    else // We have multiple hits for target (i.e. read) sequences:
    {
      // We have more than one hit.
      // This situation is expected if we have competing references for the same gene
      // or long reads and multiple genes.

      for (std::vector<vulgar*>::size_type hit_index_this_target=0; hit_index_this_target<vec_of_vulgar_hits_for_this_read_target_key.size(); ++hit_index_this_target)
      {
        vulgar *p                = vec_of_vulgar_hits_for_this_read_target_key[hit_index_this_target];
        vulgar &vul              = *p;
        short ref_index_this_hit = vul.get_queryIndex(); // Get index of reference sequence

        // Address NULL       indicates: No hit found so far.
        // Address (void *)-1 indicates: Double hit found for this reference already.

        // Is this the first hit against this reference we have found so far for this target?
        if (hits_assigned_to_references_for_this_target[ref_index_this_hit] == NULL )
        {
          hits_assigned_to_references_for_this_target[ref_index_this_hit] = p;
        }
        else // Found more than one hit against this target. This is to be handled.
        {
          if (global_include_double_hits && hits_assigned_to_references_for_this_target[ref_index_this_hit] !=  (void *)-1 )
          {
            // We include double hits, but only the best one of multiple hits:
            // Is this score better than the one of the previous hit?
            if (vul.get_score() > hits_assigned_to_references_for_this_target[ref_index_this_hit]->get_score() )
              hits_assigned_to_references_for_this_target[ref_index_this_hit] = p;
            // We used only the better of two hits, so we skipped the other one.
            vector_of_hit_stats_for_query_reference[ref_index_this_hit].increment_skipped_double_vulgar();
          }
          else // We have a double hit, but do not want to include them.
               // So we block this reference for this target.
          {
            // Best hits for this read with respect to the different references are stored in hits_to_use_for_target
            hits_assigned_to_references_for_this_target[ref_index_this_hit] = (vulgar *)-1; // We only add a hit if pointer is 0.
                                                              //        ++skipped_double_vulgar;                        // -> So we block this reference.
                                                              // We skip two hits, one that has already been entered and the one we just tried to enter.
            vector_of_hit_stats_for_query_reference[ref_index_this_hit].increment_skipped_double_vulgar();
            vector_of_hit_stats_for_query_reference[ref_index_this_hit].increment_skipped_double_vulgar();  // This might be incremented too often for triple hits.
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

      // References with pointer (void *)-1) are blocked due to double hits.
      for (std::vector<vulgar*>::size_type ref_loop_index=0; ref_loop_index<hits_assigned_to_references_for_this_target.size(); ++ref_loop_index)
      {
        if (hits_assigned_to_references_for_this_target[ref_loop_index] == (vulgar *)-1 ) // This was only used to indicate that the reference is blocked.
          hits_assigned_to_references_for_this_target[ref_loop_index] = NULL;             // Now we indicate that no hit is assigned to this reference.
      }

      // In this case, we only keep the best hit, so references are compeeting for the hits.
      if (!global_treat_references_as_individual)
      {
        unsigned best_score = 0;

        // Find best score:
        for (std::vector<vulgar*>::size_type ref_loop_index=0; ref_loop_index<hits_assigned_to_references_for_this_target.size(); ++ref_loop_index)
        {
          if (hits_assigned_to_references_for_this_target[ref_loop_index] != NULL &&
              hits_assigned_to_references_for_this_target[ref_loop_index]->get_score() > best_score )
          {
            best_score = hits_assigned_to_references_for_this_target[ref_loop_index]->get_score();
          }
        }

        // Remove all hits that do not have the best score:
        for (std::vector<vulgar*>::size_type ref_loop_index=0; ref_loop_index < hits_assigned_to_references_for_this_target.size(); ++ref_loop_index)
        {
          if (hits_assigned_to_references_for_this_target[ref_loop_index] != NULL &&
              hits_assigned_to_references_for_this_target[ref_loop_index]->get_score() < best_score)
          {
            hits_assigned_to_references_for_this_target[ref_loop_index] = NULL;
            vector_of_hit_stats_for_query_reference[ref_loop_index].increment_not_best_score();
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
    for (std::vector<vulgar*>::size_type ref_loop_index=0; ref_loop_index<hits_assigned_to_references_for_this_target.size(); ++ref_loop_index)
    {
//      vulgar     *pv         = hits_to_use_for_target[i];
      vulgar *p_vul       = hits_assigned_to_references_for_this_target[ref_loop_index];

      if (p_vul == NULL)
        continue;

      // if p_vul != NULL:
      ++seq_count;
      vulgar     &vul        = *p_vul;
      short ref_index = vul.get_queryIndex(); // Get index of reference sequence

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
                                   all_maps_of_gap_sites_in_queries[ref_loop_index],
                                   all_maps_of_insertion_sites_suggested_in_reference[ref_loop_index]);

        if (global_verbosity >= 5)
          cerr << "PROGRESS: Determining alignment for DNA seqID (target) and ref_index: " << target_key << ", " << ref_index << endl;

        result_alignments_for_references[ref_index]->add_seq_to_dataset(SeqType_dna, target_key, newseqDNA, 'N');
        vector_of_hit_stats_for_query_reference[ref_index].increment_good_hits_added_to_alignment();
        // Add DNA seq to seqs_included_in_alignments_for_references:
        //          seqs_included_in_alignments_for_references[i].insert(theSeq);

        //          seqs_DNA_result.add_seq_to_alignment(CSequence_Mol::dna, target_key, newseqDNA, 'N');
        //          os << ">" << target_key << "\n" << newseqDNA << "\n";
      }
    }

    //      m_it = m_exonerate_results.find(target_key);

    //       if (m_it == m_exonerate_results.end() )
    //       {
    //   cerr << "INTERNAL ERROR: No sequence name found for map key (sequence ID): " << target_key << endl;
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
  for (std::vector<faststring>::size_type ref_ind=0; ref_ind<seqnames_of_references.size(); ++ref_ind)
  {
    faststring ref_name = seqnames_of_references[ref_ind];
    size_t query_prot_length = map_query_prot_lengths[ref_name];

    // Write some stats about the hits:
    if (global_verbosity >= 1)
    {
      cout << "\n\nINFO: ---------------------------------\n";
      cout << "INFO: Some stats for the exonerate results for reference: " << seqnames_of_references[ref_ind] << '\n';
      cout << "INFO: ---------------------------------\n";
      cout.setf(ios::left);
      //  unsigned s1 = vec_of_hits_as_in_file.size; // vec_exonerate_results.size();
      cout << setw(70) << "INFO: Number of input sequences:" << seqs_DNA_reads_target.GetTaxaNum() << '\n';
      cout << "INFO: Number of input sequences successful aligned with exonerate\n";
      cout << left << setw(70) << "INFO: to the amino acid sequence, see vulgar file: "
      <<  vector_of_hit_stats_for_query_reference[ref_ind].hits_in_vulgar_file << endl;
      //  unsigned s2 = map_exonerate_count_results.size();
      cout << "INFO: Number of successfully aligned sequences without skipped hits.\n"
      << setw(70) << "INFO: This is the number of aligned reads: "
      << vector_of_hit_stats_for_query_reference[ref_ind].good_hits_added_to_alignment
      << endl;
      //  cout << setw(70) << "Number of sequences in result file: " << seqs_DNA_result.GetTaxaNum() << endl;

      cout << "INFO: -------------------------------------------\n";
      if (global_gap_frameshift_mode == 0)
        cout << "INFO: Skipped reads due to gappy alignment:            "
        << vector_of_hit_stats_for_query_reference[ref_ind].skipped_G << endl;
      if (global_gap_frameshift_mode == 2)
        cout << "INFO: Skipped reads due to no gaps in alignment:       "
        << vector_of_hit_stats_for_query_reference[ref_ind].skipped_no_G << endl;

      cout << "INFO: Skipped reads due to frameshifts in alignment:     "
      << vector_of_hit_stats_for_query_reference[ref_ind].skipped_F << endl;

      cout << "INFO: Skipped reads due to multiple hits:                "
      << vector_of_hit_stats_for_query_reference[ref_ind].skipped_double_vulgar  << endl;

      cout << "INFO: Skipped reads due to low rel. score:               "
      << vector_of_hit_stats_for_query_reference[ref_ind].skipped_relative_score << endl;

      cout << "INFO: Added to different reference due to better score:  "
      << vector_of_hit_stats_for_query_reference[ref_ind].not_best_score_for_query << endl;

      cout << "INFO: Total number of hits not considered:               "
      << vector_of_hit_stats_for_query_reference[ref_ind].not_considered << endl;

      cout << "INFO: -------------------------------------------\n";
      cout << "INFO: Gap insertion sites and lengths suggested from reads for the reference:\n";
      cout << "INFO: -------------------------------------------\n";
      if (all_maps_of_insertion_sites_suggested_in_reference[ref_ind].size() > 0)
        print_Mymap(cout, all_maps_of_insertion_sites_suggested_in_reference[ref_ind]);
      else
        cout << "None\n";

      cout << "INFO: -------------------------------------------\n";
      cout << "INFO: Gap sites in query sequences:\n";
      cout << "INFO: -------------------------------------------\n";
      if ( all_maps_of_gap_sites_in_queries[ref_ind].size() > 0)
        print_Mymap(cout,  all_maps_of_gap_sites_in_queries[ref_ind]);
      else
        cout << "None\n";
      cout << '\n';

      if (global_verbosity >= 2)
      {
        unsigned **coverage_profile = result_alignments_for_references[ref_ind]->get_DNA_site_profile();
        cout << "INFO: Alignment coverage for this gene:\n";
        print_DNA_profile(cout, coverage_profile, result_alignments_for_references[ref_ind]->GetPosNum());
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
      result_alignments_for_references[ref_ind]->ExportSequences(os, 'f', UINT_MAX, faststring(), UINT_MAX, false);
      os.close();
    }

    if (!global_consensus_sequence_output_filename.empty())
    {
      faststring       consensus_seq;
      size_t           replaced_tilde_gaps;
      vector<unsigned> coverage_vec;

      result_alignments_for_references[ref_ind]->WeightedConsensusSequence_DNA(consensus_seq,
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

      if (global_verbosity >= 100)
      {
         cerr << "DEBUG(100): Coverage values obtained while determining the consensus sequence:\n";
         for (unsigned i=0; i<coverage_vec.size(); ++i)
         cerr << i << "\t" << coverage_vec[i] << endl;
         cerr << "=====================\n";
      }


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
        unsigned mi=0, ma=0;
        double mean=0, sd=0, Q1=0, Q2=0, Q3=0;
        vec_min_max(coverage_vec, mi, ma);
        vec_mean_sd(coverage_vec, mean, sd);
        // Be careful: coverage_vec will be sorted after calling this function:
        vec_median_quartile_sort_method3(coverage_vec, Q1, Q2, Q3);

        cout << "INFO: -------------------------------------------\n";
        cout << "INFO: Alignment coverage stats:\n";
        cout.precision(2);
        cout.setf(ios::fixed);
        cout << "INFO: -------------------------------------------\n";
        cout << setw(50) << "INFO: Length of alignment: " << query_prot_length*3 << '\n';
        cout << setw(50) << "INFO: Coverage minimum:    " << mi   << '\n';
        cout << setw(50) << "INFO: Coverage maximum:    " << ma   << '\n';
        cout << setw(50) << "INFO: Coverage mean:       " << mean << '\n';
        cout << setw(50) << "INFO: Coverage median:     " << Q2   << '\n';
        cout << "INFO: -------------------------------------------\n";
      } // END  if (global_verbosity >= 1) // Write some final coverage stats:
    } // END if (!global_consensus_sequence_output_filename.empty())
  } // for (int i=0; i<seqnames_of_references.size(); ++i) // Final output loop:


  return 0;
} // END: align_sequences_using_vulgar()







int main(int argc, char **argv)
{
  cout << "Welcome to the " PROGNAME " program, version " VERSION << endl;

  // In exonerate the protein sequence is by convention the query sequence.
  // This makes things a bit confusing here:

  read_and_init_parameters(argc, argv);
  if (global_verbosity >= 1)
  {
    print_parameters(cout, "INFO: ");
    cout << "INFO: Progress:\n";
  }

  bool combined_input_sequence_file_created = false;
  CSequences3 seqs_DNA_reads_target(SeqType_dna);
  CSequences3 seqs_prot_query(SeqType_protein);

  //  CSequences2 seqs_DNA_result(CSequence_Mol::dna);
  // unsigned skipped_double_vulgar=0, skipped_G=0, skipped_F=0, skipped_relative_score=0, skipped_no_G=0;
  //  vector<unordered_set<CSequence_Mol*> > seqs_included_in_alignments_for_references;

  string global_input_dna_fasta_file;

  //****************************
  //*** BEGIN: concat if multiple files
  //****************************

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
    if (global_verbosity >= 1000)
      cerr << "DEBUG(1000): Temporary filename position 0: " << tmpstring << endl;

    char  *tmpstr = new char [global_tmp_directory.length() + 40];
    strcpy (tmpstr, tmpstring.c_str());
    if (global_verbosity >= 1000)
      cerr << "DEBUG(1000): Temporary filename position 1: " << tmpstr << endl;
    int   fd = mkstemp(tmpstr);

    if (global_verbosity >= 1)
      cerr << "PROGRESS: Temporary fasta file has been created: " << tmpstr << endl;

    FILE  *ofp = fdopen(fd, "w");
    combined_input_sequence_file_created = true;
    if (fd == -1 || ofp == 0)
    {
      cerr << "ERROR: Could not open the temporary fasta file that should be "
      "passed to the exonerate program for writing." << endl;
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

  //****************************
  //*** END: concat if multiple files
  //****************************


  //**************************
  // BEGIN: Read reference sequences
  //**************************

  CFile query_prot_file;

  if (global_verbosity >= 3)
    cout << "PROGRESS: Reading fasta file with protein reference sequence: " << global_input_prot_reference_sequence << endl;

  query_prot_file.ffopen(global_input_prot_reference_sequence.c_str());
  if (query_prot_file.fail())
  {
    cerr << "ERROR: Could not open specified file: " << global_input_prot_reference_sequence << endl;
    good_bye_and_exit(-1);
  }
  // cerr << "Line:   " <<  file.line() << endl;
  // cerr << "Status: " <<  (int)file.status() << endl;

  CSequence_Mol::processing_flag pflag(CSequence_Mol::convert_toupper);
  int error_prot_read = seqs_prot_query.read_from_Fasta_File(query_prot_file, pflag, 0, -1, false);
  query_prot_file.ffclose();

  if (error_prot_read < 0)
  {
    cerr << "ERROR: Exiting: An error occurred while reading the protein sequence file. "
    "Double check the input type or the file you specified."
    << endl;
    good_bye_and_exit(-3);
  }

  if (!seqs_prot_query.are_short_names_unique())
  {
    cerr << "ERROR: The sequence in the protein reference file are not unique, when trimmed after the first space."
    "Please rename the reference sequences to ensure they are unique before the first space. "
    "Often it is sufficient to rename sequences by replacing spaces with e.g. underscores.\n";
    good_bye_and_exit(-3);
  }

  //**************************
  // More reference sequence data strutures:
  //**************************

  vector<faststring>   seqnames_of_references;
  seqs_prot_query.get_short_sequence_names(seqnames_of_references);

  if (global_verbosity >= 100)
    cout << "DEBUG(100): Lenth of seqnames_of_references vector: " << seqnames_of_references.size() << endl;

  map<faststring, size_t> map_query_prot_lengths;
  for (std::vector<faststring>::size_type i=0; i < seqnames_of_references.size(); ++i)
  {
    map_query_prot_lengths[seqnames_of_references[i]] = seqs_prot_query.get_seq_by_index(i)->length();

    // Add names to static member vulgar::add_query_to_query2short_map_or_get_index
    vulgar::add_query_to_query2short_map_or_get_index(seqnames_of_references[i]);
  }

  if (global_verbosity >= 100)
  {
    cerr << "DEBUG(100): Sequence names found in reference file:\n";
    for (std::vector<faststring>::size_type i=0; i<seqnames_of_references.size(); ++i)
      cerr << seqnames_of_references[i] << endl;
  }

  //**************************
  // END: Read reference sequences
  //**************************



  //**************************
  // Vulgar data strutures:
  //**************************
  vector<vulgar *>               vec_of_hits_as_in_file;
  multimap<faststring, vulgar *> map_of_vulgar_hits_targets_as_keys;

  //**************************
  // Alignment data strutures:
  //**************************
  vector<CSequences3*> result_alignments_for_references;
  for (std::vector<faststring>::size_type ref_ind=0; ref_ind < seqnames_of_references.size(); ++ref_ind)
  {
    result_alignments_for_references.push_back(new CSequences3(SeqType_dna) );
  }

  //**************************
  // Stats data strutures:
  //**************************
  vector<map<pair<unsigned, unsigned>, unsigned> > all_maps_of_insertion_sites_suggested_in_reference(seqnames_of_references.size());
  vector<map<pair<unsigned, unsigned>, unsigned> > all_maps_of_gap_sites_in_queries(seqnames_of_references.size());
  vector<stats_for_given_target>                   vector_of_hit_stats_for_query_reference(seqnames_of_references.size() );

// Loop over all input files:
// global_input_dna_fasta_filenames && global_input_dna_fastq_filenames


  //**************************
  // Read query sequences
  //**************************
  CFile query_file;

  if (global_verbosity >= 3)
  {
    cout << "PROGRESS: Reading input DNA sequences from fasta file: " << global_input_dna_fasta_file << endl;
  }

  query_file.ffopen(global_input_dna_fasta_file.c_str());
  if (query_file.fail())
  {
    cerr << "ERROR: Could not open specified file: " << global_input_dna_fasta_file << endl;
    good_bye_and_exit(-1);
  }
  // cerr << "Line:   " <<  file.line() << endl;
  // cerr << "Status: " <<  (int)file.status() << endl;


  int error = seqs_DNA_reads_target.read_from_Fasta_File(query_file, pflag, 0, -1, false);
  query_file.ffclose();

  if (error < 0)
  {
    cerr << "ERROR: Exiting: An error occurred while reading the input file: " << global_input_dna_fasta_file
         << ". Double check that the input file exists and is in fasta format. "
            "Note that fastq files have to be specified with the -q option not with the -d option.\n";
    good_bye_and_exit(-3);
  }

  // TODO:
  // Instead of printing this error message and exiting, we could fix this
  // problem at this point, at least if we compute the vulgar file ourselves.
  // But this should be an extremely  rare case. NGS reads should always be
  // unique. I do not see the scenario in which this would occur.

  if (!seqs_DNA_reads_target.are_short_names_unique())
  {
    cerr << "ERROR: The sequence names in the DNA input file (typically the read sequences) are not unique when "
    "trimming them at the first space. This trimming is done by many programs, including exonerate. Please verify your "
    "input data and correct this. Often it is sufficient to rename sequences by replacing spaces with e.g. underscores.\n";
    good_bye_and_exit(-3);
  }

  if (global_verbosity >= 3)
  {
    cout << "PROGRESS: Finished reading input DNA sequences. Found this number of input DNA sequences: "
         << seqs_DNA_reads_target.GetTaxaNum()  << endl;
  }

  // Print memory usage details:
  if (global_verbosity > 2000)
  {
    unsigned long mem1, mem2, mem3;
    float mean2;
    seqs_DNA_reads_target.memory_usage(mem1, mem2, mem3, mean2);

    {
      cout << "DEBUG(2000): Estimated memory usage of sequences in memory: " << endl
      << "class CSequences2:           " << mem1 << endl
      << "seqData vector content:      " << mem2 << endl
      << "seqData vector content: (GB) " << (float)mem2/1024/1024/1024 << endl
      << "sn_map content:              " << mem3 << endl
      << "mean per seq object:         " << mean2 << endl;
    }
  }


  //  unsigned query_prot_length = seqs_prot_query.GetPosNum();


  if (global_verbosity >= 3)
  {
    cout << "PROGRESS: Finished reading fasta file with protein reference sequence." << endl;
  }

  //**************************
  // TODO: (Not implemented yet. No plans to implement this in the near future. User has to deal with bad input files.) Export query if sequence names had to be modified.
  //**************************


  exonerate_parameters ep(global_exonerate_binary.c_str(),
                          global_input_prot_reference_sequence.c_str(),
                          global_input_dna_fasta_file.c_str(),
                          global_vulgar_directory.c_str(),
                          global_genetic_code_number,
                          global_frameshift_penalty,
                          10, // Maximum number of trials
                          global_run_mode);

  int run_result = run_and_parse_exonerate_for_single_input_file(ep, vec_of_hits_as_in_file, map_of_vulgar_hits_targets_as_keys, vector_of_hit_stats_for_query_reference);

  //  cout << "RUN RESULT: " << run_result << endl;

  //*********




  // Remove exonerate input file if this was created in this program

/*
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
*/

  // Remove the vulgar file, if the user has not specified a file name to save the file permanently.
/*
  if (global_vulgar_file == "tmp-vulgar.txt")
    remove("tmp-vulgar.txt");
*/

  // Debug output:
  if (/* DISABLES CODE */ (0))
  {
    std::cout << "Finished reading exonerate results. Please hit any key to continue.\n" << std::endl;
    std::cin.get();
  }

  if (global_verbosity >= 3)
  {
    size_t s1 = vec_of_hits_as_in_file.size();
    cout << "PROGRESS: Number of exonerate hits (including multiple hits) to the amino "
    "acid sequence\nfound in input file:, without skipped hits: " << s1 << endl;
    //    cout << "Number of hits that have been filtered (gaps, frameshift, relative score): " << count_not_considered << endl;
    size_t s2 = map_of_vulgar_hits_targets_as_keys.size();
    // TODO: Check what is the difference to s1. The text has been changed. I think s1 used to be without double hits?
    //       If s1 == s2, we can remove one line or indicate.
    cout << "PROGRESS: Number of successful exonerate alignments with multiple hits: " << s2 << endl;
    //    cout << "Number of queryID, targetID hit pairs: " << size_query_target_set(map_of_vulgar_hits_targets_as_keys);
  }

  // Debug output:
  if (global_verbosity >= 100)
  {
    cout << "DEBUG (100): Exonerate results.\n";

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
    seqs_DNA_reads_target.get_full_sequence_names(names);

    for (unsigned i=0; i<names.size(); ++i)
    {
      cout << "DEBUG(100): Name " << i << " " << names[i] << endl;
    }
  }

  /*
   ofstream os;
   os.open(global_alignment_output_file.c_str());
   */


//// Refeactor alignment block to function:

align_sequences_using_vulgar(seqnames_of_references,
                             map_of_vulgar_hits_targets_as_keys,
                             vector_of_hit_stats_for_query_reference,
                             seqs_DNA_reads_target,
                             map_query_prot_lengths,
                             result_alignments_for_references,
                             all_maps_of_insertion_sites_suggested_in_reference,
                             all_maps_of_gap_sites_in_queries);




  if (num_WARNINGS > 0)
  {
    if (global_verbosity >= 1)
      cout << "INFO: There have been " << num_WARNINGS << " warnings or notes. See above for details." << endl;
    else
      cout << "INFO: Due to the verbosity parameter set to 0, you missed "
      << num_WARNINGS
      << " warnings or notes. Increase the verbosity if you want to see the details and run again." << endl;
  }

  cout << PROGNAME " finished successfully." << endl << endl;
} // END: main(...)
