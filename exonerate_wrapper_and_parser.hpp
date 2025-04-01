// CHECK WHETHER WE NEED THEM ALL!!!!!
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
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
#include <filesystem>

inline bool fileExists(std::string filename)
{
  std::ifstream is(filename.c_str());
  if (is.fail())
    return false;
  else
    return true;
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


struct exonerate_parameters
{
  faststring exonerate_binary;
  faststring protein_reference;
  faststring fasta_reads;
  faststring vulgar_folder;
  char genetic_code_number;
  int frameshift_penalty;
  int max_trials;
  char mode;

  exonerate_parameters(faststring exonerate_binary,
                       faststring protein_reference,
                       faststring fasta_reads,
                       faststring vulgar_folder,
                       char genetic_code_number,
                       int frameshift_penalty,
                       int max_trials,
                       char mode=0):exonerate_binary(exonerate_binary),
                                    protein_reference(protein_reference),
                                    fasta_reads(fasta_reads),
                                    vulgar_folder(vulgar_folder),
                                    genetic_code_number(genetic_code_number),
                                    frameshift_penalty(frameshift_penalty),
                                    max_trials(max_trials),
                                    mode(mode)
  {}

  faststring get_fasta_basename_without_extension() const
  {
    return fasta_reads.filename_basename_without_extension();
  }

  faststring get_prot_basename_without_extension() const
  {
    return protein_reference.filename_basename_without_extension();
  }

  faststring get_exonerate_vulgar_filename() const
  {
    faststring vulgarname = get_prot_basename_without_extension() + "_VS_" + get_fasta_basename_without_extension();
    return vulgar_folder + vulgarname + ".vulgar";
  }

  faststring get_exonerate_log_filename() const
  {
    faststring vulgarname = get_prot_basename_without_extension() + "_VS_" + get_fasta_basename_without_extension();
    return vulgar_folder + vulgarname + ".log";
  }

};


class vulgar
{
public:
  static std::map<faststring, short>  queryID2Index;

  faststring queryID;          // Exonerate query sequences are the protein reference sequences, e.g. AA-COI consensus seq. exonerate aligns against.
  short      queryIndex;       // Every query sequence (reference sequence) has a query number.
                               // This avoids many lookups for index.
  unsigned   query_start;      // 0 based numbers.
  unsigned   query_end;        // The first position after the specified range.
  char       query_strand;
  faststring targetID;         // Exonerate targets are the DNA input sequences, typically reads.
  unsigned   target_start;     // 0 based numbers.
  unsigned   target_end;       // The first position after the specified range.
  char       target_strand;
  unsigned   score;
  std::vector<Ctriple<char, unsigned, unsigned> > attributes;


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
    std::map<faststring, short>::iterator find_it = queryID2Index.find(q);
    if (find_it != queryID2Index.end())
      return find_it->second;

    short s = (short) queryID2Index.size();
    queryID2Index[q] = s;
    return s;
  }

  short get_index_in_query2short_map(faststring q)
  {
    std::map<faststring, short>::iterator find_it = queryID2Index.find(q);
    if (find_it != queryID2Index.end())
      return find_it->second;
    else
    {
      std::cerr << "Reference sequence appears in Exonerate output, but not in "
      "reference sequences. This should not be possible.\n";
      good_bye_and_exit(-127);
      return -1; // We never get here. Only introduced to silence warnings about no return for non-void function.
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

    std::vector<faststring> l;
    split(l, str);

    if (l.size() < 13) // 1+9+3
    {

      std::cerr << "ERROR: Bad vulgar string encountered with a wrong number of elements: " << str << std::endl;
      good_bye_and_exit(-4);
    }

    if (l[0] != "vulgar:")
    {
      std::cerr << "ERROR: Bad vulgar string: The vulgar string is expected to start with \"vulgar:\" but found: " << str << std::endl;
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

    size_t m = l.size()-10;
    if (m%3 != 0)
    {
      std::cerr << "ERROR: Bad vulgar string encountered. The number of fields must "
      "by 10+3*n where n is an integer. But the number is "
      << m << "." << std::endl;
      std::cerr << "       Vulgar string: " << str << std::endl;
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
  void print(std::ostream &os) const
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

  size_t num_attributes() const
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




int run_exonerate(const faststring &, const faststring &, const faststring &, const faststring &, const faststring &, char , int , int , char =0);
int run_and_parse_exonerate_for_single_input_file(
          const exonerate_parameters &ep,
          std::vector<vulgar *>  &vec_of_hits_as_in_file,
          std::multimap<faststring, vulgar *> &map_of_vulgar_hits_targets_as_keys,
          std::vector<stats_for_given_target> &vector_of_hit_stats_for_query_reference);
