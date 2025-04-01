// TODO: Check all new whether they have a corresponding delete

// IMPORTANT:
// Some functions only work if the short names are unique. These are the sequence
// names before the first space or the full name if the name contains no spaces.
// Reading and writing work even if the short names are not unique. Accessing
// sequences by the short name does not work.

// Double check results when using this for library for sequences of different
// length. Only very basic tasks are allowed!

//////////////////////////////////////////////////////////////////////
// CSequences.h: interface for the CSequences class.
//
//////////////////////////////////////////////////////////////////////

// Next general goal: Use faststring::size_type instead of unsigned int for array sizes.

#ifndef  CSEQUENCES3_H
#define  CSEQUENCES3_H
#include <vector>
#include <map>
#include "faststring3.h"
//#include "fast-realloc-vector.h"
#include "CSplit2.h"
#include "CSequence_Mol3.1.h"
#include <cassert>
#include <cmath>
#include <utility>
#include "basic-DNA-RNA-AA-routines.h"
#include <climits>
#include <cstdint>
#include <algorithm>

//#include "basic-DNA-RNA-AA-routines.h"

#define PrintMessage_cerr(s1)                fputs(s1, stderr)
#define ErrorMessage(s1)                     fputs(s1, stderr)
#define flush_cerr()                         fflush(stderr)


////////////////////////////////////////////////////////////////////
// Head file for DNA, RNA or PROTEIN molecular sequences.
//
// Characteristics:
// - most member functions assume that sequences have identical lengths, though this is not true for all.
//   If sequences are not of equal length, it is recommended to revise the code in this head file first.
//
////////////////////////////////////////////////////////////////////

// Global Functions
//////////////////////////////////

// Criteria for reference taxa:
// - Original reference taxa have 2 pipes (3 sections), other sequences have 3 pipes (4 sections) in sequence names.
// - Hamstered reference taxa are those for which section 2 and 3 are identical in the first 4 characters.



//template<class T>
inline void add_or_count(std::map<faststring, unsigned> &m, const faststring &x)
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

inline unsigned add_or_count_return(std::map<faststring, unsigned> &m, const faststring &x)
{
  std::map<faststring,unsigned>::iterator it;

  it = m.find(x);
  if (it == m.end() )
  {
    return m[x] = 1;
  }
  else
  {
    return (++it->second);
  }
}


inline void vector_of_faststring_shorten_all(std::vector<faststring> &v, faststring::size_type len)
{
  size_t N = v.size();
  for (size_t i=0; i< N; ++i)
  {
    v[i].shorten(len);
  }
}

// Shortens strings if they are longer then len, appends copies of fill if they are to short.
inline void vector_of_faststring_resize_all(std::vector<faststring> &v, faststring::size_type len, char fill)
{
  size_t N = v.size();
  for (size_t i=0; i< N; ++i)
  {
    v[i].resize(len, fill);
  }
}

inline void vector_of_faststring_replace_characters(std::vector<faststring> &v, faststring &string_of_replacement_symbols, char skip_char)
{
  unsigned char all_symbols_lookup[256];

  for (unsigned j=0; j<256; ++j)
  {
    all_symbols_lookup[j] = (unsigned char)j;
  }

  faststring::size_type n = string_of_replacement_symbols.size();
  for (faststring::size_type i=0; i<n; i+=2)
  {
    all_symbols_lookup[(unsigned char)string_of_replacement_symbols[i]] = string_of_replacement_symbols[i+1];
  }

  for (size_t i=0; i < v.size(); ++i)
  {
    char       *pos     = v[i].begin();
    char       *pos_end = v[i].end();
    char       *pos_write;
    char       c;
    unsigned   removed  = 0;
    size_t     N = v[i].length();

    pos_write = pos;

    while (pos < pos_end)
    {
      c = all_symbols_lookup[(unsigned char)*pos];
      if (c != skip_char)
      {
        *pos_write = c;
        ++pos_write;
      }
      else
      {
        ++removed;
      }
      ++pos;
    }
    v[i].shorten(N-removed);
  }
}


inline void vector_of_faststring_unique(std::vector<faststring> &v, faststring::size_type len)
{
  unsigned digits = (unsigned)log10(v.size()) + 1;
  unsigned mlen;

  if (digits+1 >= len)
    mlen = 0;
  else
    mlen = (unsigned)(len-digits-1);

  std::map<faststring, unsigned>  m_names, m_names_short;
  std::map<faststring, unsigned>::iterator  f_it, f_it2;

  faststring   tmp;

  for (size_t i=0; i < v.size(); ++i)
  {
    add_or_count(m_names, v[i]);
  }

  for (size_t i=0; i < v.size(); ++i)
  {
    if (m_names[v[i]] > 1)
    {
      v[i].shorten(mlen);
    }
    add_or_count(m_names_short, v[i]);
  }
  m_names.clear();

  char s = '_';
  if (mlen == 0)
    s = 'S';

  unsigned n;
  for (unsigned i=0; i < v.size(); ++i)
  {
    if (m_names_short[v[i]] > 1)
    {
      n = add_or_count_return(m_names, v[i]);
      v[i] += s + faststring(n, '0', digits);
    }
  }
}


inline void vector_of_faststring_trimSpaces_front_back(std::vector<faststring> &v)
{
  unsigned N = (unsigned)v.size();
  for (unsigned i=0; i< N; ++i)
  {
    v[i].removeSpacesFront();
    v[i].removeSpacesBack();
  }
}



//////////////////////////////////
// Data types
//////////////////////////////////

//////////////////////////////////
// Class CSequences
//////////////////////////////////

// Headerfile with class for aligned and not aligned DNA, RNA or PROTEIN molecular sequences.

class CSequences3
{
public:
  typedef char chartype;

private:
  SequenceDataTypesEnum	          datatype;            // DNA, Protein, or other
  unsigned                           taxaNum;
  char                               ambig_char;
  bool                               originalPosNumbers_supplied;
  faststring::size_type	             posNum;  // Stores SIZE_MAX if sequences have unequal lengths.
  // Checked for fasta. Phylip input guarantees
  // equal lengths.
  std::vector<faststring::size_type> originalPosNumbers;  // In case we excluded positions from the alignment, we want to remember the original positions.
  std::vector<CSequence_Mol*>        seqData;             // The vector of pointers to sequences

  // Map of short names. Only makes sense if the short names are unique!
  std::map<faststring, CSequence_Mol*> sn_map;       // Obtain pointer to the sequence by sequence name. Is map of short names.
  // Short names are normal names
  bool                      short_names_unique;  // The short_name is the name before the first space or the full name if the name has no spaces.

  void add_seq(CSequence_Mol* seq)
  {
    //    std::cerr << "Adding: " << seq->getName() << '\n';
    if (sn_map.find(seq->getName()) != sn_map.end() )
      short_names_unique = false;

    seqData.push_back(seq);
    sn_map[seq->getName()] = seq;

    if (taxaNum == 0) // This is the old taxaNum, before adding the seq.
      posNum = seq->length();
    else if (posNum > 0 && posNum != seq->length()) // If sequences have different lengths, we set posNum to 0.
      posNum = SIZE_MAX;

    ++taxaNum;
  }


  void remove_seq_del(unsigned seq_num)
  {
    CSequence_Mol *p = seqData[seq_num];
    faststring sname = p->getName();

    --taxaNum;
    seqData.erase(seqData.begin()+seq_num);
    std::map<faststring, CSequence_Mol*>::iterator find_it = sn_map.find(sname);
    if (find_it != sn_map.end() )
      sn_map.erase(find_it);
    delete p;

    if (taxaNum > 0 && equal_length_of_all_sequences())
      posNum = seqData[0]->length();
    else
      posNum = SIZE_MAX;
  }

  CSequence_Mol* remove_seq(unsigned seq_num)
  {
    CSequence_Mol *p = seqData[seq_num];
    faststring sname = p->getName();

    --taxaNum;
    seqData.erase(seqData.begin()+seq_num);
    std::map<faststring, CSequence_Mol*>::iterator find_it = sn_map.find(sname);
    if (find_it != sn_map.end() )
      sn_map.erase(find_it);
    if (taxaNum > 0 && equal_length_of_all_sequences())
      posNum = seqData[0]->length();
    else
      posNum = SIZE_MAX;
    return p;
  }


  void determine_map_of_sequence_names()
  {
    short_names_unique = true;
    // Creates/recreates the map of short names.
    if (sn_map.size() == 0)
    {
      std::vector<CSequence_Mol*>::size_type N=seqData.size();
      for(std::vector<CSequence_Mol*>::size_type i=0; i<N; ++i)
      {
        if (sn_map.find(seqData[i]->getName()) != sn_map.end() )
          short_names_unique = false;
        sn_map[seqData[i]->getName()] = seqData[i];
      }
    }
  }

  void recompute_are_shortnames_unique()
  {
    sn_map.clear();
    determine_map_of_sequence_names();
  }

public:
  // Minimal constructor: Empty sequences object
  CSequences3(SequenceDataTypesEnum Dt, char set_ambig_char='?'):
  datatype(Dt), taxaNum(0), ambig_char(set_ambig_char), originalPosNumbers_supplied(false),posNum(0), short_names_unique(true)
  {}

  // Constructor for a set of empty sequences with names and length.
  CSequences3(SequenceDataTypesEnum Dt, std::vector<faststring> names, faststring::size_type len):
  datatype(Dt), taxaNum(0), ambig_char('?'), originalPosNumbers_supplied(false),posNum(0), short_names_unique(true)
  {
    int n = (int)names.size();

    seqData.reserve(n);
    CSequence_Mol   *seq;

    for (unsigned i=0; i<n; ++i)
    {
      seq = new CSequence_Mol (SeqType_dna, names[i], len, ambig_char);
      add_seq(seq);
    }
  }


  ~CSequences3()
  {
    int n=taxaNum;

    for (int i=0; i<n; ++i)
    {
      delete seqData[i];
    }
  }


  // This constructor can be used as/to
  //  - general copy constructor
  //  - extract a range of sites
  // Coordinates:
  //    pos1 must be the first column. Is 0 based index.
  //    pos2 must be the index after the last column. Is 0 based index.
  //    pos2-pos1 must be the number of bases that are copied to this sequence.

  CSequences3(const CSequences3 &s, faststring::size_type pos1 = 0, faststring::size_type pos2 = faststring::npos):
  datatype(s.datatype), taxaNum(0), ambig_char(s.ambig_char),
  originalPosNumbers_supplied(false), posNum(0), short_names_unique(true)
  {
    unsigned        n=s.taxaNum;
    CSequence_Mol   *seq;

    if (n > 0)
    {
      // This is taken out of the loop to silence the warning that seq is used uninitialised below.
      seq = new CSequence_Mol ( *(s.seqData[0]), pos1, pos2);
      add_seq(seq);

      for (unsigned i=1; i<n; ++i)
      {
        seq = new CSequence_Mol ( *(s.seqData[i]), pos1, pos2);
        add_seq(seq);
      }
      // The following code does not check anything - remove later.
      //    if (i != n)
      //    {
      //      std::cerr << "Critical error in CSequences3 constructor: taxaNum and number of sequences found disagree.\n";
      //    }

      //      posNum = seq->length(); // Is set in add_seq
    }
  }

  bool are_short_names_unique()
  {
    return short_names_unique;
  }

  bool equal_length_of_all_sequences()
  {
    return taxaNum == 0 || posNum != SIZE_MAX;
  }

  bool equal_length_of_all_sequences_compute_again()
  {
    faststring::size_type len;

    if (taxaNum==0)
      return true;

    len = seqData[0]->length();

    for (unsigned i=1; i<taxaNum; ++i)
    {
      if (len != seqData[i]->length())
        return false;
    }
    return true;
  }


  void trim_seq_names(faststring::size_type max_len)
  {
    for (unsigned i=0; i<taxaNum; ++i)
    {
      seqData[i]->trim_seq_name(max_len);
    }
    recompute_are_shortnames_unique();
  }


  void trim_seq_names(const char *trim_at_these_symbols)
  {
    for (unsigned i=0; i<taxaNum; ++i)
    {
      seqData[i]->trim_seq_name(trim_at_these_symbols);
    }
    recompute_are_shortnames_unique();
  }


  unsigned long memory_usage(unsigned long &value1, unsigned long &value2, unsigned long &value3, float &mean2)
  {
    unsigned long tmp1=sizeof(*this);
    unsigned long tmp2=0, tmp3=0;

    unsigned long i,N;
    for (i=0, N=seqData.size(); i<N; ++i)
      tmp2 += seqData[i]->memory_usage();

    std::map<faststring, CSequence_Mol*>::iterator it;
    it = sn_map.begin();
    while (it != sn_map.end() )
    {
      tmp3 += it->first.memory_usage();
      tmp3 += sizeof(it->first);
      tmp3 += sizeof(it->second);
      ++it;
    }

    value1 = tmp1;
    value2 = tmp2;
    value3 = tmp3;
    mean2 = (float)tmp2/N;
    return tmp1 + tmp2 + tmp3;
  }


  CSequence_Mol* get_seq_by_name(faststring name) // Must be short name.
  {
    // Test code:
    {
      //      print_DEBUG(cerr, 0);
    }

    std::map<faststring, CSequence_Mol*>::iterator find_it = sn_map.find(name);
    if (find_it != sn_map.end() )
      return find_it->second;
    else
      return NULL;
  }

  CSequence_Mol* get_seq_by_index(unsigned id)
  {
    if (id >= seqData.size())
    {
      return NULL;
    }
    else
    {
      return seqData[id];
    }
  }

  // Depreciated: Use get_seq_by_index instead
  /*
   CSequence_Mol* get_seq(unsigned id)
   {
   return get_seq_by_index(id);
   }
   */

  const char* get_Seq_Data(unsigned id)
  {
    if (id >= seqData.size())
    {
      return NULL;
    }
    else
    {
      return seqData[id]->getSeqStr();
    }
  }

  const char* get_Seq_Name(unsigned id)
  {
    if (id >= seqData.size())
    {
      return NULL;
    }
    else
    {
      return seqData[id]->getName();
    }
  }

  const char* get_Seq_FullName(unsigned id)
  {
    if (id >= seqData.size())
    {
      return NULL;
    }
    else
    {
      return seqData[id]->getFullName();
    }
  }

  unsigned              GetTaxaNum() { return taxaNum;}
  faststring::size_type GetPosNum()  { return posNum;}

  SequenceDataTypesEnum    get_datatype()
  {
    return datatype;
  }

  chartype GetChar(unsigned TaxaIndex, faststring::size_type PosIndex) const
  {
    assert(TaxaIndex < taxaNum);
    CSequence_Mol *seq = seqData[TaxaIndex];
    assert(PosIndex  < seq->length());
    return seq->get_pos(PosIndex);
  }

  void print_DEBUG(std::ostream &os, unsigned flag=1)
  {
    if (flag & 1) // scalars:
    {
      os << "Debug output CSequences3 object, flag==0\n";
      os << "Data type:  " << (int)datatype<< '\n';
      os << "taxaNum:    " << taxaNum      << '\n';
      os << "posNum:     " << posNum       << '\n';
      os << "ambig_char: " << ambig_char   << '\n';
      os << "originalPosNumbers_supplied " << (int)originalPosNumbers_supplied << '\n';
      os << "size of originalPosNumbers vector: " << originalPosNumbers.size() << '\n';
      os << "size of seqData                  : " << seqData.size()            << '\n';
      os << "size of sn_map                   : " << sn_map.size()             << '\n';
    }

    if (flag & 2)
    {
      faststring::size_type n=seqData.size();

      os << "Data types of sequences:\n";
      for (unsigned i=0; i<n; ++i)
      {
        os << i << ": " << seqData[i]->get_datatype() << " " << seqData[i]->type_as_string() << '\n';
      }
    }
  }



  /*    void             SetChar(unsigned TaxaIndex, unsigned PosIndex, chartype mychar) */
  /*    { */
  /*      assert(TaxaIndex < taxaNum); */
  /*      assert(PosIndex < posNum); */
  /*      seqData[TaxaIndex]->set_pos(PosIndex, mychar); */
  /*    } */

  void SetOriginalPosNumber(faststring::size_type index, faststring::size_type theOriginalPosNumber)
  {
    originalPosNumbers_supplied = true;
    assert(index < posNum);
    originalPosNumbers[index] = theOriginalPosNumber;
  }

  faststring::size_type GetOriginalPosNumber(faststring::size_type index) const
  {
    assert(index < posNum);
    if (originalPosNumbers_supplied)
      return originalPosNumbers[index];
    else
      return UINT_MAX;
  }

  // Moves taxon with index index to top, preserving the order of all other taxa.
  // Of course this changes all indices with the only exception that this taxon is
  // already at the top of the vector.
  // The index is of course 0 based.
  void reorder_move_seq_to_top(unsigned index)
  {
    if (index < taxaNum && index > 0)
    {
      seqData.insert(seqData.begin(), seqData[index]);
      seqData.erase(seqData.begin()+index+1); // We have to add 1 since we have one additional entry at the beginning
    }
  }

  // TODO: Not very efficient!!
  // Return true if taxon_name has been found, false otherwise
  bool reorder_move_seq_to_top(faststring &taxon_name)
  {
    for (unsigned i=0; i< taxaNum; ++i)
    {
      faststring iname = seqData[i]->getFullName();
      if (seqData[i]->getFullName() == taxon_name)
      {
        std::cerr << "Move to top: " << i << '\n';
        reorder_move_seq_to_top(i);
        return true;
      }
    }
    return false;
  }

  // Reorder the sequences such the those in the given vector are at the top
  // in the order given by the vector.
  // TODO: Not very efficient! - Very simple and thus secure implementation.
  bool reorder_sort_by(std::vector<faststring> &names)
  {
    std::vector<faststring*>::size_type  n=names.size();
    bool success = true;

    for (int i=(int)n-1; i>=0; --i)
    {
      success = success & reorder_move_seq_to_top(names[i]);
    }
    return success;
  }

  // Reorder the sequences according to the
  // in the order given by the vector.
  // TODO: Not very efficient! - Very simple and thus secure implementation.

  // Reorder sequences: flag
  // 0: Do not change
  // 1: Alphabetical
  // 2: Alphabetical case insensitive
  // 3: By first unsigned in sequence names
  void reorder_sequences(short flag)
  {
    if (flag == 1)
      sort(seqData.begin(), seqData.end(), less_than_full_seqname_using_pointer);
    else if (flag == 2)
      sort(seqData.begin(), seqData.end(), less_than_full_seqname_caseinsensitive_using_pointer);
    else if (flag == 3)
      sort(seqData.begin(), seqData.end(), less_than_full_seqname_first_unsigned_using_pointer);
  }


  void reorder_sort_by_first_number_in_seqquence_names()
  {
    sort(seqData.begin(), seqData.end(), less_than_full_seqname_first_unsigned_using_pointer);
  }

  bool get_Originalposnumbers_Supplied()
  {
    return originalPosNumbers_supplied;
  }

  char get_ambiguity_character()
  {
    return ambig_char;
  }

  void trimSeqNamesAtSymbolsInString(const char *symbol_list)
  {
    for (unsigned i=0; i < taxaNum; ++i)
    {
      seqData[i]->trim_seq_name(symbol_list);
    }
  }

  // Returns the sequence names upto the first space
  void get_short_sequence_names(std::vector<faststring> &snames)
  {
    snames.clear();
    snames.reserve(taxaNum);
    //     unsigned  max=0;

    for (unsigned i=0; i < taxaNum; ++i)
    {
      snames.push_back(seqData[i]->getName_faststring());
    }
    //     return max;
  }

  // Returns the maximum sequence, length
  void get_full_sequence_names(std::vector<faststring> &snames)
  {
    snames.clear();
    snames.reserve(taxaNum);
    //     unsigned  max=0;

    for (unsigned i=0; i < taxaNum; ++i)
    {
      snames.push_back(seqData[i]->getFullName());
      /*        if (snames[i].size() > max) */
      /*        { */
      /* 	 max = snames[i].size(); */
      /*        } */
    }
    //     return max;
  }


  // Only works for recoded sequences of equal length!!!!
  void     ConsensusSequence(faststring&           conSeq,
                             const CSplit&         setOfTaxa,
                             unsigned char         numDistCharacters, // Distinguishable characters. Recoded characters have indices 0, ..., numDistCharacters-1
                             unsigned char         numSymbols,
                             double                consensusThreshold )
  {
    unsigned char           maxindex, equalindex;
    unsigned                taxaCount = 0;
    unsigned                consensusSymMinimum;

    std::vector<unsigned>  counterSymbolsVec(numSymbols,0);

    if (posNum == SIZE_MAX)
      return;

    for (faststring::size_type pos=0; pos < posNum; ++pos)
    {
      // Initialise variable
      for (unsigned char i=0; i < numDistCharacters; ++i)
        counterSymbolsVec[i] = 0;

      // Count number of occuring symbols for this position over all taxa
      for (unsigned taxon = 0; taxon < taxaNum; ++taxon)
      {
        if (setOfTaxa.test(taxon))
        {
          ++taxaCount;
          ++counterSymbolsVec[seqData[taxon]->get_pos(pos)];
        }
      }

      consensusSymMinimum = (unsigned) std::ceil(consensusThreshold * taxaCount);

      maxindex = 0;
      equalindex = numDistCharacters;
      for (unsigned char i = 1; i < numDistCharacters; ++i)
      {
        if (counterSymbolsVec[i] >= counterSymbolsVec[maxindex])
        {
          if (counterSymbolsVec[i] == counterSymbolsVec[maxindex])
            equalindex = i;
          maxindex = i;
        }
      }
      if (counterSymbolsVec[maxindex] >= consensusSymMinimum &&
          maxindex != equalindex)
        conSeq[pos] = maxindex;
      else // Default value, in case consensus cannot be determined
        conSeq[pos] = numDistCharacters;
    }
  }

  // This version works for any kind of sequence, not only for recoded sequences.
  // See also special version for DNA sequences with additional features.
  void     ConsensusSequence(faststring&           conSeq,
                             double                consensusThreshold,
                             char                  default_char = '\0')
  {
    unsigned          maxindex, equalindex;
    unsigned          taxaCount = 0;
    unsigned          consensusSymMinimum;

    const char        **seq_strs;
    seq_strs = new const char* [taxaNum];
    for (unsigned seq_index=0; seq_index<taxaNum; ++seq_index)
    {
      seq_strs[seq_index] =  get_Seq_Data(seq_index);
    }

    conSeq.clear();

    if (posNum == SIZE_MAX)
      return;

    if (default_char == '\0')
    {
      if (get_datatype() == SeqType_dna && default_char == '\0')
      {
        default_char = 'N';
      }
      else if (get_datatype() == SeqType_protein)
      {
        default_char = 'X';
      }
    }

    //     std::cout << "Default char is: " << default_char << '\n';

    unsigned  *counterSymbolsVec = new unsigned[256];

    for (faststring::size_type pos=0; pos < posNum; ++pos)
    {
      //       std::cout << "--- \n";

      // Initialise variable
      taxaCount = 0;
      std::memset(counterSymbolsVec, 0, 256*sizeof(unsigned));

      unsigned char c;

      // Count number of occuring symbols for this position over all taxa
      for (unsigned taxon = 0; taxon < taxaNum; ++taxon)
      {
        //	char c = toupper_lookup[(unsigned char)seqData[taxon]->get_pos(pos)];
        c = (unsigned char)seq_strs[taxon][pos];
        ++counterSymbolsVec[(unsigned char)c];
        //	 std::cout << "Counting << " << c << " in this column " << counterSymbolsVec[c] << '\n';
      }

      // All gaps
      if (taxaNum == counterSymbolsVec[(unsigned char)'-'])
      {
        conSeq.push_back('-');
      }
      else
      {
        taxaCount = taxaNum - counterSymbolsVec[(unsigned char)'-'];
        consensusSymMinimum = (unsigned) std::ceil(consensusThreshold * taxaCount);

        maxindex = 0;
        equalindex = UINT_MAX;

        for (unsigned char i = 0; i < 254; ++i)
        {
          if (counterSymbolsVec[i] >= counterSymbolsVec[maxindex] && i != '-')
          {
            if (counterSymbolsVec[i] == counterSymbolsVec[maxindex])
              equalindex = i;
            maxindex = i;
          }
        }
        //	 std::cout << "Pos: " << pos << " maxindex: " << (char) maxindex << " equalindex " << (char)equalindex << '\n';
        if (counterSymbolsVec[maxindex] >= consensusSymMinimum && maxindex != equalindex)
          conSeq.push_back(maxindex);
        else // Default value, in case consensus cannot be determined
          conSeq.push_back(default_char);
      }
      //       std::cout << "conSeq: " << conSeq << '\n';
    }
    delete [] counterSymbolsVec;
  }


  void     ConsensusSequence_DNA(faststring&           conSeq,
                                 double                consensusThreshold,
                                 unsigned              minimum_uppercase_coverage,
                                 unsigned              minimum_total_coverage,
                                 char                  default_char = '\0')

  {
    ConsensusSequence_DNA(conSeq, consensusThreshold, minimum_uppercase_coverage, minimum_total_coverage, NULL, default_char='\0');
  }


  void     ConsensusSequence_DNA(faststring&           conSeq,
                                 double                consensusThreshold,
                                 unsigned              minimum_uppercase_coverage,
                                 unsigned              minimum_total_coverage,
                                 std::vector<unsigned> *Pointer_coverage_vec,
                                 char                  default_char = '\0')
  {
    unsigned          taxaCount = 0;
    unsigned          consensusSymMinimum;

    const char        **seq_strs;
    seq_strs = new const char* [taxaNum];
    for (unsigned seq_index=0; seq_index<taxaNum; ++seq_index)
    {
      seq_strs[seq_index] =  get_Seq_Data(seq_index);
    }

    conSeq.clear();
    if (posNum == SIZE_MAX)
      return;

    if (Pointer_coverage_vec)
      Pointer_coverage_vec->clear();

    if (default_char == '\0')
    {
      default_char = 'N';
    }

    //     std::cout << "Default char is: " << default_char << '\n';

    unsigned  *counterSymbolsVec = new unsigned[256];
    std::memset(counterSymbolsVec, 0, 256*sizeof(unsigned));

    for (faststring::size_type pos=0; pos < posNum; ++pos)
    {
      //       std::cout << "--- \n";

      // Initialise variable
      taxaCount = 0;

      unsigned char c;

      // Count number of occuring symbols for this position over all taxa
      for (unsigned taxon = 0; taxon < taxaNum; ++taxon)
      {
        //	char c = toupper_lookup[(unsigned char)seqData[taxon]->get_pos(pos)];
        c = (unsigned char)seq_strs[taxon][pos];
        ++counterSymbolsVec[(unsigned char)c];
        //	 std::cout << "Counting << " << c << " in this column " << counterSymbolsVec[c] << '\n';
      }

      // All gaps
      if (taxaNum == counterSymbolsVec[(unsigned char)'-'])
      {
        conSeq.push_back('-');
      }
      else
      {
        unsigned cA     = counterSymbolsVec[(unsigned char)'A']+counterSymbolsVec[(unsigned char)'a'];
        unsigned cC     = counterSymbolsVec[(unsigned char)'C']+counterSymbolsVec[(unsigned char)'c'];
        unsigned cG     = counterSymbolsVec[(unsigned char)'G']+counterSymbolsVec[(unsigned char)'g'];
        unsigned cT     = counterSymbolsVec[(unsigned char)'T']+counterSymbolsVec[(unsigned char)'t'];
        unsigned ctilde = counterSymbolsVec[(unsigned char)'~'];
        //	 unsigned cN = counterSymbolsVec[(unsigned char)'N']+counterSymbolsVec[(unsigned char)'n'];

        /*
         if (verbosity > 0)
         {
         std::cout << "Consensus at index: " << pos << '\n';
         std::cout << cA << " " << cC << " " << cG << " " << cT << '\n';
         std::cout << counterSymbolsVec[(unsigned char)'A'] << " " << counterSymbolsVec[(unsigned char)'a'] << '\n';
         std::cout << counterSymbolsVec[(unsigned char)'C'] << " " << counterSymbolsVec[(unsigned char)'c'] << '\n';
         std::cout << counterSymbolsVec[(unsigned char)'G'] << " " << counterSymbolsVec[(unsigned char)'g'] << '\n';
         std::cout << counterSymbolsVec[(unsigned char)'T'] << " " << counterSymbolsVec[(unsigned char)'t'] << '\n';
         }
         */

        taxaCount = cA + cC + cG + cT + ctilde;
        consensusSymMinimum = (unsigned) std::ceil(consensusThreshold * taxaCount);

        if (consensusSymMinimum < minimum_total_coverage)
          consensusSymMinimum = minimum_total_coverage;

        if (cA > cC && cA > cG && cA > cT && cA > ctilde)
        {
          if (cA >= consensusSymMinimum && counterSymbolsVec[(unsigned char)'A'] >= minimum_uppercase_coverage)
            conSeq.push_back('A');
          else
            conSeq.push_back(default_char);
        }
        else if (cC > cA && cC > cG && cC > cT && cC > ctilde)
        {
          if (cC >= consensusSymMinimum && counterSymbolsVec[(unsigned char)'C'] >= minimum_uppercase_coverage)
            conSeq.push_back('C');
          else
            conSeq.push_back(default_char);
        }
        else if (cG > cA && cG > cC && cG > cT && cG > ctilde)
        {
          if (cG >= consensusSymMinimum && counterSymbolsVec[(unsigned char)'G'] >= minimum_uppercase_coverage)
            conSeq.push_back('G');
          else
            conSeq.push_back(default_char);
        }
        else if (cT > cA && cT > cC && cT > cG && cT > ctilde)
        {
          if (cT >= consensusSymMinimum && counterSymbolsVec[(unsigned char)'T'] >= minimum_uppercase_coverage)
            conSeq.push_back('T');
          else
            conSeq.push_back(default_char);
        }
        else if (ctilde > cA && ctilde > cC && ctilde > cG && ctilde > cT)
        {
          if (ctilde >= consensusSymMinimum && ctilde >= minimum_uppercase_coverage)
            conSeq.push_back('~');
          else
            conSeq.push_back(default_char);
        }
        else // No consensus, the larger two or more are equal
        {
          conSeq.push_back(default_char);
        }
        /*
         if (verbosity > 0)
         {
         std::cout << "Decision: " << conSeq.back() << '\n';
         }
         */
      } // END else of all gaps
      //       std::cout << "conSeq: " << conSeq << '\n';

      if (Pointer_coverage_vec)
        Pointer_coverage_vec->push_back(taxaCount);

      // Reset the counter array:
      counterSymbolsVec[(unsigned char)'A'] = 0;
      counterSymbolsVec[(unsigned char)'a'] = 0;
      counterSymbolsVec[(unsigned char)'C'] = 0;
      counterSymbolsVec[(unsigned char)'c'] = 0;
      counterSymbolsVec[(unsigned char)'G'] = 0;
      counterSymbolsVec[(unsigned char)'g'] = 0;
      counterSymbolsVec[(unsigned char)'T'] = 0;
      counterSymbolsVec[(unsigned char)'t'] = 0;
      counterSymbolsVec[(unsigned char)'-'] = 0;
      counterSymbolsVec[(unsigned char)'~'] = 0;

    } // END for (pos=0; pos < posNum; ++pos)
    delete [] counterSymbolsVec;
  } // END  void     ConsensusSequence_DNA(...)


  void     WeightedConsensusSequence_DNA(faststring&           conSeq,
                                         float                 consensusThreshold,
                                         unsigned              minimum_uppercase_coverage,
                                         unsigned              minimum_total_coverage,
                                         std::vector<unsigned> *Pointer_coverage_vec,
                                         unsigned char         ends_width_with_lower_weight,
                                         unsigned char         fraction_ends_weight,
                                         char                  default_char = '\0')
  {
    //    size_type            pos;
    //    size_type            taxon;
    unsigned  taxaCount = 0;
    float                  weightCount = 0;
    float                  consensusWeightMinimum;

    char                   *weights;
    const char             **seq_strs;
    int                    count_ends_nuc;

    float                  one_over_fractiion_ends_weight = 1.0/fraction_ends_weight;

    conSeq.clear();
    if (posNum == SIZE_MAX)
      return;

    unsigned  tNum = taxaNum;
    faststring::size_type  pNum = posNum;
    faststring::size_type  weights_size = tNum*pNum;

    weights = new char [weights_size];
    std::memset(weights, 1, weights_size);

    seq_strs = new const char* [tNum];
    for (unsigned seq_index=0; seq_index<tNum; ++seq_index)
    {
      seq_strs[seq_index] =  get_Seq_Data(seq_index);

      count_ends_nuc=0;
      for (faststring::size_type pos=0; pos < pNum; ++pos)
      {
        weights[seq_index * pNum + pos] = fraction_ends_weight;
        if (seq_strs[seq_index][pos] != '-')
          ++count_ends_nuc;
        if (count_ends_nuc == ends_width_with_lower_weight)
          break;
      }

      count_ends_nuc=0;
      for (faststring::size_type pos=pNum-1; pos > 0; --pos) // pos > 0 is correct. We will break before pos == 0.
      {
        weights[seq_index * pNum + pos] = fraction_ends_weight;
        if (seq_strs[seq_index][pos] != '-')
          ++count_ends_nuc;
        if (count_ends_nuc == ends_width_with_lower_weight)
          break;
      }
    }

    if (0)
    {
      std::cerr << "DEBUG OUTPUT: Weight matrix:\n";
      for (unsigned seq_index=0; seq_index<tNum; ++seq_index)
      {
        for (faststring::size_type pos=0; pos < pNum; ++pos)
        {
          std::cerr << ((int)weights[seq_index*pNum+pos])%2;
        }
        std::cerr << '\n';
      }
      std::cerr << '\n';
    }

    if (Pointer_coverage_vec)
      Pointer_coverage_vec->clear();

    if (default_char == '\0')
    {
      default_char = 'N';
    }

    //     std::cout << "Default char is: " << default_char << '\n';

    unsigned  *counterSymbolsVec   = new unsigned[256];
    float     *weightedCoverageVec = new float[256];
    std::memset(counterSymbolsVec,   0, 256*sizeof(unsigned));
    std::memset(weightedCoverageVec, 0, 256*sizeof(float));

    for (faststring::size_type pos=0; pos < pNum; ++pos)
    {
      //       std::cout << "--- \n";

      // Initialise variable
      taxaCount = 0;
      weightCount = 0;

      unsigned char c;

      // Count number of occuring symbols for this position over all taxa
      for (unsigned taxon = 0; taxon < tNum; ++taxon)
      {
        //	char c = (unsigned char)seqData[taxon]->get_pos(pos);
        c = (unsigned char)seq_strs[taxon][pos];
        ++counterSymbolsVec[(unsigned char)c];
        if (weights[taxon * pNum + pos] == fraction_ends_weight) // Use the weights:
          weightedCoverageVec[(unsigned char)c] += one_over_fractiion_ends_weight;
        else
          weightedCoverageVec[(unsigned char)c] += 1;
        //	 std::cout << "Counting << " << c << " in this column " << counterSymbolsVec[c] << '\n';
      }

      // All gaps
      if (taxaNum == counterSymbolsVec[(unsigned char)'-'])
      {
        conSeq.push_back('-');
      }
      else
      {
        unsigned cA     = counterSymbolsVec[(unsigned char)'A']+counterSymbolsVec[(unsigned char)'a'];
        unsigned cC     = counterSymbolsVec[(unsigned char)'C']+counterSymbolsVec[(unsigned char)'c'];
        unsigned cG     = counterSymbolsVec[(unsigned char)'G']+counterSymbolsVec[(unsigned char)'g'];
        unsigned cT     = counterSymbolsVec[(unsigned char)'T']+counterSymbolsVec[(unsigned char)'t'];
        unsigned ctilde = counterSymbolsVec[(unsigned char)'~'];

        float    wA     = weightedCoverageVec[(unsigned char)'A']+counterSymbolsVec[(unsigned char)'a'];
        float    wC     = weightedCoverageVec[(unsigned char)'C']+counterSymbolsVec[(unsigned char)'c'];
        float    wG     = weightedCoverageVec[(unsigned char)'G']+counterSymbolsVec[(unsigned char)'g'];
        float    wT     = weightedCoverageVec[(unsigned char)'T']+counterSymbolsVec[(unsigned char)'t'];
        float    wtilde = weightedCoverageVec[(unsigned char)'~'];

        //	 unsigned cN = counterSymbolsVec[(unsigned char)'N']+counterSymbolsVec[(unsigned char)'n'];
        /*
         if (verbosity > 0)
         {
         std::cout << "Consensus at index: " << pos << '\n';
         std::cout << cA << " " << cC << " " << cG << " " << cT << '\n';
         std::cout << counterSymbolsVec[(unsigned char)'A'] << " " << counterSymbolsVec[(unsigned char)'a'] << '\n';
         std::cout << counterSymbolsVec[(unsigned char)'C'] << " " << counterSymbolsVec[(unsigned char)'c'] << '\n';
         std::cout << counterSymbolsVec[(unsigned char)'G'] << " " << counterSymbolsVec[(unsigned char)'g'] << '\n';
         std::cout << counterSymbolsVec[(unsigned char)'T'] << " " << counterSymbolsVec[(unsigned char)'t'] << '\n';
         }
         */

        taxaCount   = cA + cC + cG + cT + ctilde;
        weightCount = wA + wC + wG + wT + wtilde;

        consensusWeightMinimum = consensusThreshold * weightCount;

        //	if (consensusWeightMinimum < minimum_total_coverage)
        //	  consensusWeightMinimum = minimum_total_coverage;

        if (0)
        {
          std::cerr << "DEBUG call conensus at postions: " << pos << '\n';
          std::cerr << "cA:     " << cA      << '\n';
          std::cerr << "cC:     " << cC      << '\n';
          std::cerr << "cG:     " << cG      << '\n';
          std::cerr << "cT:     " << cT      << '\n';
          std::cerr << "ctilde: " << ctilde  << '\n';
          std::cerr << "wA:     " << wA      << '\n';
          std::cerr << "wC:     " << wC      << '\n';
          std::cerr << "wG:     " << wG      << '\n';
          std::cerr << "wT:     " << wT      << '\n';
          std::cerr << "wtilde: " << wtilde  << '\n';
          std::cerr << "taxaCount:                  " << taxaCount           << '\n';
          std::cerr << "weightCount:                " << weightCount         << '\n';
          std::cerr << "consensusWeightMinimum:     " << consensusWeightMinimum << '\n';
          std::cerr << "minimum_uppercase_coverage: " << minimum_uppercase_coverage << '\n';
          std::cerr << "minimum_total_coverage:     " << minimum_total_coverage << '\n';
        }


        if (wA > wC && wA > wG && wA > wT && wA > wtilde)
        {
          if (wA >= consensusWeightMinimum && counterSymbolsVec[(unsigned char)'A'] >= minimum_uppercase_coverage && cA >= minimum_total_coverage)
            conSeq.push_back('A');
          else
            conSeq.push_back(default_char);
        }
        else if (wC > wA && wC > wG && wC > wT && wC > wtilde)
        {
          if (wC >= consensusWeightMinimum && counterSymbolsVec[(unsigned char)'C'] >= minimum_uppercase_coverage && cC >= minimum_total_coverage)
            conSeq.push_back('C');
          else
            conSeq.push_back(default_char);
        }
        else if (wG > wA && wG > wC && wG > wT && wG > wtilde)
        {
          if (wG >= consensusWeightMinimum && counterSymbolsVec[(unsigned char)'G'] >= minimum_uppercase_coverage && cG >= minimum_total_coverage)
            conSeq.push_back('G');
          else
            conSeq.push_back(default_char);
        }
        else if (wT > wA && wT > wC && wT > wG && wT > wtilde)
        {
          if (wT >= consensusWeightMinimum && counterSymbolsVec[(unsigned char)'T'] >= minimum_uppercase_coverage && cT >= minimum_total_coverage)
            conSeq.push_back('T');
          else
            conSeq.push_back(default_char);
        }
        else if (wtilde > wA && wtilde > wC && wtilde > wG && wtilde > wT)
        {
          if (wtilde >= consensusWeightMinimum && ctilde >= minimum_uppercase_coverage && ctilde >= minimum_total_coverage)
            conSeq.push_back('~');
          else
            conSeq.push_back(default_char);
        }
        else // No consensus, the larger two or more are equal
        {
          conSeq.push_back(default_char);
        }
        /*
         if (verbosity > 0)
         {
         std::cout << "Decision: " << conSeq.back() << '\n';
         }
         */
      } // END else of all gaps
      //       std::cout << "conSeq: " << conSeq << '\n';

      if (Pointer_coverage_vec)
        Pointer_coverage_vec->push_back(taxaCount);

      // Reset the counter array:
      counterSymbolsVec[(unsigned char)'A'] = 0;
      counterSymbolsVec[(unsigned char)'a'] = 0;
      counterSymbolsVec[(unsigned char)'C'] = 0;
      counterSymbolsVec[(unsigned char)'c'] = 0;
      counterSymbolsVec[(unsigned char)'G'] = 0;
      counterSymbolsVec[(unsigned char)'g'] = 0;
      counterSymbolsVec[(unsigned char)'T'] = 0;
      counterSymbolsVec[(unsigned char)'t'] = 0;
      counterSymbolsVec[(unsigned char)'-'] = 0;
      counterSymbolsVec[(unsigned char)'~'] = 0;

      weightedCoverageVec[(unsigned char)'A'] = 0;
      weightedCoverageVec[(unsigned char)'a'] = 0;
      weightedCoverageVec[(unsigned char)'C'] = 0;
      weightedCoverageVec[(unsigned char)'c'] = 0;
      weightedCoverageVec[(unsigned char)'G'] = 0;
      weightedCoverageVec[(unsigned char)'g'] = 0;
      weightedCoverageVec[(unsigned char)'T'] = 0;
      weightedCoverageVec[(unsigned char)'t'] = 0;
      weightedCoverageVec[(unsigned char)'-'] = 0;
      weightedCoverageVec[(unsigned char)'~'] = 0;

    } // END for (pos=0; pos < pNum; ++pos)
    delete [] counterSymbolsVec;
    delete [] weightedCoverageVec;
  } // END  void     ConsensusSequence_DNA(...)


  // For unknown reasons, this version is slower. The differences were introduced with the intention to speed up
  // the loop. In the original (faster) version, we used memset to clear the counterSymbolsVec.
  // Since we loop over the array anyway and compare all values once anyway, it could have been faster
  // to reset only the small number of values that are !=0, but this does not seem to be the case.
  // This version works for any kind of sequence, not only for recoded sequences.
  void     ConsensusSequence_slow(faststring&           conSeq,
                                  double                consensusThreshold,
                                  char                  default_char = '\0')
  {
    unsigned          maxindex, equalindex, maxvalue;
    unsigned          taxaCount = 0;
    unsigned          consensusSymMinimum;

    conSeq.clear();
    if (posNum == SIZE_MAX)
      return;

    if (default_char == '\0')
    {
      if (get_datatype() == SeqType_dna && default_char == '\0')
      {
        default_char = 'N';
      }
      else if (get_datatype() == SeqType_protein)
      {
        default_char = 'X';
      }
    }

    //     std::cout << "Default char is: " << default_char << '\n';

    unsigned  *counterSymbolsVec = new unsigned[256];
    std::memset(counterSymbolsVec, 0, 256*sizeof(unsigned));

    for (faststring::size_type pos=0; pos < posNum; ++pos)
    {
      //       std::cout << "--- \n";

      // Initialise variable
      taxaCount = 0;

      // Count number of occuring symbols for this position over all taxa
      for (unsigned taxon = 0; taxon < taxaNum; ++taxon)
      {
        char c = toupper_lookup[(unsigned char)seqData[taxon]->get_pos(pos)];
        ++counterSymbolsVec[(unsigned char)c];
        //	 std::cout << "Counting << " << c << " in this column " << counterSymbolsVec[c] << '\n';
      }

      // All gaps
      if (taxaNum == counterSymbolsVec[(unsigned char)'-'])
      {
        conSeq.push_back('-');
      }
      else
      {
        taxaCount = taxaNum - counterSymbolsVec[(unsigned char)'-'];
        consensusSymMinimum = (unsigned) std::ceil(consensusThreshold * taxaCount);

        maxindex = 0;
        maxvalue = 0;

        equalindex = UINT_MAX;

        for (unsigned char i = 0; i < 255; ++i)
        {
          if (counterSymbolsVec[i] != 0)
          {
            if (counterSymbolsVec[i] >= maxvalue && i != '-')
            {
              if (counterSymbolsVec[i] == maxvalue)
                equalindex = i;
              maxindex = i;
              maxvalue = counterSymbolsVec[i];
            }
            counterSymbolsVec[i] = 0;
          }
        }
        //	 std::cout << "Pos: " << pos << " maxindex: " << (char) maxindex << " equalindex " << (char)equalindex << '\n';
        if (maxvalue >= consensusSymMinimum && maxindex != equalindex)
          conSeq.push_back(maxindex);
        else // Default value, in case consensus cannot be determined
          conSeq.push_back(default_char);
      }
      //       std::cout << "conSeq: " << conSeq << '\n';
    }
    delete [] counterSymbolsVec;
  }


  void GetSymbolFrequenciesAtPosition(unsigned pos, const CSplit &taxaSet,
                                      unsigned numSymbols,
                                      std::vector<unsigned>& frequencies)
  {
    //  assert(frequencies.size() <= numDistCharacters);

    for (unsigned i=0; i<numSymbols; ++i)
      frequencies[i]=0;

    for (unsigned i=0; i<taxaNum; ++i)
    {
      if (taxaSet.test(i))
      {
        ++frequencies[seqData[i]->get_pos(pos)];
      }
    }
  }

  // This function can be used to read a whole fasta file,
  // or to read a file incrementally. The first_seq and last_seq
  // parameters are always with respect to the current file position.
  // By calling the function again and again with first_seq=0 and
  // last_seq = N, will read N sequences, or less if at the end of the file.
  int read_from_Fasta_File(CFile& infile,
                           CSequence_Mol::processing_flag pflag,
                           unsigned first_seq,   // 1 based
                           unsigned last_seq,    // 1 based inclusive
                           bool report_status_to_cerr=true
                           )
  {
    CSequence_Mol   *seq;
    unsigned        count_seq = 0;         // Effectively one based counter of sequences that have been handled.
    unsigned        count_seqs_stored = 0;


    // TODO
    (void) report_status_to_cerr;

    // Remove all sequences from this object - we want to read a new data set.
    clear(datatype);

    while (!infile.eof())
    {
      ++count_seq;
      if (count_seq > last_seq)
      {
        PrintMessage_cerr("Skipping all remaining sequences in this file.");
        flush_cerr();
	break;
      }

      seq = new CSequence_Mol (datatype);
      seq->readFastaSequence_generic(infile, pflag);

      // PrintMessage_cerr(seq.getName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getFullName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getDescription());PrintMessage_cerr("-\n");
      if ( infile.fail() )
      {
        if (infile.fail_reason1() )
        {
          PrintMessage_cerr("\n\n");
          faststring errormsg;
          errormsg   =  "An error occurred while reading the input file. The data type of the sequence is DNA, but the sequence that is read contains non-DNA symbols.\n";
          errormsg  +=  "File position: line ";
          errormsg  +=  faststring(infile.line());
          ErrorMessage(errormsg.c_str());
          PrintMessage_cerr("\n");
          flush_cerr();
          return -24;
        }
        else if (infile.fail_reason2() )
        {
          PrintMessage_cerr("\n\n");
          faststring errormsg;
          errormsg   =  "An error occurred while reading the input file. The data type of the sequence is protein, but the sequence that is read from the file contains non valid amino acid symbols.\n";
          errormsg  +=  "File position: line ";
          errormsg  +=  faststring(infile.line());
          ErrorMessage(errormsg.c_str());
          PrintMessage_cerr("\n");
          flush_cerr();
          return -24;
        }
        else
        {
          PrintMessage_cerr("\n\n");
          faststring errormsg;
          errormsg   =  "An error occurred while reading the input file. It might not be a valid fasta file.\n";
          errormsg  +=  "File position: line ";
          errormsg  +=  faststring(infile.line());
          ErrorMessage(errormsg.c_str());
          PrintMessage_cerr("\n");
          flush_cerr();
          return -25;
        }
      }

      if ( count_seq < first_seq)
      {
        if (report_status_to_cerr)
        {
          PrintMessage_cerr("Skipping sequence ");
          PrintMessage_cerr(faststring(count_seq).c_str());
          PrintMessage_cerr(": ");
          PrintMessage_cerr(seq->getName());
          PrintMessage_cerr("\n");
          flush_cerr();
        }
        continue;
      }

      if (report_status_to_cerr)
      {
        PrintMessage_cerr("Processing sequence ");
        PrintMessage_cerr(faststring(count_seq).c_str());
        PrintMessage_cerr(": ");
        PrintMessage_cerr(seq->getName());
        PrintMessage_cerr("\n");
        flush_cerr();
      }

      add_seq(seq);
      ++count_seqs_stored;
    } // End while

    if (taxaNum > 0)
    {
      ambig_char = seqData[0]->get_ambiguity_character();
    }

    if (taxaNum == 0 && infile.eof())
    {
      infile.clear(CFile::__fail_flag | infile.rdstate());
      return -1;
    }


    // Check that all sequences have the same length.
    // If they have different lengths, we return with a value of 1.
    // This is not an error. 
    if (posNum == SIZE_MAX )
      return 1;
    
    return 0;
  }  // End read_from_Fasta_File


  // Can read the following format variants:
  // Can read sequential and interleaved phylip files
  // Can read strict and relaxed, i.e. 4 combinations.

  // Cannot read multi phylip files.

  // format_flag: 0 strict phylip standard. Sequence names are exactly 10 characters long.
  //              1 relaxed sequence name length. Sequence names end with first space and can have any length.
  //                Spaces are not allowed in sequence names.


  int read_from_Phylip_File(CFile& infile, unsigned format_flag = 1)
  {
    CSequence_Mol   *seq;

    long local_taxaNum = 0;
    long local_posNum = 0;
    //     unsigned        count_seq;
    //     unsigned        global_sequence_from=0, global_sequence_to=UINT_MAX;
    //     bool            global_gaps_to_Ns=false;
    std::vector<faststring*>::size_type    all_lines_N;
    unsigned                               curr_taxon;

    char read_mode;

    // The allowed symbols from the phylip manual:


    char ch;
    faststring line;
    std::vector<faststring> fsvec;
    faststring *p;

    // Clear the oject and initialize it with the supplied data type
    clear(datatype); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    ch = infile.peek();

    // If the first char is a white space or a digit, we expect this
    // is the line with the number of taxa and seq. positions.
    if ( isspace(ch) || isdigit(ch) )
    {
      infile.getline(line);
      split(fsvec, line);

      if (fsvec.size()==2 && fsvec[0].isAnUnsigned() && fsvec[1].isAnUnsigned() )
      {
        local_taxaNum = fsvec[0].ToUnsigned();
        local_posNum  = fsvec[1].ToUnsigned();

        /* 	 std::cout << "Number of taxa     " << taxaNum << '\n'; */
        /* 	 std::cout << "Number of residues " << posNum  << '\n'; */
      }
    }

    //     std::vector<faststring*> sqs;
    std::vector<faststring*> all_lines;

    /*      if (taxaNum == 0 && posNum == 0) */
    /*      { */
    /*        faststring *pp = new faststring; */
    /*        *pp = line; */
    /*        all_lines.push_back(pp); */
    /*      } */

    // We read the complete file into a vector of strings.
    p = new faststring();
    infile.getline(*p);
    p->removeSpacesFront();
    p->removeSpacesBack();

    while (!infile.fail() )
    {
      all_lines.push_back(p);
      p = new faststring();
      infile.getline(*p);
      p->removeSpacesFront();
      p->removeSpacesBack();
    }

    // The problem is: The phylip format is very difficult
    // to parse. Main problem is to distinguish between
    // sequential and interleaved format.
    // The interleaved format is distinguished as follows:
    // After each block there must be a blank line.
    // So if we find blank lines which are not the last line
    // the format is interleaved. Otherwise it is sequential.

    // Another problem is relaxed versus strict format.
    // In the relaxed format the sequence name and the sequence are separated by white spaces.
    // In the strict format the sequence name is exaclty 10 character long.
    // Example: Name123456A MILLV
    // Can be interpreted as sequence AMILLV or as sequence MILLV with the two different names.
    // In the strit format the following would be legal:
    // NAME123456MILLVWSAL LADDPKHIMV
    // The sequence name would be NAME123456 the rest is the sequence.
    // Note that the format allows spaces in the sequence.
    // The two versions cannot be distinguished unambiguously!!
    // Pragmatic approach: In the presence of one white space island, we assume the relaxed format.
    // Otherwise the strict format.

    // Since empty lines are allowed at the end of the file in all formats we remove them.
    // They cannot be indicative of the interleaved format.

    unsigned i;
    all_lines_N = all_lines.size();

    //     std::cerr << all_lines[all_lines_N-1]->length() << '\n';

    while (all_lines[all_lines_N-1]->length() == 0)
    {
      --all_lines_N;
      delete all_lines[all_lines_N];
      std::vector<faststring*>::iterator it = all_lines.end();
      --it;
      all_lines.erase(it);
    }

    unsigned number_blank_lines=0;

    for (i=0; i<all_lines_N; ++i)
    {
      if (all_lines[i]->length() == 0)
      {
        ++number_blank_lines;
      }
    }

    if (all_lines_N == 0 || number_blank_lines == all_lines_N) // Nothing to do???
      return 0;

    if (number_blank_lines == 0) // sequential format
    {
      read_mode = 's';
      // Example:
      //a          ATGATTCAAC CTCAGACCCT TTTAAATGTA GCAGATAACA GTGGAGCTCG
      //           AAAATTGATG
      //b          ATGATTCAAC CTCAGACCCA TTTAAATGTA GCGGATAACA GCGGGGCTCG
      //           AGAATTGATG

      // This can only be interpreted, if the number
      // of taxa and positions has been specified in
      // the first line.
      // What makes it more complicated
      // is the fact, that the phylip format tolerates
      // numbers and other stuff at the beginning of pure sequence lines.

      faststring  seq_taxon_name;
      faststring  seq_in_one_line;
      //       unsigned    current_length = 0;
      unsigned    curr_line_num=0;

      faststring  *currLine;

      // Return if local_taxaNum or local_posNum are 0.
      if (local_taxaNum == 0 || local_posNum == 0)
      {
        // delete all faststrings in all_lines
        for (unsigned ind=0; ind < all_lines.size(); ++ind)
          delete all_lines[ind];
        return -3;
      }


      // For all taxa
      //
      for (curr_taxon = 0; curr_taxon < local_taxaNum; ++curr_taxon)
      {
        if (curr_line_num >= all_lines_N)
        {
          std::cerr << "Parse error while reading the input file:\nUnexpected end of file in line: " << curr_line_num <<  '\n';
          std::cerr << "Trying to read " << local_taxaNum << " taxa but found only " << curr_taxon << "!\n";

          exit(-45);
        }

        currLine       = all_lines[curr_line_num];
        const char* pos = currLine->c_str();
        //	 currLineLen     = currLine->length();
        seq_taxon_name.clear();
        seq_in_one_line.clear();

        // Copy the sequence name
        if (format_flag == 0)
        {
          for (i=0; i<10 && *pos != '\0'; ++i, ++pos)
          {
            seq_taxon_name.push_back(*pos);
          }
          if (i<10) // Fill name with spaces.
            seq_taxon_name.append(10-i, ' ');
        }
        else
        {
          for (; (*pos != '\0') && (*pos != ' ' && *pos != '\t') ; ++pos)
          {
            seq_taxon_name.push_back(*pos);
          }
        }

        // Skip all spaces:
        while ( *pos != '\0' && (*pos == ' ' || *pos == '\t') )
          ++pos;

        // Read data after sequence name:
        if (*pos)
        {
          // Copy the sequence data:
          {
            // Find first character that shall be copied.
            // pos should point to the first char of the sequence.
            char c;

            // Skip all non-allowed symbols - since we skipped spaces, there should be none.
            c = *pos;
            while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
            {
              ++pos;
              c = *pos;
            }
          }
          // Now pos points to the first char that should be copied:
          while (*pos)
          {
            if (*pos != ' ' && *pos != '\t')
            {
              seq_in_one_line.push_back(*pos);
            }
            ++pos;
          }
        } // END if (*pos) AND read data after sequence name.

        // Read lines until we have the required number of residues:
        // TODO: More error checking could be done here:
        // We could check that only allowed symbols are added.
        // By this we also might be able to detect an error earlier
        while (seq_in_one_line.length() < local_posNum)
        {
          ++curr_line_num;
          if (curr_line_num >= all_lines_N)
          {
            std::cerr << "Parse error while reading the input file. The file has been interpreted as being in sequential phylip format.\n";
            std::cerr << "Several problems can trigger this error: (i) Wrong number of residues in sequences compared to the number specified\n";
            std::cerr << "in the file header. (ii) Wrong number of taxa (sequences) compared to the number specified in the file header.\n";
            std::cerr << "(iii) Missing blank line between blocks of data if this file should be in interleaved format.\n";
            exit(-44);
          }
          currLine       = all_lines[curr_line_num];
          //	   const char* pos = currLine->c_str();
          {
            // Find first character that shall be copied.
            char c;
            pos = currLine->c_str();
            c = *pos;
            while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
            {
              ++pos;
              c = *pos;
            }
          }
          // Now pos points to the first char that should be copied:
          // Here we are in the sequential section, reading all residues until
          // we have filled our sequence.
          while (*pos)
          {
            if (*pos != ' ' && *pos != '\t')
            {
              seq_in_one_line.push_back(*pos);
            }
            ++pos;
          }
        }
        // If we are here, the sequence should have been copied to
        // seq_in_one_line

        seq = new CSequence_Mol (datatype);
        seq->set_taxon_and_sequence(seq_taxon_name, seq_in_one_line, datatype);
        add_seq(seq);

        // The next line should start with a new sequence:
        ++curr_line_num;

      } // End for loop - for all taxa

      if (curr_line_num != all_lines_N)
      {
        std::cerr << "Parse error while reading the input file: This can have several reasons:\n"
        "(i) The number of taxa found in the first block does not match the number specified in the file header.\n"
        "(ii) Since only one data block was found, this file is interpreted as sequential format. If this should be a file\n"
        "in interleaved format, blank line must be added between blocks of data.\n"
        "(iii) It could also be that the number of residues in the sequences is larger than specified in the header, so that additional\n"
        "lines of data are interpreted as additional taxa in the sequential format.\n";
        exit(-67);
      }

    }    // End sequential format
    else // Interleaved format
    {
      read_mode = 'i';
      // In the first block, the first 10 characters are
      // the sequence name. Then the first base starts
      // the sequence.
      // After the first block, the sequence begins/continues at
      // the first base symbol in each line. Thus, numbers and spaces are ignored.
      // Sequence names are not allowed in interleaved blocks.

      // With format_flag == 1 we may have a different number of sequence name characters.
      // But then, no spaces are allowed in the sequence name.

      unsigned     curr_line_num=0;
      faststring   *currLine;

      //       unsigned current_length = 0;
      unsigned curr_taxon = 0;

      faststring *pseq_taxon_name   = NULL;
      faststring *pseq_in_one_line  = NULL;

      std::vector<faststring*> taxonnames;
      std::vector<faststring*> sequences;

      currLine        = all_lines[curr_line_num];
      const char* pos;

      // read first block (interleaved) - until we find a blank line:
      while ( currLine->length() != 0 ) // We read all lines of the first block. Counter: curr_line_num
      {
        pos = currLine->c_str();

        // We will store all pointer to sequences in the vectors defined above,
        // so we do not delete sequences here. We keep on creating new ones.
        pseq_taxon_name  = new faststring;
        pseq_in_one_line = new faststring;

        // Copy the sequence name
        if (format_flag == 0)
        {
          for (i=0; i<10 && *pos != '\0'; ++i, ++pos)
          {
            pseq_taxon_name->push_back(*pos);
          }
          if (i<10) // Fill name with spaces.
            pseq_taxon_name->append(10-i, ' ');
        }
        else
        {
          for (;*pos != '\0' && *pos != ' ' && *pos != '\t'; ++pos)
          {
            pseq_taxon_name->push_back(*pos);
          }
        }

        // Skip all spaces:
        while (*pos == ' ' && *pos != '\0' && *pos != '\t')
          ++pos;

        // Read data after sequence name:
        if (*pos) // we need to read the rest of the line
        {
          // Copy the sequence data:
          {
            char c;

            c = *pos;
            while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
            {
              ++pos;
              c = *pos;
            }
          }
          // Now pos points to the first char that should be copied:
          while (*pos)
          {
            if (*pos != ' ' && *pos != '\t')
              pseq_in_one_line->push_back(*pos);
            ++pos;
          }
        } // if (*pos) AND read data after sequence name.

        // Another taxon has been read: Let us store it:

        taxonnames.push_back(pseq_taxon_name);
        sequences.push_back(pseq_in_one_line);

        //	 std::cerr << "Read partial sequence\n";
        //	 std::cerr << pseq_taxon_name->c_str() << '\n';
        //	 std::cerr << pseq_in_one_line->c_str() << '\n';


        // Read next line with next taxon name
        ++curr_taxon;
        ++curr_line_num;
        if (curr_line_num >= all_lines_N)
          break;
        // If we break, currLine will be set to its new value further down
        currLine       = all_lines[curr_line_num];
      } // End while loop over all non blank lines in first block (interleaved)

      // All remaining blocks still need to be read.

// We returned from this function if either one is 0.
//      if (local_taxaNum == 0 && local_posNum == 0)
//        taxaNum = curr_taxon;

      // Check that the number of taxa that might had been specified agrees with the
      // number of taxa we found now.
      if (local_taxaNum != curr_taxon) // do we have a problem?
      {
        std::cerr << "Parse error while reading the input file:\nThe number of taxa in the first block of the input file ("
        << curr_taxon  << ") does not match the number of taxa specified in the header of the phylip file (" << taxaNum << ")\n";
        exit(-55);
      }

      // This line should be blank: Skip blank lines between the blocks
      while (curr_line_num < all_lines_N && all_lines[curr_line_num]->length()==0)
      {
        ++curr_line_num;
      }

      // Now we need to read all remaining blocks:
      while (true)  // (curr_line_num < all_lines_N)
      {
        // Read the next block:
        for (curr_taxon=0; curr_taxon < local_taxaNum; ++curr_taxon)
        {
          // The current line should be the first line of this block
          if (curr_line_num >= all_lines_N)
          {
            // We are within a block and ran out of lines:
            std::cerr << "Parse error while reading the input file:\nUnexpected end of file. More data has been expected.\n";
            exit(-46);
          }

          currLine       = all_lines[curr_line_num];

          // Find first character that shall be copied.
          char c;
          pos = currLine->c_str();
          c = *pos;
          while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
          {
            ++pos;
            c = *pos;
          }

          // Copy the sequence on this line to our data.
          while (*pos)
          {
            if (*pos != ' ' && *pos != '\t')
              sequences[curr_taxon]->push_back(*pos);
            ++pos;
          }
          ++curr_line_num;
        } // End of for loop which reads the next block

        // The next line(s) should be blank or we should be at the end of the list of lines
        if (curr_line_num < all_lines_N && all_lines[curr_line_num]->length()!=0)
        {
          std::cerr << "Parse error while reading the input file:\nThe phylip file was interpreted to be in the phylip format."
          "Either the number of lines of data in one of the data blocks is not identical to the number of taxa\n"
          "or a blank line might be missing between different data blocks in the interleaved format.\n";
          std::cerr << "This problem was detected in line " << curr_line_num+1 << " of the file, but it could have occurred earlier.\n";
          exit(-77);
        }

        while (curr_line_num < all_lines_N && all_lines[curr_line_num]->length()==0)
        {
          ++curr_line_num;
        }
        if (curr_line_num >= all_lines_N) // we parsed all lines of the data set.
        {
          break;
        }
        // This should be the first line of the next block
        currLine = all_lines[curr_line_num];
        // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      } // End while loop read all remaining blocks

//      if (taxaNum > 0 && posNum == 0)
//        posNum = sequences[0]->length();

      // Now we should have read the sequences
      // and we can delete the data we reserved.
      for (curr_taxon=0; curr_taxon < local_taxaNum; ++curr_taxon)
      {
        seq = new CSequence_Mol (datatype);
        seq->set_taxon_and_sequence(*taxonnames[curr_taxon], *sequences[curr_taxon], datatype);
        add_seq(seq);
        delete taxonnames[curr_taxon];
        delete sequences[curr_taxon];
      }
    } // End interleaved format


    // Do some final error checking

    if (taxaNum == 0)
    {
      // delete all faststrings in all_lines
      for (unsigned ind=0; ind < all_lines.size(); ++ind)
        delete all_lines[ind];
      return 0;
    }

    std::vector<CSequence_Mol*>::size_type len1 = seqData[0]->length();

    if (len1 != posNum)
    {
      std::cerr << "Parse error while reading the input file:\n"
      << "The file was interpreted in the " << (read_mode == 's'? "sequential":"interleaved") << " format.\n"
      "The number of positions specified in the file header is not identical to the number of positons\n"
      "in the first sequence.\n"
      "This can have multiple reasons:\n"
      "In the interleaved format, make sure, taxon names only occur in the first block of sequences.\n"
      "In the interleaved format a blank line must be found between all blocks of data.\n";
      std::cerr << "The number of positions in the first sequence is: " << len1
      << ", while the number specified in the header is " << posNum << ".\n";

      //       std::cerr << seqData[0]->getSeqStr() << '\n';

      exit(-66);
    }

    // Check that the number of positions is correct and equal for all sequences:
    for (curr_taxon=1; curr_taxon < taxaNum; ++curr_taxon)
    {
      seq                        = seqData[curr_taxon];
      faststring::size_type len2 = seq->length();

      if (len1 != len2)
      {
        std::cerr << "Parse error while reading the input file:\nThe number of residues in sequence " << curr_taxon+1 << " differs from the number of residues in sequence 1. The sequence lengths for sequence 1 and sequence " << curr_taxon+1
        <<  " are respectively " << len1 << " and " << len2 << ".\n";
        exit(-67);
      }
    }

    // Check that the data type is equal for all sequences and that
    SequenceDataTypesEnum dt_tmp;

    datatype = seqData[0]->get_datatype();
    for (curr_taxon=1; curr_taxon < taxaNum; ++curr_taxon)
    {
      dt_tmp = seqData[curr_taxon]->get_datatype();

      if (datatype != dt_tmp)
      {
        std::cerr << "Warning: Not all sequences seem to have the same data type.\n";
        datatype = SeqType_mixed;
      }
    }

    char ambig0, ambig_tmp;
    ambig0 = seqData[0]->get_ambiguity_character();
    for (curr_taxon=1; curr_taxon < taxaNum; ++curr_taxon)
    {
      ambig_tmp = seqData[curr_taxon]->get_ambiguity_character();

      if (ambig0 != ambig_tmp)
      {
        std::cerr << "Warning: Not all sequences seem to have the same ambiguity character.\n";
      }
    }




    /*      // PrintMessage_cerr(seq.getName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getFullName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getDescription());PrintMessage_cerr("-\n"); */
    /*      if ( infile.fail() && !infile.eof() ) */
    /*      { */
    /*        PrintMessage_cerr("\n\n"); */
    /*        faststring errormsg; */
    /*        errormsg   =  "An error occurred while reading the input file. It might not be a valid fasta file.\n"; */
    /*        errormsg  +=  "File position: line "; */
    /*        errormsg  +=  faststring(infile.line()).c_str(); */
    /*        ErrorMessage(errormsg.c_str()); */
    /*        PrintMessage_cerr("\n"); */
    /*        flush_cerr(); */
    /*        return -25; */
    /*      } */

    /*      if ( count_seq < global_sequence_from || count_seq > global_sequence_to ) */
    /*      { */
    /*        if (false) */
    /*        { */
    /* 	 PrintMessage_cerr("Skipping sequence "); */
    /* 	 PrintMessage_cerr(faststring(count_seq).c_str()); */
    /* 	 PrintMessage_cerr(": "); */
    /* 	 PrintMessage_cerr(seq->getName()); */
    /* 	 PrintMessage_cerr("\n"); */
    /* 	 flush_cerr(); */
    /*        } */
    /*        continue; */
    /*      } */

    /*      PrintMessage_cerr("Processing sequence "); */
    /*      PrintMessage_cerr(faststring(count_seq).c_str()); */
    /*      PrintMessage_cerr(": "); */
    /*      PrintMessage_cerr(seq->getName()); */
    /*      PrintMessage_cerr("\n"); */
    /*      flush_cerr(); */

    /*      seqData.push_back(seq); */
    /*    } // End while */

    for (unsigned j=0; j<all_lines.size(); ++j)
      delete all_lines[j];

    return 0;
  }  // End read_from_Phylip_File

  // This function does not fill in gaps and starts that might have been removed before.
  void ExportSequences(FILE *ofp, char format, unsigned interleaved_len)
  {
    if (taxaNum == 0)
      return;

    if (format=='f')
    {
      for (unsigned i=0; i<taxaNum; ++i)
      {
        seqData[i]->writeSequence_fasta(ofp, interleaved_len);
      }
    }
    else
    {
      std::cerr << "INTERNAL ERROR: This is a specialized export function and only the fasta format is currently supported.\n";
      exit(-111);
    }
  }


  // Backwards compatible and convenience function ExportSequences with fewer parameters.
  void ExportSequences(std::ostream &os, char format, unsigned interleaved_len)
  {
    faststring tmp;
    ExportSequences(os, format, interleaved_len, tmp, UINT_MAX, true);
  }

  // Normal export to different formats:
  //  void ExportSequences(std::ostream &os, char format, unsigned interleaved_len, faststring &replace_symbols_in_sequence_names, unsigned trim_seq_names_length, const char *trim_at_symobls, bool unique_names)
  void ExportSequences(std::ostream &os, char format, unsigned interleaved_len, faststring replace_symbols_in_sequence_names, unsigned trim_seq_names_length, bool unique_names)
  {
    if (taxaNum == 0)
      return;

    std::vector<faststring> vofs;

    unsigned i;

    std::vector<faststring> sequence_names;

    for (i=0; i< taxaNum; ++i)
    {
      vofs.push_back(faststring((faststring::size_type)100, (char)'\0'));
      vofs[i] = seqData[i]->getSeq_faststring();
      sequence_names.push_back(seqData[i]->getFullName());
    }

    // if format is strict phylip, set trim_seq_names_length to 10.
    if (format == 'p')
    {
      //       unique_names = true;
      trim_seq_names_length = 10;
    }

    if (format == 'r')
    {
      bool space_in_query = false;;

      for (faststring::size_type j=0; j<replace_symbols_in_sequence_names.size(); j+=2)
      {
        if (replace_symbols_in_sequence_names[j] == ' ')
          space_in_query = true;
        if (replace_symbols_in_sequence_names[j+1] == ' ') // This target is not allowed.
          replace_symbols_in_sequence_names[j+1] = '_';
      }
      if (!space_in_query)   // We append a replacement for ' '.
        replace_symbols_in_sequence_names += " _";
    }

    if ( !replace_symbols_in_sequence_names.empty() )
    {
      vector_of_faststring_replace_characters(sequence_names, replace_symbols_in_sequence_names, '>');
    }

    if (trim_seq_names_length < UINT_MAX)
    {
      vector_of_faststring_shorten_all(sequence_names, trim_seq_names_length);
    }

    vector_of_faststring_trimSpaces_front_back(sequence_names);

    if (unique_names)
    {
      vector_of_faststring_unique(sequence_names, trim_seq_names_length);
    }

    if (format=='r') // relaxed phylip
    {
      unsigned pos=0;
      unsigned i;

      if (posNum == SIZE_MAX)
      {
          std::cerr << "ERROR while exporting sequences. Exporting to phylip format (format \'r\')"
	               " is not allowed and not possible for sequences of unequal lengths.\n";
          exit(-23);
      }

      os << "\t" << taxaNum << "\t" << posNum << '\n';

      // write first block with sequence names:
      if (pos < posNum)
      {
        for (i=0; i<taxaNum; ++i)
        {
          os << sequence_names[i] << " ";
          os << vofs[i].substr(pos, interleaved_len) << '\n';
        }
        os << '\n'; // Separate blocks
        pos += interleaved_len;
      }

      while (pos < posNum)
      {
        for (i=0; i<taxaNum; ++i)
        {
          os << vofs[i].substr(pos, interleaved_len) << '\n';
        }
        os << '\n'; // Separate blocks
        pos += interleaved_len;
      }
    }
    else
      if (format=='p') // strict phylip
      {
	if (posNum == SIZE_MAX)
	{
          std::cerr << "ERROR while exporting sequences. Exporting to phylip format (format \'p\')"
	               " is not allowed and not possible for sequences of unequal lengths.\n";
          exit(-23);
        }

        unsigned pos=0;
        unsigned i;

        for (i=0; i<taxaNum; ++i)
        {
          sequence_names[i].fill_if_shorter(10, ' ');
        }

        os << "\t" << taxaNum << "\t" << posNum << '\n';

        // write first block with sequence names:
        if (pos < posNum)
        {
          for (i=0; i<taxaNum; ++i)
          {
            os << sequence_names[i] << " ";
            os << vofs[i].substr(pos, interleaved_len) << '\n';
          }
          os << '\n'; // Separate blocks
          pos += interleaved_len;
        }

        while (pos < posNum)
        {
          for (i=0; i<taxaNum; ++i)
          {
            os << vofs[i].substr(pos, interleaved_len) << '\n';
          }
          os << '\n'; // Separate blocks
          pos += interleaved_len;
        }
      }
      else if (format=='f')
      {
        unsigned pos=0;
        unsigned i;

        for (i=0; i<taxaNum; ++i)
        {
          os << ">" << sequence_names[i] << '\n';
          pos = 0;
          faststring::size_type N_len =  vofs[i].size();
          while (pos < N_len)
          {
            os << vofs[i].substr(pos, interleaved_len) << '\n';
            pos += interleaved_len;
          }
        }
      }
      else if (format=='n') // Nexus format with data block:
      {
        /*        Begin data; */
        /*        Dimensions ntax=4 nchar=15; */
        /*        Format datatype=dna symbols="ACTG" missing=? gap=-; */
        /*        Matrix */
        /* 	 Species1   atgctagctagctcg */
        /* 	 Species2   atgcta??tag-tag */
        /* 	 Species3   atgttagctag-tgg */
        /* 	 Species4   atgttagctag-tag            */
        /* 	 ; */
        /*        End; */

        if (posNum == SIZE_MAX)
        {
          std::cerr << "ERROR while exporting sequences to nexus format (format \'n\') which is not possible "
	               "for sequences of unequal lengths.\n";
          exit(-23);
        }

        os << "Begin data;\n";
        os << "Dimensions ntax="<< taxaNum << " nchar="<< posNum << ";" << '\n';
        os << "Format datatype="<< seqData[0]->type_as_string()
        //  << " symbols=\" "<< <<  "\""
        << " missing=" << ambig_char << " gap=-;\n";
        os << "matrix\n";

        unsigned pos=0;
        unsigned i;

        while (pos < posNum)
        {
          for (i=0; i<taxaNum; ++i)
          {
            //           os << "'" << sequence_names[i] << "'" << " ";
            os << sequence_names[i] << " ";
            os << vofs[i].substr(pos, interleaved_len) << '\n';
          }
          pos += interleaved_len;
        }
        os << ";" << '\n' << "end;\n";
      }
      else if (format == 'N') // Nexus format with taxa and characters block:
      {
        unsigned pos;
        unsigned i;

        os << "Begin taxa;\n";
        os << "Dimensions ntax="<< taxaNum << ";\n";
        os << "TAXLABELS ";

        for (i=0; i<taxaNum-1; ++i)
        {
          os << "'" << sequence_names[i] << "'" << " ";
        }
        os << "'" << sequence_names[i] << "'" << ";" << '\n' << "end;\n\n";

        os << "Begin characters;\n";
        os << "Dimensions ntax="<< taxaNum << " nchar="<< posNum << ";\n";
        os << "Format datatype="<< seqData[0]->type_as_string()
        //  << " symbols=\" "<< <<  "\""
        << " missing=" << ambig_char << " gap=-;\n";
        os << "matrix\n";

        pos=0;
        while (pos < posNum)
        {
          for (i=0; i<taxaNum; ++i)
          {
            os << "'" << sequence_names[i] << "'" << " ";
            os << vofs[i].substr(pos, interleaved_len) << '\n';
          }
          pos += interleaved_len;
        }
        os << ";" << '\n' << "end;\n";
      } // END format == 'N'
      else
      {
        std::cerr << "Internal error in ExportSequences(std::ostream, char, unsigned):\n"
        "Unknown value for format parameter. No output has been generated.\n";
      }
  }

  // Export do different formats with site-filter
  // Commented out since not finished and not tested.
  /*    void ExportSequences(std::ostream &os, char format, unsigned interleaved_len, */
  /* 			std::vector<unsigned> site_filter) */
  /*    { */
  /*      // TODO: Currently only supports fasta without interleaved_len */

  /*      if (taxaNum == 0) */
  /*        return; */

  /*      // vector of sequences with gaps and stars. */
  /*      // Why do we need a vector here?? 9.9.2912 */
  /*      std::vector<faststring> vofs; */

  /*      unsigned i; */

  /*      for (i=0; i< taxaNum; ++i) */
  /*      { */
  /*        vofs.push_back(faststring('\0', (faststring::size_type)10000)); */
  /*        seqData[i]->getSequence_fill_in_gaps_and_stars(vofs[i]); */
  /*      } */

  /*      if (format=='p') // Phylip */
  /*      { */
  /* /\*        std::vector<faststring> uphynames; *\/ */
  /* /\*        unique_maxlen_names(uphynames); *\/ */

  /* /\*        unsigned pos=0; *\/ */
  /* /\*        unsigned i; *\/ */

  /* /\*        os << "\t" << taxaNum << "\t" << posNum << '\n'; *\/ */

  /* /\*        while (pos < posNum) *\/ */
  /* /\*        { *\/ */
  /* /\* 	 for (i=0; i<taxaNum; ++i) *\/ */
  /* /\* 	 { *\/ */
  /* /\*            os << uphynames[i] << " "; *\/ */
  /* /\*            os << vofs[i].substr(pos, interleaved_len) << '\n'; *\/ */
  /* /\* 	 } *\/ */
  /* /\* 	 os << '\n'; *\/ */
  /* /\* 	 pos += interleaved_len; *\/ */
  /* /\*        } *\/ */
  /*      } */
  /*      else if (format=='f') */
  /*      { */
  /*        unsigned pos=0; */
  /*        unsigned i; */

  /*        for (i=0; i<taxaNum; ++i) */
  /*        { */
  /*          os << ">" << seqData[i]->getFullName() << '\n'; */
  /* 	 pos = 0; */
  /* 	 faststring &seq_ref = vofs[i]; */
  /* 	 while (pos < posNum) */
  /* 	 { */
  /* 	   if (site_filter[pos] != 0) */
  /* 	     os << seq_ref[pos] << '\n'; */
  /* 	   ++pos; */
  /* 	 } */
  /*        } */
  /*      } */
  /*      else if (format=='n') */
  /*      { */
  /* /\*        Begin data; *\/ */
  /* /\*        Dimensions ntax=4 nchar=15; *\/ */
  /* /\*        Format datatype=dna symbols="ACTG" missing=? gap=-; *\/ */
  /* /\*        Matrix *\/ */
  /* /\* 	 Species1   atgctagctagctcg *\/ */
  /* /\* 	 Species2   atgcta??tag-tag *\/ */
  /* /\* 	 Species3   atgttagctag-tgg *\/ */
  /* /\* 	 Species4   atgttagctag-tag            *\/ */
  /* /\* 	 ; *\/ */
  /* /\*        End; *\/ */

  /* /\*        os << "Begin data;" << '\n'; *\/ */
  /* /\*        os << "Dimensions ntax="<< taxaNum << " nchar="<< posNum << ";" << '\n'; *\/ */
  /* /\*        os << "Format datatype="<< seqData[0]->type_as_string() *\/ */
  /* /\* 	 //  << " symbols=\" "<< <<  "\""  *\/ */
  /* /\* 	 << " missing=" << ambig_char << " gap=-;" << '\n'; *\/ */
  /* /\*        os << "matrix" << '\n'; *\/ */

  /* /\*        unsigned pos=0; *\/ */
  /* /\*        unsigned i; *\/ */

  /* /\*        while (pos < posNum) *\/ */
  /* /\*        { *\/ */
  /* /\* 	 for (i=0; i<taxaNum; ++i) *\/ */
  /* /\* 	 { *\/ */
  /* /\*            os << "'" << seqData[i]->getName() << "'" << " "; *\/ */
  /* /\*            os << vofs[i].substr(pos, interleaved_len) << '\n'; *\/ */
  /* /\* 	 } *\/ */
  /* /\* 	 pos += interleaved_len; *\/ */
  /* /\*        } *\/ */
  /* /\*        os << ";" << '\n' << "end;" << '\n'; *\/ */

  /*      } */

  /*    } */


  void ExportSequences_no_fill_in_ext_phylip_range(std::ostream &os, unsigned begin_exp, unsigned end_exp)
  {
    if (taxaNum == 0)
      return;

    unsigned i;
    faststring::size_type len;

    os << "\t" << taxaNum << "\t" << end_exp-begin_exp << '\n';

    for (i=0; i<taxaNum; ++i)
    {
      len = end_exp-begin_exp;
      if (len > seqData[i]->length() )
        len = seqData[i]->length();
      os.write(seqData[i]->getName(), seqData[i]->getName_faststring().length());
      os.put(' ');
      os.write(seqData[i]->getSeqBegin()+begin_exp, len);
      os.put('\n');
    }
  }

  // Needs to be revised:
  /*
   void ExportSequences_no_fill_in_ext_phylip_range(std::ostream &os, unsigned *coords, unsigned num_coord_pairs)
   {
   if (taxaNum == 0)
   return;

   // if (end_exp > posNum || end_exp <= begin_exp)
   //       return;

   unsigned i,j, N=0;
   unsigned begin_exp;
   unsigned end_exp;


   for (j=0; j<num_coord_pairs; ++j)
   {
   begin_exp = coords[2*j];
   end_exp   = coords[2*j+1];
   if (end_exp > posNum)
   end_exp = posNum;
   if (begin_exp < end_exp)
   N += end_exp-begin_exp;
   }

   os << "\t" << taxaNum << "\t" << N << '\n';

   for (i=0; i<taxaNum; ++i)
   {
   os.write(seqData[i]->getName(), seqData[i]->getName_faststring().length());
   os.put(' ');
   for (j=0; j<num_coord_pairs; ++j)
   {
   begin_exp = coords[2*j];
   end_exp   = coords[2*j+1];
   if (end_exp > posNum)
   end_exp = posNum;
   if (begin_exp < end_exp)
   os.write(seqData[i]->getSeqBegin()+begin_exp, end_exp-begin_exp);
   }
   os.put('\n');
   }
   }
   */


  // Not yet well tested!!!!!
  /*
   void Export_Single_Sequence(std::ostream &os, const faststring& fullname,
   char format, unsigned interleaved_len)
   {
   if (taxaNum == 0)
   return;

   vector of sequences with gaps and stars.
   Why do we need a vector here?? 9.9.2912
   std::vector<faststring> vofs;

   unsigned i;

   for (i=0; i< taxaNum; ++i)
   {
   vofs.push_back(faststring());
   if ( seqData[i]->getFullName() == fullname)
   seqData[i]->getSequence_fill_in_gaps_and_stars(vofs[i]);
   }

   if (format=='p') Phylip
   {
   std::vector<faststring> uphynames;
   unique_maxlen_names(uphynames);

   unsigned pos=0;
   unsigned i;

   for (i=0; i<taxaNum; ++i)
   uphynames[i].fill_if_shorter(10, ' ');

   os << "\t" << taxaNum << "\t" << posNum << '\n';

   while (pos < posNum)
   {
   for (i=0; i<taxaNum; ++i)
   {
   if (!vofs[i].empty())
   {
   os << uphynames[i] << " ";
   os << vofs[i].substr(pos, interleaved_len) << '\n';
   }
   }
   os << '\n';
   pos += interleaved_len;
   }
   }
   else if (format=='f')
   {
   unsigned pos=0;
   unsigned i;

   for (i=0; i<taxaNum; ++i)
   {
   if (!vofs[i].empty())
   {
   os << ">" << seqData[i]->getFullName() << '\n';
   pos = 0;
   while (pos < posNum)
   {
   os << vofs[i].substr(pos, interleaved_len) << '\n';
   pos += interleaved_len;
   }
   }
   }
   }
   else if (format=='n')
   {
   Begin data;
   Dimensions ntax=4 nchar=15;
   Format datatype=dna symbols="ACTG" missing=? gap=-;
   Matrix
   Species1   atgctagctagctcg
   Species2   atgcta??tag-tag
   Species3   atgttagctag-tgg
   Species4   atgttagctag-tag
   ;
   End;

   os << "Begin data;" << '\n';
   os << "Dimensions ntax="<< 1 << " nchar="<< posNum << ";" << '\n';
   os << "Format datatype="<< seqData[0]->type_as_string()
   << " symbols=\" "<< <<  "\""
   << " missing=" << ambig_char << " gap=-;" << '\n';
   os << "matrix" << '\n';

   unsigned pos=0;
   unsigned i;

   while (pos < posNum)
   {
   for (i=0; i<taxaNum; ++i)
   {
   if (!vofs[i].empty())
   {
   os << seqData[i]->getName() << " ";
   os << vofs[i].substr(pos, interleaved_len) << '\n';
   }
   }
   pos += interleaved_len;
   }
   os << ";" << '\n' << "end;" << '\n';

   }


   }

   */



  //   void            GetBaseFrequencies(unsigned, std::vector<unsigned>);

  //   void            initSplits();
  //   void            computeSplits();
  //   CSplits         *splits;
  //   CBasefreq       *basefreq;

  // Changed 29.5.2019: Renamed and removed the clear since it prevents an initialisation of the map for all patterns.
  void compute_site_pattern_map_no_clear(std::map<faststring, unsigned> &m, bool remove_gap_ambig_sites=true)
  {
    //     m.clear();

    unsigned i, j;

    // Here we will store pointers to all sequences for fast access:
    const char **tax = new const char* [taxaNum];

    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    // Store pointers to all sequences in the tax array - this will be faster
    for (j=0; j < taxaNum; ++j)
    {
      tax[j] = seqData[j]->getSeqStr();
    }

    faststring tmp;

    bool contains_gap_ambig;
    char c;

    // For all positions in alignment
    for (i=0; i<posNum; ++i)
    {
      tmp.clear();
      contains_gap_ambig = false;

      // Old and slow version:
      // Copy all symbols to site pattern vector
      /*        for (j=0; j < taxaNum; ++j) */
      /*        { */
      /* 	 c = tax[j][i]; */
      /*  	 if (remove_gap_ambig_sites) */
      /*  	 { */
      /* 	   // Check whether this is a position we do not want to consider. */
      /*  	   if (datatype == CSequence_Mol::dna && (is_DNA_ambig(c) || c == '-') ) */
      /*  	   { */
      /* 	     contains_gap_ambig = true; */
      /*  	   } */
      /* 	   else */
      /* 	   { */
      /* 	     if (datatype == CSequence_Mol::protein && (is_aa_ambig(c) || c == '-') ) */
      /* 	     { */
      /* 	       contains_gap_ambig = true; */
      /* 	     } */
      /* 	   } */
      /*  	 } // END if (remove_gap_ambig_sites) */
      /* 	 tmp.push_back(c); */
      /*        } // END: For all taxa */

      // Faster version:
      if (remove_gap_ambig_sites && datatype == SeqType_dna)
      {
        for (j=0; j < taxaNum; ++j) // Copy all symbols to site pattern vector
        {
          c = tax[j][i];
          if (is_DNA_ambig(c) || c == '-')
          {
            contains_gap_ambig = true;
          }
          tmp.push_back(c);
        }
      }
      else if (remove_gap_ambig_sites && datatype == SeqType_protein)
      {
        for (j=0; j < taxaNum; ++j) // Copy all symbols to site pattern vector
        {
          c = tax[j][i];
          if (is_aa_ambig(c) || c == '-')
          {
            contains_gap_ambig = true;
          }
          tmp.push_back(c);
        }
      }
      else
      {
        for (j=0; j < taxaNum; ++j) // Copy all symbols to site pattern vector
        {
          c = tax[j][i];
          tmp.push_back(c);
        }
      }
      if (!(remove_gap_ambig_sites && contains_gap_ambig))
        add_or_count(m, tmp);
    } // END: for all positions in alignment
  }


  double compute_site_pattern_stat(bool remove_gap_ambig_sites=true)
  {
    std::map<faststring, unsigned> m;
    compute_site_pattern_map_no_clear(m, remove_gap_ambig_sites);

    std::map<faststring, unsigned>::iterator it, it_end;
    it     = m.begin();
    it_end = m.end();

    double   res = 0;
    double   x;
    unsigned N   = 0;

    // Formula from "Statistical Tests of Models of DNA Substitution", Journal of Molecular Evolution, Springer-Verlag NewYork Inc. 1993
    // Nick Goldman
    // Formula: \sum   N_i*ln(N_i) - N*ln(N), where N_i are the pattern frequencies and N is the total number of sites.

    while (it != it_end)
    {
      x = it->second;
      N += x;
      if (x > 0)
        res += x*log(x);
      //       std::cout << x << " " << res << '\n';
      ++it;
    }
    res -= N*log(N);
    //     std::cout << x << " " << res << '\n';

    return res;
  } // END compute_site_pattern_stat(bool remove_gap_ambig_sites=true)


  double compute_site_pattern_stat(std::map<faststring, unsigned> &m, bool remove_gap_ambig_sites=true)
  {
    m.clear();
    compute_site_pattern_map_no_clear(m, remove_gap_ambig_sites);

    std::map<faststring, unsigned>::iterator it, it_end;
    it     = m.begin();
    it_end = m.end();

    double res=0;
    double x;
    unsigned N=0;

    while (it != it_end)
    {
      x = it->second;
      N += x;
      if (x > 0)
        res += x*log(x);
      //       std::cout << x << " " << res << '\n';
      ++it;
    }
    res -= N*log(N);
    //     std::cout << x << " " << res << '\n';

    return res;
  }

  // Determines unfiltered SNP sites:
  void determine_SNP_site(std::map<unsigned, std::pair<int, faststring> > &m, bool remove_gap_ambig_sites=true, int ref_seq=0)
  {
    m.clear();
    unsigned i, j;
    int      pos_in_ref;

//    unsigned gap_ambig_site_count = 0;
//    unsigned non_gap_ambig_site_count = 0;

    // Here we will store pointers to all sequences for fast access:
    const char **tax = new const char* [taxaNum];

    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    // Store pointers to all sequences in the tax array - this will be faster
    for (j=0; j < taxaNum; ++j)
    {
      tax[j] = seqData[j]->getSeqStr();
    }

    faststring pattern;
    bool contains_gap_ambig;
    bool is_SNP_site;
    char c0, c;
    bool is_ref_gap_site;

    // For all positions in alignment
    for (i=0, pos_in_ref=1; i<posNum; ++i)
    {
      //      pattern.clear();
      contains_gap_ambig = false;
      is_SNP_site        = false;
      is_ref_gap_site    = false;

      // Faster version:
      if (datatype == SeqType_dna)
      {
        c0 = tax[0][i];
        if (tax[ref_seq][i] != '-')
          ++pos_in_ref;
        else
          is_ref_gap_site = true;

        if (is_DNA_ambig(c0) || c0 == '-')
        {
          contains_gap_ambig = true;
        }

        for (j=1; j < taxaNum; ++j)
        {
          c = tax[j][i];
          if (is_DNA_ambig(c) || c == '-')
          {
            contains_gap_ambig = true;
          }
          if (c0 != c)
          {
            is_SNP_site = true;
          }
        }
      }
      else if (datatype == SeqType_protein)
      {
        c0 = tax[0][i];
        if (tax[ref_seq][i] != '-')
          ++pos_in_ref;
        else
          is_ref_gap_site = true;

        if (is_DNA_ambig(c0) || c0 == '-')
        {
          contains_gap_ambig = true;
        }

        for (j=1; j < taxaNum; ++j)
        {
          c = tax[j][i];
          if (is_aa_ambig(c) || c == '-')
          {
            contains_gap_ambig = true;
          }
          if (c0 != c)
          {
            is_SNP_site = true;
          }
        }
      }
      else
      {
        std::cout << "ERROR: Incompatible data type in determine_SNP_site member function of "
        "the CSequences class. Exiting.\n";
        std::exit(0);
      }

      std::pair<int, faststring> index_pattern_pair;
      if (is_SNP_site && !(contains_gap_ambig && remove_gap_ambig_sites) )
      {
        pattern.clear();
        for (j=0; j < taxaNum; ++j) // Copy all symbols to site pattern vector
        {
          c = tax[j][i];
          pattern.push_back(c);
        }

        // Map entries are:
        // If position has no gap in reference:
        // (0-based index in alignment, (         1-based index in reference, pattern))
        // If position has gap in reference:
        // (0-based index in alignment, (negative 1-based index in reference, pattern))
        // 1-based index <=0 indicate gap in reference. 0 for first SNPs at beginning with gap in ref.

        if (is_ref_gap_site)
          index_pattern_pair = std::make_pair(-(pos_in_ref-1), pattern);
        else
          index_pattern_pair = std::make_pair(pos_in_ref-1, pattern);
        m.insert(m.end(), std::make_pair(i, index_pattern_pair));
      }
    } // END: for all positions in alignment
  }

  void get_partial_pattern(unsigned pos, unsigned n, faststring &pat)
  {
    unsigned i;

    if (n > taxaNum)
      n = taxaNum;

    pat.clear();
    for(i=0; i<n; ++i)
      pat.push_back((seqData[i]->getSeqStr())[pos]);
  }


  // Reports the vector of frequencies. Does not allow to get information on the pattern for the frequency.
  void get_vector_of_site_pattern_frequencies(std::vector<unsigned> &spf, bool remove_gap_ambig_sites)
  {
    std::map<faststring, unsigned> m;
    compute_site_pattern_map_no_clear(m, remove_gap_ambig_sites);

    std::map<faststring, unsigned>::iterator it, it_end;
    it     = m.begin();
    it_end = m.end();

    while (it != it_end)
    {
      spf.push_back(it->second);
      ++it;
    }
  }

  void clear(SequenceDataTypesEnum datatype_param = SeqType_unknown,
             char ambig_char_param = '?')
  {
    taxaNum    = 0;
    posNum     = 0;
    ambig_char = ambig_char_param;
    datatype   = datatype_param;

    originalPosNumbers.clear();
    originalPosNumbers_supplied = false;

    int i, n= (int)seqData.size();

    for (i=0; i<n; ++i)
      delete seqData[i];

    seqData.clear();
    sn_map.clear();
  }

  void remove_gap_only_sequences()
  {
    int i, n= (int)seqData.size();

    // We do this in reversed order, since this avoids adapting i and n after
    // an element has been removed.
    for (i=n-1; i>=0; --i)
    {
      if (seqData[i]->has_only_gaps() )
      {
        remove_seq_del(i);
      }
    }
  }

  void remove_ambig_only_sequences()
  {
    int i, n = (int)seqData.size();

    // We do this in reversed order, since this avoids adapting i and n after
    // an element has been removed.
    for (i=n-1; i>=0; --i)
    {
      if (seqData[i]->has_only_ambigs() )
      {
        remove_seq_del(i);
      }
    }
  }

  void remove_gap_or_ambig_only_sequences()
  {
    int i, n = (int)seqData.size();

    // We do this in reversed order, since this avoids adapting i and n after
    // an element has been removed.
    for (i=n-1; i>=0; --i)
    {
      //              std::cout << "Checking sequence: " << i+1 << '\n';
      if (seqData[i]->has_only_ambigs_or_gap() )
      {
        remove_seq_del(i);
      }
    }
  }

  void remove_gap_or_ambig_only_sites()
  {
    std::vector<bool> not_gap_ambig_only(posNum, false);

    faststring::size_type i, j;
    const char *the_seq;
    char c;

    if (datatype == SeqType_dna)
    {
      for (i=0; i< taxaNum; ++i)
      {
        the_seq = seqData[i]->getSeqStr();
        for(j=0; j<posNum; ++j)
        {
          c = the_seq[j];

          //     std::cout << "is_DNA_ambig(c):               " << (is_DNA_ambig(c)               ? "yes":"no") << '\n';

          if ( !(is_DNA_ambig(c) || c == '-') )
          {
            not_gap_ambig_only[j] = true;
          }
        }
      }
    }
    else if (datatype == SeqType_protein)
    {
      for (i=0; i< taxaNum; ++i)
      {
        the_seq = seqData[i]->getSeqStr();
        for(j=0; j<posNum; ++j)
        {
          c = the_seq[j];
          if ( !(is_aa_ambig(c) || c == '-') )
          {
            not_gap_ambig_only[j] = false;
          }
        }
      }
    }
    else
    {
      std::cerr << "Call to function remove_gap_or_ambig_only_sites() is only possible of the data"
      " type is either dna or protein.\n";
      exit(-33);
    }

    for (i=0; i< taxaNum; ++i)
    {
      seqData[i]->remove_with_vector(not_gap_ambig_only);
    }
  }


  void recode_with_replacement_string(faststring replace)
  {
    for (unsigned i=0; i < seqData.size(); ++i)
    {
      seqData[i]->recode_with_replacement_string(replace);
    }
  }

  void degap_all_sequences()
  {
    for (unsigned i=0; i < seqData.size(); ++i)
    {
      seqData[i]->degap();
    }
  }


  void steal_sequence(CSequences3 &seqs, unsigned i)
  {
    if (i < seqs.taxaNum)
    {
      //      ++taxaNum; // This is done in add_seq()
      add_seq(seqs.seqData[i]);
      seqs.remove_seq(i); // The sequence itself is not deleted, since it is stored in the *this element.
    }
  }


  void add_seq_to_dataset(SequenceDataTypesEnum datatype_param,
                          const faststring                  &fullname_param,
                          const faststring                  &seq_data_param,
                          char                              ambig_param='\0' // default auto
  )
  {
    if (datatype != datatype_param)
    {
      std::cerr << "Warning: Not all sequences seem to have the same data type.\n";
      datatype = SeqType_mixed;
    }

    //     std::cout << "datatype_param: " << datatype_param  << '\n';

    // Handle general ambig TODO XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    //     if (general_ambig == CSequence_Mol::unknown)

    CSequence_Mol   *seq;
    seq = new CSequence_Mol (datatype_param, ambig_param);
    //     std::cout << "Datatype after constructor: " << seq->get_datatype() << '\n';
    seq->set_taxon_and_sequence(fullname_param, seq_data_param, datatype_param);
    //     std::cout << "Datatype after constructor+set_taxon_and_sequence: " << seq->get_datatype() << '\n';

    add_seq(seq);

    // Further auto detect stuff???????
  }

  void get_vectors_of_site_coverages_and_site_entropies_DNA(std::vector<float> &coverage, std::vector<float> &entropies)
  {
    unsigned symbols[256];

    coverage.clear();
    entropies.clear();

    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    unsigned     site_index, seq_index;
    const char   **seq_strs;

    seq_strs = new const char* [taxaNum];

    for (seq_index=0; seq_index<taxaNum; ++seq_index)
    {
      seq_strs[seq_index] =  get_Seq_Data(seq_index);
    }

    for (site_index = 0; site_index < posNum; ++site_index)
    {
      std::memset(symbols, 0, 256*sizeof(unsigned));

      for (seq_index = 0; seq_index < taxaNum; ++seq_index)
      {
        unsigned char c = seq_strs[seq_index][site_index];
        c = toupper_lookup[c];
        ++symbols[c];
      }

      int A    = symbols[(unsigned char)'A'];
      int C    = symbols[(unsigned char)'C'];
      int G    = symbols[(unsigned char)'G'];
      int T    = symbols[(unsigned char)'T']+symbols[(unsigned char)'U'];
      //       int gaps = symbols[(int)'-'];
      int amb  =   symbols[(unsigned char)'R']+symbols[(unsigned char)'Y']+symbols[(unsigned char)'S']
      +symbols[(unsigned char)'W']+symbols[(unsigned char)'K']+symbols[(unsigned char)'M']
      +symbols[(unsigned char)'B']+symbols[(unsigned char)'D']+symbols[(unsigned char)'H']
      +symbols[(unsigned char)'V']+symbols[(unsigned char)'N']+symbols[(unsigned char)'?'];

      float all = A+C+G+T;
      float a   = (float)A/all;
      float c   = (float)C/all;
      float g   = (float)G/all;
      float t   = (float)T/all;

      float cov = (all+amb)/taxaNum;
      float ent = 0;
      if (a>0)
        ent -= a*log(a);
      if (c>0)
        ent -= c*log(c);
      if (g>0)
        ent -= g*log(g);
      if (t>0)
        ent -= t*log(t);

      coverage.push_back(cov);
      entropies.push_back(ent);
    }
    delete [] seq_strs;
  }


  void get_vectors_of_site_coverages_values_DNA(std::vector<unsigned> &coverage, bool count_ambig, unsigned start=0, faststring::size_type end=faststring::npos)
  {
    unsigned symbols[256];

    coverage.clear();
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    unsigned     site_index, seq_index;
    const char   **seq_strs;

    seq_strs = new const char* [taxaNum];

    for (seq_index=0; seq_index<taxaNum; ++seq_index)
    {
      seq_strs[seq_index] =  get_Seq_Data(seq_index);
    }

    if (end > posNum)
      end = posNum;

    for (site_index = start; site_index < end; ++site_index)
    {
      std::memset(symbols, 0, 256*sizeof(unsigned));

      for (seq_index = 0; seq_index < taxaNum; ++seq_index)
      {
        unsigned char c = seq_strs[seq_index][site_index];
        c = toupper_lookup[c];
        ++symbols[c];
      }

      unsigned A    = symbols[(unsigned char)'A'];
      unsigned C    = symbols[(unsigned char)'C'];
      unsigned G    = symbols[(unsigned char)'G'];
      unsigned T    = symbols[(unsigned char)'T']+symbols[(unsigned char)'U'];
      // unsigned gaps = symbols[(int)'-'];

      unsigned all = A+C+G+T;

      if (count_ambig)
      {
        unsigned amb  = symbols[(unsigned char)'R']+symbols[(unsigned char)'Y']+symbols[(unsigned char)'S']+symbols[(unsigned char)'W']
        +symbols[(unsigned char)'K']+symbols[(unsigned char)'M']+symbols[(unsigned char)'B']+symbols[(unsigned char)'D']
        +symbols[(unsigned char)'H']+symbols[(unsigned char)'V']+symbols[(unsigned char)'N']+symbols[(unsigned char)'?'];
        all += amb;
      }
      coverage.push_back(all);
    }
    delete [] seq_strs;
  }


  //  Returns a pointer to a 2D matrix with length of sequence number of columns and rows containing the number counts of: A, C, G, (T or U), (gap), ambig, '~', other i.e. 8 rows.
  //  The 2D array has a size of (posNum,8)
  unsigned **get_DNA_site_profile()
  {
    unsigned   **profile = new unsigned*[posNum];
    unsigned   symbols[256];
    const char **seq_strs;

    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return NULL;

    seq_strs = new const char* [taxaNum];

    for (unsigned seq_index=0; seq_index<taxaNum; ++seq_index)
    {
      seq_strs[seq_index] =  get_Seq_Data(seq_index);
    }

    for (unsigned site_index = 0; site_index < posNum; ++site_index)
    {
      std::memset(symbols, 0, 256*sizeof(unsigned));

      for (unsigned seq_index = 0; seq_index < taxaNum; ++seq_index)
      {
        unsigned char c = seq_strs[seq_index][site_index];
        c = toupper_lookup[c];
        ++symbols[c];
      }

      unsigned A     = symbols[(unsigned char)'A'];
      unsigned C     = symbols[(unsigned char)'C'];
      unsigned G     = symbols[(unsigned char)'G'];
      unsigned T     = symbols[(unsigned char)'T']+symbols[(unsigned char)'U'];
      unsigned amb   = symbols[(unsigned char)'R']+symbols[(unsigned char)'Y']+symbols[(unsigned char)'S']+symbols[(unsigned char)'W']
      +symbols[(unsigned char)'K']+symbols[(unsigned char)'M']+symbols[(unsigned char)'B']+symbols[(unsigned char)'D']
      +symbols[(unsigned char)'H']+symbols[(unsigned char)'V']+symbols[(unsigned char)'N']+symbols[(unsigned char)'?'];
      unsigned gaps  = symbols[(int)'-'];
      unsigned tilde = symbols[(int)'~'];
      //      unsigned all_DNA = A+C+G+T;

      profile[site_index] = new unsigned [7];

      profile[site_index][0] = A;
      profile[site_index][1] = C;
      profile[site_index][2] = G;
      profile[site_index][3] = T;
      profile[site_index][4] = gaps;
      profile[site_index][5] = amb;
      profile[site_index][6] = tilde;
    }
    delete [] seq_strs;
    return profile;
  } // END unsigned **get_DNA_site_profile()



  // Gets the data from *this object and stores the new sequences in seqs.
  void alignment_without_gap_ambig_only_positions(CSequences3 &seqs)
  {
    seqs.clear(datatype, seqs.ambig_char);

    // Action does not make sense in these cases:
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    std::vector<bool> gap_ambig_positions(posNum, true);

    faststring::size_type i, j;
    const char *the_seq;
    char c;

    if (datatype == SeqType_dna)
    {
      for (i=0; i< taxaNum; ++i)
      {
        the_seq = seqData[i]->getSeqStr();
        for(j=0; j<posNum; ++j)
        {
          c = the_seq[j];

          //	   std::cout << "is_DNA_ambig(c):               " << (is_DNA_ambig(c)               ? "yes":"no") << '\n';

          if ( !(is_DNA_ambig(c) || c == '-') )
          {
            gap_ambig_positions[j] = false;
          }
        }
      }
    }
    else if (datatype == SeqType_protein)
    {
      for (i=0; i< taxaNum; ++i)
      {
        the_seq = seqData[i]->getSeqStr();
        for(j=0; j<posNum; ++j)
        {
          c = the_seq[j];
          if ( !(is_aa_ambig(c) || c == '-') )
          {
            gap_ambig_positions[j] = false;
          }
        }
      }
    }
    else
    {
      std::cerr << "Call to function alignment_without_gap_ambig_only_positions() is only possible of the data"
      " type is either dna or protein.\n";
      exit(-33);
    }

    faststring new_seq;

    for (i=0; i< taxaNum; ++i)
    {
      new_seq.clear();
      the_seq = seqData[i]->getSeqStr();

      for(j=0; j<posNum; ++j)
      {
        if (!gap_ambig_positions[j])
          new_seq.push_back(the_seq[j]);
      }
      seqs.add_seq_to_dataset(datatype, seqData[i]->getFullName(), new_seq, seqData[i]->get_ambiguity_character() );
    }
  }


  void get_vector_gap_ambig_only_positions(std::vector<bool> &gap_ambig_positions)
  {
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    gap_ambig_positions.assign(posNum, true);

    faststring::size_type i, j;
    const char *the_seq;
    char c;

    if (datatype == SeqType_dna)
    {
      for (i=0; i< taxaNum; ++i)
      {
        the_seq = seqData[i]->getSeqStr();
        for(j=0; j<posNum; ++j)
        {
          c = the_seq[j];

          //	   std::cout << "is_DNA_ambig(c):               " << (is_DNA_ambig(c)               ? "yes":"no") << '\n';

          if ( !(is_DNA_ambig(c) || c == '-') )
          {
            gap_ambig_positions[j] = false;
          }
        }
      }
    }
     else if (datatype == SeqType_protein)
    {
      for (i=0; i< taxaNum; ++i)
      {
        the_seq = seqData[i]->getSeqStr();
        for(j=0; j<posNum; ++j)
        {
          c = the_seq[j];
          if ( !(is_aa_ambig(c) || c == '-') )
          {
            gap_ambig_positions[j] = false;
          }
        }
      }
    }
    else
    {
      std::cerr << "Call to function alignment_without_gap_ambig_only_positions() is only possible of the data type is either dna or protein.\n";
      exit(-33);
    }
  }

  // Gets the data from *this object and stores the new sequences in seqs.
  void alignment_without_gap_ambig_lowerCase_only_positions(CSequences3 &seqs)
  {
    seqs.clear(datatype, seqs.ambig_char);
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    std::vector<bool> gap_ambig_LowerCase_positions(posNum, true);

    faststring::size_type i, j;
    const char *the_seq;
    char c;

    if (datatype == SeqType_dna)
    {
      // TODO: Change loops, use break if one ... has been found.
      //       std::cout << "DNA\n";
      for (i=0; i< taxaNum; ++i)
      {
        the_seq = seqData[i]->getSeqStr();
        for(j=0; j<posNum; ++j)
        {
          c = the_seq[j];
          if (!(is_DNA_ambig(c) || c == '-' || (c >= 'a' && c <= 'z')) )
          {
            gap_ambig_LowerCase_positions[j] = false;
          }
        }
      }
    }
    else if (datatype == SeqType_protein)
    {
      //       std::cout << "PROTEIN\n";
      for (i=0; i< taxaNum; ++i)
      {
        the_seq = seqData[i]->getSeqStr();
        for(j=0; j<posNum; ++j)
        {
          c = the_seq[j];
          if ( !(is_aa_ambig(c) || c == '-' || (c >= 'a' && c <= 'z')) )
          {
            gap_ambig_LowerCase_positions[j] = false;

            //	     std::cout << "FOUND: " << j+1 << " " << i+1 << " " << c << " "<< gap_ambig_LowerCase_positions[j] << '\n';
          }
        }
      }
    }
    else
    {
      std::cerr << "Call to function alignment_without_gap_ambig_lowerCase_positions(...) is only possible of the data type is either dna or protein.\n";
      exit(-33);
    }

    faststring new_seq;

    for (i=0; i< taxaNum; ++i)
    {
      new_seq.clear();
      the_seq = seqData[i]->getSeqStr();

      for(j=0; j<posNum; ++j)
      {
        //	 std::cout << "Position: " << j <<  " " << gap_ambig_LowerCase_positions[j] << '\n';
        if (!gap_ambig_LowerCase_positions[j])
          new_seq.push_back(the_seq[j]);
      }
      seqs.add_seq_to_dataset(datatype, seqData[i]->getFullName(), new_seq, seqData[i]->get_ambiguity_character() );
    }
  }


  // Gets the data from *this object and stores the new sequences in seqs.
  void get_vector_gap_ambig_lowerCase_only_positions(std::vector<bool> &gap_ambig_lowercase_positions)
  {
    gap_ambig_lowercase_positions.assign(posNum, true);
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    faststring::size_type i, j;
    const char *the_seq;
    char c;

    if (datatype == SeqType_dna)
    {
      //       std::cout << "DNA" << '\n';
      // TODO: Change loops, use break if one ... has been found.
      //       std::cout << "DNA" << '\n';
      for (i=0; i< taxaNum; ++i)
      {
        the_seq = seqData[i]->getSeqStr();
        for(j=0; j<posNum; ++j)
        {
          c = the_seq[j];
          if (!(is_DNA_ambig(c) || c == '-' || (c >= 'a' && c <= 'z')) )
          {
            gap_ambig_lowercase_positions[j] = false;
          }
        }
      }
    }
    else if (datatype == SeqType_protein)
    {
      //       std::cout << "PROTEIN" << '\n';
      for (i=0; i< taxaNum; ++i)
      {
        the_seq = seqData[i]->getSeqStr();
        for(j=0; j<posNum; ++j)
        {
          c = the_seq[j];
          if ( !(is_aa_ambig(c) || c == '-' || (c >= 'a' && c <= 'z')) )
          {
            gap_ambig_lowercase_positions[j] = false;

            //	     std::cout << "FOUND: " << j+1 << " " << i+1 << " " << c << " " << gap_ambig_LowerCase_positions[j] << '\n';
          }
        }
      }
    }
    else
    {
      std::cerr << "Call to function alignment_without_gap_ambig_lowerCase_positions(...) is only possible of the data type is either dna or protein.\n";
      exit(-33);
    }
  }

  // Gets the data from *this object and stores the new sequences in seqs.
  void alignment_without_gap_only_positions(CSequences3 &seqs)
  {
    //     std::cout << "Entering alignment_without_gap_only_positions" << '\n';

    seqs.clear(datatype, ambig_char);
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    // Count the number of gaps at all positions.
    // Gap only positions are those that have taxaNum gaps.
    std::vector<bool> gap_positions(posNum, true);

    faststring::size_type i, j;
    const char *the_seq;
    char c;

    for (i=0; i< taxaNum; ++i)
    {
      the_seq = seqData[i]->getSeqStr();
      for(j=0; j<posNum; ++j)
      {
        c = the_seq[j];
        if (c != '-')
        {
          gap_positions[j] = false;
        }
      }
    }

    faststring new_seq_str;

    for (i=0; i< taxaNum; ++i)
    {
      new_seq_str.clear();
      the_seq = seqData[i]->getSeqStr();

      for(j=0; j<posNum; ++j)
      {
        if (!gap_positions[j]) // not only gaps
          new_seq_str.push_back(the_seq[j]);
      }
      seqs.add_seq_to_dataset(datatype,seqData[i]->getFullName(),new_seq_str,
                                seqData[i]->get_ambiguity_character());
    }
    /*      std::cout << "DEBUGOUT-alignment_without_gap_only_positions" << '\n'; */
    /*      seqs.print_DEBUG(std::cout,3); */
  }

  // Gets the data from *this object and stores the new sequences in seqs.
  void get_vector_gap_only_positions(std::vector<bool> &gap_positions)
  {
    // Count the number of gaps at all positions.
    // Gap only positions are those that have taxaNum gaps.
    gap_positions.assign(posNum, true);
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    const char *the_seq;
    char c;

    // TODO: Swap loops and use break to make this more efficient.
    for (faststring::size_type i=0; i< taxaNum; ++i)
    {
      the_seq = seqData[i]->getSeqStr();
      for(faststring::size_type j=0; j<posNum; ++j)
      {
        c = the_seq[j];
        if (c != '-')
        {
          gap_positions[j] = false;
        }
      }
    }
  }


  // Copy data from this to seqs, but without the positions specified in the vector
  void alignment_without_specified_positions(CSequences3 &seqs, std::vector<bool> without)
  {
    const char *the_seq;
    seqs.clear(datatype, seqs.ambig_char);
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    faststring new_seq;
    for (faststring::size_type i=0; i< taxaNum; ++i)
    {
      new_seq.clear();
      the_seq = seqData[i]->getSeqStr();

      for(faststring::size_type j=0; j<posNum; ++j)
      {
        //	 std::cout << "Position: " << j <<  " " << gap_ambig_LowerCase_positions[j] << '\n';
        if (!without[j])
          new_seq.push_back(the_seq[j]);
      }
      seqs.add_seq_to_dataset(datatype, seqData[i]->getFullName(), new_seq, seqData[i]->get_ambiguity_character() );
    }
  }


  // Copy data from this to seqs, but without the positions specified in the vector
  // assumeing the specified positions come from an amino acid alignment. Positions are removed
  // in a nucleotid alignment.
  void alignment_without_specified_positions_prot2nuc(CSequences3 &seqs, std::vector<bool> without)
  {
    faststring::size_type j, k, h;
    const char *the_seq;
    seqs.clear(datatype, seqs.ambig_char);
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
      return;

    faststring new_seq;
    for (faststring::size_type i=0; i< taxaNum; ++i)
    {
      new_seq.clear();
      the_seq = seqData[i]->getSeqStr();

      for (j=0, k=0, h=0; k<posNum; ++k, ++h)
      {
        // j is the aa index, i.e. the index in the without array
        // k is the seq index
        // h is a helper index
        if (h==3)
        {
          h = 0;
          ++j;
        }
        //	 std::cout << "Position: " << j <<  " " << gap_ambig_LowerCase_positions[j] << '\n';
        if (!without[j])
          new_seq.push_back(the_seq[k]);
      }
      seqs.add_seq_to_dataset(datatype, seqData[i]->getFullName(), new_seq, seqData[i]->get_ambiguity_character() );
    }
  }


  void replace_sequence_interval(const CSequences3 &s, faststring::size_type pos1=0, faststring::size_type pos2=faststring::npos)
  {
    determine_map_of_sequence_names();

    faststring name;
    faststring::size_type i, N = s.seqData.size();
    std::map<faststring, CSequence_Mol*>::iterator find_it;

    if (taxaNum != s.taxaNum)
    {
      std::cerr << "Critical error in replace_sequence_interval: Number of taxa are not equal\n";
      exit(-1);
    }

    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
    {
      std::cerr << "ERROR in replace_sequence_interval. Sequences have unequal or 0 lenths\n";
      exit(-24);
    }

    for (i=0; i<N; ++i)
    {
      name = s.seqData[i]->getName();
      find_it = sn_map.find(name);
      if (find_it == sn_map.end())
      {
        std::cerr << "Critical error in replace_sequence_interval: Sequence "
        << name << " does not exist.\n";
        exit(-1);
      }
      find_it->second->replace_part_of_sequence(*(s.seqData[i]), pos1, pos2);
    }
    if (taxaNum > 0 && equal_length_of_all_sequences())
      posNum = seqData[0]->length();
    else
      posNum = SIZE_MAX;
  }

  void append_to_sequences(const CSequences3 &s)
  {
    determine_map_of_sequence_names();

    faststring name;
    faststring::size_type i,N=s.seqData.size();
    std::map<faststring, CSequence_Mol*>::iterator find_it;

    if (taxaNum != s.taxaNum)
    {
      std::cerr << "Critical error in append_sequences: Number of taxa are not equal. Exiting.\n";
      exit(-1);
    }

    for (i=0; i<N; ++i)
    {
      name = s.seqData[i]->getName();
      find_it = sn_map.find(name);
      if (find_it == sn_map.end())
      {
        std::cerr << "Critical error in append_sequences: Sequence "
        << name << " does not exist.\n";
        exit(-1);
      }
      find_it->second->append_sequence(*(s.seqData[i]));
    }
    if (taxaNum > 0 && equal_length_of_all_sequences())
      posNum = seqData[0]->length();
    else
      posNum = 0;
  }

  // Append a column to an alignment.
  void append_to_sequences(const faststring &site_pattern)
  {
    unsigned i;

    if ( taxaNum != site_pattern.size() )
    {
      std::cerr << "Error when calling the append_to_sequences member function:\n"
      "the number of sites in the pattern that shall be\n"
      "appended the sequences object differs from the\n"
      "number of sequences. This pattern will not be appended.\n";
      exit(-22);
    }
    for (i=0; i<taxaNum; ++i)
    {
      seqData[i]->append_residue_unchecked(site_pattern[i]);
    }
    if (taxaNum > 0 && equal_length_of_all_sequences())
      posNum = seqData[0]->length();
    else
      posNum = 0;
  }


  bool is_insertion_column(unsigned col_pos)
  {
    unsigned i_tax;
    char c;

    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
    {
      std::cerr << "ERROR in element function is_insertion_column. Sequences have unequal or 0 lenths\n";
      exit(-24);
    }

    for (i_tax=0; i_tax<taxaNum; ++i_tax)
    {
      c = seqData[i_tax]->get_pos(col_pos);
      if (c == '.' || (c >='a' && c <= 'z') )
        return true;
    }
    return false;
  }


  bool is_mainly_gap_pos(unsigned col_pos, double min_prop)
  {
    unsigned i_tax;
    char c;
    unsigned count_gap= 0;

    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
    {
      std::cerr << "ERROR in element function is_mainly_gap_pos. Sequences have unequal or 0 lenths\n";
      exit(-24);
    }

    for (i_tax=0; i_tax<taxaNum; ++i_tax)
    {
      c = seqData[i_tax]->get_pos(col_pos);
      if (c == '-' || c == '.')
        ++count_gap;
    }
    if ((float)count_gap/taxaNum >= min_prop)
      return true;
    else
      return false;
  }


  // Assumes that insertion positions have lower case characters
  // and they could have dots '.' instead of gaps.
  // The letter is not necessary, since insertion positions should
  // at least have one lower case character.
  // So dots could have been replaced by gaps already.
  unsigned get_next_insertion_column(unsigned col_pos)
  {
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
    {
      std::cerr << "ERROR in element function get_next_insertion_column. Sequences have unequal or 0 lenths\n";
      exit(-24);
    }

    for (; col_pos < posNum; ++col_pos)
      if (is_insertion_column(col_pos))
        return col_pos;
    return col_pos;
  }

  unsigned get_next_insertion_or_mainly_gap_column(unsigned col_pos, float min_gap_prop)
  {
    for (; col_pos < posNum; ++col_pos)
      if (is_insertion_column(col_pos) || is_mainly_gap_pos(col_pos, min_gap_prop))
        return col_pos;
    return col_pos;
  }

  // Assumes that non insertion positions have upper case characters.
  // If dots have already been converted to gaps, they are not
  // indicative of non insertion positions!!
  unsigned get_next_non_insertion_column(unsigned col_pos)
  {
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
    {
      std::cerr << "ERROR in element function get_next_non_insertion_column. Sequences have unequal or 0 lenths\n";
      exit(-24);
    }

    for (; col_pos < posNum; ++col_pos)
      if (!is_insertion_column(col_pos))
        return col_pos;
    //	 if (c == '.' || c >='a' && c <= 'z') // Not a non-insertion column
    return col_pos;
  }

  unsigned get_next_non_insertion_and_non_mainly_gap_column(unsigned col_pos, float  min_gap_prop)
  {
    if (taxaNum == 0 || posNum == 0 || posNum == SIZE_MAX)
    {
      std::cerr << "ERROR in get_next_non_insertion_and_non_mainly_gap_column. Sequences have unequal or 0 lenths\n";
      exit(-24);
    }

    for (; col_pos < posNum; ++col_pos)
      if (!(is_insertion_column(col_pos) || is_mainly_gap_pos(col_pos, min_gap_prop)) )
        return col_pos;
    //	 if (c == '.' || c >='a' && c <= 'z') // Not a non-insertion column
    return col_pos;
  }

  // RETURN VALUE NEEDS TO BE SET. DO NOT USE IN THIS FORM.
  /*    void get_bp_entropies_nucleotide(double *ent) */
  /*    { */
  /*      const char **tax = new const char* [taxaNum]; */
  /*      unsigned  j; */

  /*      if (taxaNum == 0 || posNum == 0) */
  /*        return; */

  /*      // Store pointers to all sequences in the tax array - this will be faster */
  /*      for (j=0; j < taxaNum; ++j) */
  /*      { */
  /*        tax[j] = seqData[j]->getSeqStr(); */
  /*      } */
  /*    } */
  // Computes relative amount of residues that are non gap and non ambig in this range:
  /*
  double get_non_ambig_non_gap_coverage_in_range(unsigned pos, unsigned pos_end, int flag=0)
  {
    const char **tax = new const char* [taxaNum];
    unsigned i, j;
    unsigned N = 0;


    if (taxaNum == 0 || posNum == 0)
      return -1; // Indicates an error

    // Store pointers to all sequences in the tax array - this will be faster
    if (flag == 0 || flag == 1) // 1kite all or general mode
    {
      N = taxaNum;
      for (j=0; j < taxaNum; ++j)
      {
        tax[j] = seqData[j]->getSeqStr();
      }
    }
    else if (flag == 13) // Special 1kite flag: Filter hamstered reference taxa:
    {
#ifdef MY_PROJECT_DEBUG
      std::cerr << "Called: CSequences3.h:get_non_ambig_non_gap_coverage_in_range with flag==1\n";
#endif

      unsigned jj=0;

      for (j=0; j < taxaNum; ++j)
      {
        // In case we compute coverage values, we only have hamstered reference taxa left in 1kite:
        if (is_1kite_hamstered_reference_seq_id(seqData[j]->getFullName() ) )
        {
          tax[jj] = seqData[j]->getSeqStr();
          ++jj;
        }
      }
      N = jj;
    }
    else
    {
      return -1; // Indicates an error
    }

    unsigned count_all=0;
    unsigned count_non_ambig_non_gap=0;
    double   res;
    char     c;

    if (datatype ==  CSequence_Mol::protein)
    {
      for (i =pos; i< pos_end && i < posNum; ++i)
      {
        for (j=0; j< N; ++j)
        {
          c = tax[j][i];
          ++count_all;
          if ( !(is_aa_ambig(c) || c == '-') )
            ++count_non_ambig_non_gap;
        }
      }
    }
    else if (datatype ==  CSequence_Mol::dna)
    {
      for (i =pos; i< pos_end && i < posNum; ++i)
      {
        for (j=0; j< taxaNum; ++j)
        {
          c = tax[j][i];
          ++count_all;
          if ( !(is_DNA_ambig(c) || c == '-' ) )
            ++count_non_ambig_non_gap;
        }
      }
    }
    res = count_non_ambig_non_gap;
    res /= count_all;

    //          std::cerr << count_non_ambig_non_gap << '\n';
    //          std::cerr << count_all               << '\n';

    return res;
  }
  */

/*
  // Computes relative amount of gaps in sequence
  double get_gap_proportion(unsigned seq_index)
  {
    return  seqData[seq_index]->get_gap_proportion();
  }
*/

  /*    int move_1kite_core_of_other_seqs_to_top(const CSequences3 &other_seqs) */
  /*    { */
  /* #ifdef MY_PROJECT_DEBUG */
  /*      std::cerr << "Called: CSequences3.h:move_1kite_core_of_other_seqs_to_top" << '\n'; */
  /* #endif */

  /*      std::vector<faststring> vec; */

  /*      other_seqs.get_1kite_cores_vec(vec); */

  /*      int i; */

  /*      for (i=vec.size()-1; i>= 0; --i) */
  /*      { */
  /*        reorder_move_seq_to_top(vec[i]); */
  /*      } */
  /*      return vec.size(); */
  /*    } */



  // Currently only implemented for datatype DNA
  //  Allowed flag values:
  // 0:  normal mode
  // 1:  Gap regions always get a distance of -3
  // 2:  Compute gap distance, i.e. count gap versus non-gap nucleotides.
  // The result is in the range from 0..1.
  // Ns are treated as mismatch. If this is not desired, use the other overwrite of this function below.


  double sequence_distance(unsigned s1, unsigned s2, size_t begin_range, size_t end_range, unsigned flag=0)
  {
    if (s1 >= taxaNum || s2 >= taxaNum)
      return -1;

    CSequence_Mol* seq1 = seqData[s1];
    CSequence_Mol* seq2 = seqData[s2];

    return sequence_distance(seq1, seq2, begin_range, end_range, flag=0);
  }


  double sequence_distance(CSequence_Mol* seq1, CSequence_Mol* seq2, faststring::size_type begin_range, faststring::size_type end_range, unsigned flag=0)
  {
    char   c1, c2;

    // Do some range checks:
    if (end_range < begin_range) // This is considered a major error and the caller should correct this problem.
    {
      std::cerr << "Error in sequence_distance: end_range < begin_range\n";
      exit(-33);
    }

    // Out of range checks are corrected silently.
    // It is not clear this is always the best method to deal with this.
    if ((end_range > seq1->length()))
      end_range = seq1->length();
    if ((end_range > seq2->length()))
      end_range = seq2->length();
    if (begin_range > end_range)
      begin_range = end_range;

    faststring::size_type i = begin_range;

    unsigned len  = 0;
    unsigned diff = 0;

    if (flag < 2)
    {
      if (seq1->get_datatype() != SeqType_dna)
        return -2;

      while (i != end_range)
      {
        // TODO: DO WE NEED TOUPPER??
        c1 = toupper(seq1->get_pos(i));
        c2 = toupper(seq2->get_pos(i));

        if (c1 != '-' && c2 != '-')
        {
          ++len;
          if (c1 != c2 || c1 == 'N') // If they differ or if both are Ns we count this as a difference.
          {
            ++diff;
          }
        }
        else
        {
          if (flag == 1)
            return -3.0;
        }
        ++i;
      }
      if (len == 0)
      {
        diff = 2;
        len  = 1;
      }
      return (double)diff / (double)len;
    }
    else if (flag == 2)
    {
      //       len = end_range - begin_range;

      while (i != end_range)
      {
        // TODO: DO WE NEED TOUPPER??
        c1 = toupper(seq1->get_pos(i));
        c2 = toupper(seq2->get_pos(i));

        // We only count differences in the gap structur: Non-gap versus gap.
        if ( (c1 == '-' && c2 != '-') || (c1 != '-' && c2 == '-') )
        {
          ++diff;
        }
        ++i;
      }
      return diff;
    }
    else
    {
      std::cerr << "Flag value undefined in CSequences3:sequence_distance(..);\n";
      exit(-127);
    }
  }

  // Counts overlap as the number of positions in which both sequences are unequal to gaps.
  void sequences_overlap_mismatches(CSequence_Mol* seq1, CSequence_Mol* seq2, unsigned begin_range, unsigned end_range, unsigned &overlap, unsigned &mismatches, unsigned &DNA_range_start, unsigned &DNA_range_end)
  {
    const char  *s1 = seq1->getSeqStr();
    const char  *s2 = seq2->getSeqStr();

    char   c1, c2;

    overlap    = 0;
    mismatches = 0;
    DNA_range_start = 0;
    DNA_range_end   = UINT_MAX;

    // Do some range checks:
    if (end_range < begin_range) // This is considered a major error and the caller should correct this problem.
      return;

    // Out of range problems are corrected silently.
    // It is not clear this is always the best method to deal with this.
    if ((end_range > seq1->length()))
      end_range = (unsigned) seq1->length();
    if ((end_range > seq2->length()))
      end_range = (unsigned) seq2->length();
    if (begin_range > end_range)
      begin_range = end_range;

    // Determine effective begin range:
    unsigned i = begin_range;
    c1 = s1[i];
    c2 = s2[i];
    while (i < end_range && c1 == '-' && c2 == '-')
    {
      ++i;
      c1 = s1[i];
      c2 = s2[i];
    }
    begin_range = i;

    // Determine effective end range:
    i = end_range - 1;
    c1 = s1[i];
    c2 = s2[i];
    while (begin_range < i && c1 == '-' && c2 == '-')
    {
      --i;
      c1 = s1[i];
      c2 = s2[i];
    }
    end_range = i+1;

    unsigned local_overlap    = 0;
    unsigned local_mismatches = 0;

    i = begin_range;
    while (i != end_range)
    {
      c1 = s1[i];
      c2 = s2[i];
      if (c1 != '-' && c2 != '-')
      {
        ++local_overlap;
        if (c1 != c2 || is_DNA_ambig(c1) || is_DNA_ambig(c2) )
          ++local_mismatches;
      }
      ++i;
    }
    mismatches = local_mismatches;
    overlap    = local_overlap;
    DNA_range_start = begin_range;
    DNA_range_end   = end_range;
  }



  // Currently only implemented for datatype DNA
  // Note: N positions are not counted if the other sequence contains a gap.
  // In this case only the gap and not the N is counted, since the sequence positions
  // are not comparable.
  // The function does to compute the distance but the length, differences, gaps
  // and Ns positions.
  // Note: length_comparable_nucleotieds + gaps_extensions + Ns must add up to
  //       the length of the range.
  void sequence_distance(unsigned s1, unsigned s2, size_t begin_range, size_t end_range,
                         unsigned &differences,  unsigned &length_comparable_nucleotieds,
                         unsigned &gap_openings, unsigned &gaps_extensions, unsigned &Ns)
  {
    size_t i = begin_range;
    char     c1, c2;

    if (s1 >= taxaNum || s2 >= taxaNum)
    {
      differences = UINT_MAX;
      return;
    }

    CSequence_Mol* seq1 = seqData[s1];
    CSequence_Mol* seq2 = seqData[s2];

    if (end_range < begin_range)
    {
      std::cerr << "Error in sequence_distance: end_range < begin_range\n";
      exit(-33);
    }

    if ((end_range - begin_range > seq1->length()) || (end_range - begin_range > seq2->length()) )
    {
      end_range = begin_range + seq1->length();
      if (end_range - begin_range > seq2->length())
        end_range = begin_range + seq2->length();
    }

    unsigned local_len           = 0;
    unsigned local_numgapsopen   = 0;
    unsigned local_numgapext     = 0;
    unsigned local_differences   = 0;
    unsigned local_numNs         = 0;
    bool     in_gap        = false;

    if (seq1->get_datatype() != SeqType_dna)
    {
      differences = -2;
      return;
    }


    while (i != end_range)
    {
      // TODO: DO WE NEED TOUPPER??
      c1 = toupper(seq1->get_pos(i));
      c2 = toupper(seq2->get_pos(i));

      if (c1 == '-' || c2 == '-')
      {
        if (!in_gap)
        {
          ++local_numgapsopen;
          in_gap = true;
        }
        ++local_numgapext;
      }
      else
      {
        in_gap = false;
        if (c1 == 'N' || c2 == 'N')
          ++local_numNs;
        else
        {
          ++local_len;
          if (c1 != c2)
            ++local_differences;
        }
      }
      ++i;
    }

    differences                   = local_differences;
    length_comparable_nucleotieds = local_len;
    gap_openings                  = local_numgapsopen;
    gaps_extensions               = local_numgapext;
    Ns                            = local_numNs;
  }


  // Currently only implemented for datatype DNA
  //  Allowed flag values:
  // 0: normal mode.
  // 1: Gap regions always get a distance of -3.
  void get_sequences_distances(size_t begin_range, size_t end_range, void (*call_back_distance)(short, short, double), unsigned flag=0)
  {
    // Compute all pairwise distances:
    unsigned N = taxaNum;
    double dist;

    for (unsigned i=0; i<N; ++i)
    {
      for (unsigned j=i+1; j<N; ++j)
      {
        dist = sequence_distance(i, j, begin_range, end_range, flag);
        call_back_distance(i,j,dist);
      }
    }
  }

  // Same as above, but for some distances, those that have a false value in vector
  // we assign an invalidly large value to the distance.
  // Here a value of 100 is assigned to those distances
  void get_sequences_distances(faststring::size_type begin_range, faststring::size_type end_range, void (*call_back_distance)(short, short, double), const std::vector<bool> &sequences_valid, unsigned flag=0)
  {
    // Compute all pairwise distances:
    unsigned N = taxaNum;
    double dist;

    if (sequences_valid.size() != N)
    {
      std::cerr << "Critical internal error: In function get_sequences_distances the size of the vector sequences_valid must be equal to the number of sequences in the sequences object.\n";
      exit(-3);
    }

    for (unsigned i=0; i<N; ++i)
    {
      for (unsigned j=i+1; j<N; ++j)
      {
        if (!sequences_valid[i] || !sequences_valid[j])
        {
          dist = 100;
        }
        else
        {
          dist = sequence_distance(i, j, begin_range, end_range, flag);
        }
        call_back_distance(i,j,dist);
      }
    }
  }


  CSequence_Mol* get_sequence_with_sequence_name_match(faststring &partial_name)
  {
    for (unsigned i=0; i < taxaNum; ++i)
    {
      if ( seqData[i]->getFullName_faststring().find(partial_name) != faststring::npos )
      {
        return seqData[i];
      }
    }
    return NULL;
  }

  void get_alignment_as_copy(std::vector<faststring> &v)
  {
    v.clear();
    v.reserve(taxaNum+1);

    for (unsigned i=0; i<taxaNum; ++i)
    {
      v.push_back(seqData[i]->getSeq_faststring());
      v[i].c_str();  // Make buffer 0-terminated.
    }
  }

};



// inline CSequences::CSequences():posNum(0),datatype(0),splits(NULL){}


//inline string    CSequences::GetTaxonLabel(unsigned i)
//                             {return taxa.GetTaxonLabel(i);}
//inline void      CSequences::SetTaxonLabel(unsigned i, string s)
//                             {taxa.SetTaxonLabel(i,s);}






/* inline void PrintSpecial_TaxaSeq(ostream& os, */
/* 			  const std::vector<unsigned>& taxa_vec, */
/* 			  const std::vector<unsigned>& pos_vec, */
/* 			  const CTaxa       *ptaxa, */
/* 			  const CSequences3 *pseq, */
/* 			  const char*      symbolsList, */
/* 			  bool       printTaxonLabels, */
/* 			  unsigned   taxonLabelPrintWidth, */
/* 			  unsigned   maxTaxonLabelLength) */
/* { */
/*   std::vector<unsigned>::const_iterator taxa_vec_it; */
/*   std::vector<unsigned>::const_iterator taxa_vec_end = taxa_vec.end(); */

/*   std::vector<unsigned>::const_iterator pos_vec_it; */
/*   std::vector<unsigned>::const_iterator pos_vec_end = pos_vec.end(); */

/*   unsigned theTaxonNumber; */

/*   std::ios_base::fmtflags orig_ios_flag = os.flags(); */
/*   os.flags(orig_ios_flag | ios::left); */

/*   for (taxa_vec_it=taxa_vec.begin(); taxa_vec_it!=taxa_vec_end; ++taxa_vec_it) */
/*   { */
/*     theTaxonNumber = *taxa_vec_it; */

/*     if (printTaxonLabels) */
/*     { */
/*       os.width(maxTaxonLabelLength); */
/*       faststring tl = ptaxa->GetTaxonLabel(theTaxonNumber); */
/*       tl.shorten(maxTaxonLabelLength); */
/*       os << tl; */
/*       int morespaces = (int)taxonLabelPrintWidth-(int)maxTaxonLabelLength; */

/*       if ( morespaces > 0 ) */
/* 	{ */
/* 	  os.width(morespaces); */
/* 	  os << " "; */
/* 	} */
/*     } */

/*     for (pos_vec_it = pos_vec.begin(); pos_vec_it != pos_vec_end; ++pos_vec_it) */
/*     { */
/*       os << symbolsList[pseq->GetChar(theTaxonNumber,*pos_vec_it)]; */
/*     } */
/*     os << '\n'; */
/*   } */
/*   os.flags(orig_ios_flag); */
/* } */

/* inline  */
/* void PrintSpecial_PosNumbers(ostream& os, */
/* 			     std::vector<unsigned> pos_vec, // we get a real copy */
/* 			     const CSequences3 *pseq, */
/* 			     unsigned         taxonLabelPrintWidth, */
/* 			     bool             UseNexusComments) */
/* { */
/*   unsigned i; */
/*   unsigned numPos           = pos_vec.size(); */
/*   unsigned printNoZerosTill = numPos; */
/*   bool     printNoZero_set  = false; */

/*   // We determine the original position numbers: */
/*   for (i=0; i<numPos; ++i) */
/*     pos_vec[i] = pseq->GetOriginalPosNumber(pos_vec[i]); */

/*   unsigned maxDigits = ((unsigned)log10((double)pos_vec[numPos-1]))+1; */
/*   unsigned factor    = 1; */
/*   unsigned digit; */

/*   for (i=1; i<maxDigits; ++i) */
/*     factor *= 10; */

/*   for (; factor; factor /= 10, printNoZero_set = false) */
/*   { */
/*     if (UseNexusComments) */
/*     { */
/*       os.width(taxonLabelPrintWidth); */
/*       os << "[ "; */
/*     } */
/*     else */
/*     { */
/*       os.width(taxonLabelPrintWidth); */
/*       os << ' '; */
/*     } */
/*     for (i=0; i < numPos; ++i) */
/*     { */
/*       digit = pos_vec[i]/factor; */

/*       if ( i >= printNoZerosTill || digit > 0) */
/*       { */
/* 	if (!printNoZero_set) */
/* 	{ */
/* 	  printNoZerosTill = i; */
/* 	  printNoZero_set = true; */
/* 	} */
/* 	os << (char)(digit+'0'); */
/* 	pos_vec[i] -= factor*digit; */
/*       } */
/*       else */
/*       { */
/* 	os << ' '; */
/*       } */
/*     } */
/*     if (UseNexusComments) */
/*     { */
/*       os << "]" << '\n'; */
/*     } */
/*     else */
/*     { */
/*       os << '\n'; */
/*     } */
/*   } */
/* } */


/* inline void      CSequences::GetBaseFrequencies(unsigned ti, */
/* 						std::vector<unsigned> bf) */
/* { */
/*   unsigned pos; */
/*   for (pos = 0; pos < posNum; ++pos) */
/*     ++bf[seqData[ti][pos]] */
/* } */


#endif
