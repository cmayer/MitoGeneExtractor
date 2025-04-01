#ifndef CSEQUENCE_MOL3_H
#define CSEQUENCE_MOL3_H

// Changes:
// 22.08.2009: According to the fasta format specifications, the sequence name
//             is only the string between the ">" and the first space. The rest
//             of the line contains a sequence description or annotation
//             information.
//             This is now reflected in the CSequence_Mol class. The variables
//             full_name, name, and description contain this information.
//             The sequence name variables are now of type faststring (not string
//             as before).
// 07.01.2011: Migrated to fastring2.h and CFile2_1.h
// 30.04.2011: Sliding window CG-AT analysis
// 25.12.2020: Changed comment symbol for fasta files from '#' to ';', which is the standard comment symbol.
// 25.12.2020: Unset fail flag in public functions reading data, if >0 residues
//             have been read. Could break old code!!!!! Double check if used
//             with old code!!!!
//
// 10.2022:    Migrated to faststring3.h (64bit unsigned), CFile2_3.h,

// 12.2024:    Removed all stats. Sequence stats code has been refactored and
//             and cleaned. Stats are only needed in very few applications
//             and should only be computed in these applications.
//
// 03.2025:    remove_gap_ambig_only_sites     


// TODO: Data type handling is/was not optimal.
//       In function set_taxon_and_sequence the auto_detect_type always superset the
//       value the user specified when creating the CSequence_Mol object.
//       This caused unwanted changes of the data type.


// CFile fail states:
// The CFile fail-states are used to hand back reading failure to users of the CSequence_Mol class.
// fail(): General failure to read a sequence.
// fail_reason1(): Not a valid DNA symbol found while reading a file.
// fail_reason2(): Not a valid Protein/AA symbol found while reading a file.

// Functions provided to read sequences:
//  void readSeqName_ignore_Sequence_data(CFile& infile)
//  void readRawFastaSequence(CFile& infile)                         // No character checks
//  void readRawFastaSequence_toupper(CFile& infile)                 // No character checks
//  void readFastaSequence_toupper_ambig2N_removegaps(CFile& infile) // 
//  void readFastaSequence_toupper_ambig2N_gaps2ambig(CFile& infile) //  
//  void readFastaSequence_toupper_removegaps(CFile& infile)         // In Arbeit
//  void readFastaSequence_toupper_ambig2N_gaps2N(CFile& infile)     // deprecated, backward compatible function name.
//  void readFastaSequence_generic(CFile& infile, processing_flag pflag)  // Should be used for most cases.

#include <vector>
#include <iostream>
#include "faststring3.h"
//#include "CFile/CFile3b.h"
#include "CFile/CFile2_3.h"
#include <cstdio>
#include "basic-DNA-RNA-AA-routines.h"
#include <cstring>
#include <cstdio>

class CSequence_Mol;

bool less_than_full_seqname_using_pointer(CSequence_Mol *, CSequence_Mol *);
bool less_than_full_seqname_caseinsensitive_using_pointer(CSequence_Mol *, CSequence_Mol *);
bool less_than_full_seqname_first_unsigned_using_pointer(CSequence_Mol *, CSequence_Mol *);

enum SequenceDataTypesEnum {
   SeqType_dna,      /* DNA sequences (states A, C, G, T) */
   SeqType_rna,      /* RNA sequences (states A, C, G, U) */
   SeqType_protein,      /* amino acid sequences         */
   SeqType_molecular,          /* AA, dna, or rna              */
   SeqType_unknown,
   SeqType_auto_detect_type,
   SeqType_mixed
};


struct basic_sequence_stats
{
   SequenceDataTypesEnum datatype;

   size_t  number_of_AAs;           // Does not include Xs and ?
   size_t  number_of_DNARNA_bases;  // Includes Us but not ambigs
   size_t  number_of_DNARNA_amgibs; // Includes Ns and '?'
   size_t  number_of_gaps;
   size_t  number_of_stars;
   size_t  number_of_Us;

   size_t  number_of_Xs;
   size_t  number_of_Ns;
   size_t  number_raw_original_length; //no spaces, but gaps and stars

   void reset()
   {
     number_of_AAs = 0;
     number_of_DNARNA_bases = 0;
     number_of_DNARNA_amgibs = 0;
     number_of_gaps = 0;
     number_of_stars = 0;
     number_of_Us = 0;
     number_of_Xs = 0;
     number_of_Ns = 0;
     number_raw_original_length = 0;
   }
};



class CSequence_Mol
{
public:
  enum processing_flag {
    convert_toupper                  = 1,
    convert_ambig2N                  = 2,
    autodetecttype                   = 4,
    removegaps                       = 8,
    removestars                      = 16,
    compute_numbers_of_residues_flag = 32,
    gaps_to_ambig                    = 64
  };

private:
  faststring      full_name;   // with description after first space
  faststring      short_name;  // without description after first space
  faststring      description; // description only

  faststring              data;
  SequenceDataTypesEnum   type;
  char                    general_ambig;

  // General rule: If we convert symbols, we count the new symbols we have in memory,
  // not the symbols in the file. This is more consistent.
  //
  //  basic_sequence_stats basic_seq_stats;

  bool bitcode_recoded;

private:
  void find_shortname_and_description()
  {
    full_name.removeSpacesBack();
    full_name.removeSpacesFront();

    faststring::size_type pos_description = full_name.find(' ');
    faststring::size_type full_name_len   = full_name.size();

    if (pos_description == faststring::npos) // no space -> no sequence description
    {
      pos_description = full_name.size();
      short_name.assign(full_name.begin(), full_name.begin() + pos_description);
    }
    else
    {
      short_name.assign(full_name.begin(), full_name.begin() + pos_description);
      // Skip all spaces
      while (pos_description < full_name_len &&
             full_name.get_unckecked(pos_description) == ' ')
      {
        //	std::cerr << "-------------DEBUG: ++pos_description" << std::endl;
        ++pos_description;
      }

      //      std::cerr << "-------------DEBUG: " << full_name.end() - (full_name.begin() + pos_description) << std::endl;
      //      if (full_name.end() - (full_name.begin() + pos_description) > 0) // Should not occur any more!!!!
      //      {
      description.assign(full_name.begin() + pos_description, full_name.end() );
      //      }
    }
  }

  void toupper_this_sequence()
  {
    data.toupper();
  }

  

  // TODO: Keep case of ambig.
  void change_ambigs()
  {
    if (type == SeqType_dna || type == SeqType_rna)
    {
      char *it     = data.begin();
      char *it_end = data.end();

      while (it != it_end)
      {
        char c       = *it;
        if (c == 'n' || c == '?')
        {
          // Convert all ambiguity codes to N's
          *it = 'N';
        }
        ++it;
      }
    }
    else if (type == SeqType_protein)
    {
      char *it     = data.begin();
      char *it_end = data.end();

      while (it != it_end)
      {
        char c       = *it;

        if (c == 'x' || c == '?' )
        {
          // Convert all ambiguity codes to X's
          *it = 'X';
        }
        ++it;
      }
    }

    // TODO: Adjust the numbers
  }

  void remove_gaps_stars_this_sequence(bool rg, bool rs)
  {
    //    char     *it     = data.begin();
    //    char     *it_end = data.end();
    
    
    // TODO

    if (rg)
    {}

    if (rs)
    {}
  }

  void gaps_to_ambig_this_sequence()
  {
    char     *it     = data.begin();
    char     *it_end = data.end();

     if (type == SeqType_dna || type == SeqType_rna)
     {
        while (it != it_end)
        {
           if (*it =='-')
           {
              *it = 'N';
           }
           ++it;
        }
     }
     else if (type == SeqType_protein)
     {
        while (it != it_end)
        {
           if (*it =='-')
           {
              *it = 'X';
           }
           ++it;
        }
     }

  }

  // This function was introduced since it was not clear whether
  // other code should also be added. Since this never happend, the function
  // could be removed and the line could be pasted to where needed?
  void addCharToData(char c)
  {
    data.push_back(c);
  }


  // Reads line of data and removes only white spaces
  // TODO: Make this more efficient by readling a whole line at once.
  void getRawLine(CFile& infile)
  {
    char            c;

    c = infile.getchar();
    while (!infile.fail() && c != '\n')
    {
      if ( isspace(c) )   // We remove white spaces.
      {
        c = infile.getchar();
        continue;
      }
      addCharToData(c);
      c = infile.getchar();
    }
  }

  // Reads line of data and removes only white spaces
  // Reports back the number of non space symbols read.
  // TODO: Make this more efficient by readling a whole line at once.
  void getRawLine(CFile& infile, faststring::size_type &num)
  {
    char               c;
    faststring::size_type local_num = 0;
    
    c = infile.getchar();
    while (!infile.fail() && c != '\n')
    {
      if ( isspace(c) )   // We remove white spaces.
      {
        c = infile.getchar();
        continue;
      }
      addCharToData(c);
      c = infile.getchar();
      ++local_num;
    }
    num = local_num;
  }

  /*
   // Reads line of data and removes only white spaces
   void getRawLine_toupper(CFile& infile)
   {
   char            c;

   c = infile.getchar();
   while (!infile.fail() && c != '\n')
   {
   if ( isspace(c) )   // We remove white spaces.
   {
   c = infile.getchar();
   continue;
   }
   //      addCharToData(toupper(c));
   c = infile.getchar();
   }
   toupper_this_sequence();
   }
   */

  // Called e.g. indirectly from Phobos
  void getLineOfDNA_toupper_ambig2N_removegaps (CFile& infile, faststring::size_type &pos_in_seq)
  {
    char                  c;
    faststring::size_type pos = pos_in_seq; // Lokal variable used to speed up the code.
    char                  local_general_ambig = toupper_char(general_ambig);

    //     ++pos;
    c = infile.getchar();
    while (!infile.fail() && c != '\n')
    {
      c = toupper_char(c);

      if (is_DNA_base(c) || c == 'U')
      {
        addCharToData(c);
      }
      else if ( is_DNA_or_DNA_ambig(c) )
      {
        // Convert all ambiguity codes to N's
        addCharToData(local_general_ambig);
      }
      else if ( isspace(c) || c=='-')
      {
         // Skip spaces, gaps
      }
      else // Not a valid symbol
      {
        // we do not count this -just to be consistent before we bail out
        infile.clear(CFile::__fail_reason1 | CFile::__fail_flag | infile.rdstate()); // Set both fail flags. || 20.7.2021: add rdstate, so that eof flags is not overwritten
        break;
      }
      c = infile.getchar();
      ++pos;
    }
    pos_in_seq = pos;
  }


  void getLineOfDNA_generic (CFile& infile, faststring::size_type &pos_in_seq,
                             processing_flag pflag)
  {
    char                c;
    char                c_upper;
    char                local_general_ambig = general_ambig;
    faststring::size_type  pos = pos_in_seq;
    // signed char     delta;
    bool                convert_toupper_bool  = (pflag & convert_toupper);   // handled
    bool                convert_ambig2N_bool  = (pflag & convert_ambig2N);   // handled
    bool                removegaps_bool       = (pflag & removegaps);        // handled
//    bool                removestars_bool      = (pflag & removestars);       // not relevant for DNA. Stars in the sequence are an error.
    bool                gaps_to_ambig_bool    = (pflag & gaps_to_ambig);     // handled

    // In the routine calling this function it should be checked, that the general_ambig
    // character is upper or lower case, as specified.

    if (convert_toupper_bool)
      local_general_ambig = toupper(general_ambig);

    //     ++pos;
    c = infile.getchar();
    c_upper = toupper_char(c);
    while (!infile.fail() && c != '\n')
    {
      if (convert_toupper_bool)
        c = c_upper;

      if ( is_DNA_base(c) || c_upper == 'U')
      {
        addCharToData(c);
      }
      else if (is_DNA_or_DNA_ambig(c))
      {
        if (c == '?')
        {
           addCharToData('N');
        }
        else if (convert_ambig2N_bool)
        {
          // Convert all ambiguity codes to N's, but we preserve the case. 'w' -> 'n' and 'W' -> 'N' etc.
          if (c <='Z') // i.e. upper case.
            addCharToData('N');
          else
            addCharToData('n');
        }
        else
        {
           addCharToData(c);
        }
      }
      else if ( isspace(c) )
      {
      }                           // Spaces are removed
      else if ( c=='-')
      {
        if ( gaps_to_ambig_bool)
        {
          addCharToData(local_general_ambig);
        }
        else if ( removegaps_bool )
        {
        }
        else
        {
          addCharToData('-');
        }
      }
      else
      {
        // we do not count this -just to be consistent before we bail out
        infile.clear(CFile::__fail_reason1 | CFile::__fail_flag | infile.rdstate()); // Set both fail flags. || 20.7.2021: add rdstate, so that eof flags is not overwritten
        break;
      }
      c = infile.getchar();
      c_upper = toupper_char(c);
      ++pos;
    }
    pos_in_seq = pos;
  }


  void getLineOfProtein_generic (CFile& infile, faststring::size_type &pos_in_seq,
                                 processing_flag pflag)
  {
    char                c;
    char                c_upper;
    char                local_general_ambig = general_ambig;

    faststring::size_type  pos = pos_in_seq;
    bool    convert_toupper_bool  = (pflag & convert_toupper);   // handled
                                                                 //    bool            convert_ambig2N_bool  = (pflag & convert_ambig2N);   // handled
    bool    removegaps_bool       = (pflag & removegaps);        // handled
                                                                 //    bool            removestars_bool      = (pflag & removestars);       // handled
    bool    gaps_to_ambig_bool    = (pflag & gaps_to_ambig);// handled

    // In the routine calling this function it should be checked, that the general_ambig
    // character is upper or lower case, as specified.

    if (convert_toupper_bool)
      local_general_ambig = toupper_char(general_ambig);

    //     ++pos;
    c = infile.getchar();
    c_upper = toupper_char(c);
    while (!infile.fail() && c != '\n')
    {
      if (convert_toupper_bool)
        c = c_upper;

      if (c == '?' || c_upper == 'X' )
      {
        // c is '?' or 'X' or 'x' and general_ambig is 'X' or 'x'
        {
          if (c == '?')
          {
            addCharToData(local_general_ambig);
          }
          else // c is 'X' 'x', general_ambig is 'X' 'x'
            addCharToData(c); // Append X or x and preserve case !!! to be consistent with the flags.
                              // Note: In case of convert to upper, we have already converted c to upper.
        }
      }
      else // The following code assumes we do not have X or ?. Otherwise Xs are duplicated in the data.
        if ( is_aa_or_aa_ambig_extended(c) ) // includes other ambigs (!'X', !'x' ,!'?')
        {
          addCharToData(c);
        }
        else if ( isspace(c) )
        {
        }                              // Spaces are removed
        else if ( c=='-' || c=='*')    // Lets handle gaps and *s
        {

          if (c=='-')
          {
            if ( gaps_to_ambig_bool)
            {
              addCharToData(local_general_ambig);
            }
            else if ( removegaps_bool )
            {
            }
            else
            {
              addCharToData('-');
            }
          }

          /* 	else if (removestars_bool) */
          /* 	{ */
          /* 	  ++basic_seq_stats.number_of_stars; */
          /* 	} */
          //	c = infile.getchar();
          //	++pos;
        }
        else // Unknown symbol
        {
          // we do not count this -just to be consistent before we bail out
          infile.clear(CFile::__fail_reason2 | CFile::__fail_flag  | infile.rdstate()); // Set both fail flags. || 20.7.2021: add rdstate, so that eof flags is not overwritten
          break;
        }
      c = infile.getchar();
      c_upper = toupper_char(c);
      ++pos;
    }
    pos_in_seq = pos;
  }


  void getLineOfDNA_toupper_removegaps (CFile& infile, faststring::size_type &pos_in_seq)
  {
    char                c;
    faststring::size_type  pos = pos_in_seq;

    //     ++pos;
    c = infile.getchar();
    while (!infile.fail() && c != '\n')
    {
      c = toupper_char(c);

      if ( is_DNA_base(c) || c == 'U')
      {
        addCharToData(c);
      }
      else if (c == '?')
      {
        addCharToData('N');
      }
      else if (is_DNA_or_DNA_ambig(c))
      {
         addCharToData(c);
      }
      else if ( isspace(c) )          // Spaces are removed
      {
      }
      else if ( c=='-' || c=='*')     // What do we do with gaps and *s - Here they are removed, but their positions are stored in order to reconstruct the original coordinates.
      {                               // I don't really know why we tolerate '*' ?
        if (c=='-')
        {
        }
        else
        {
        }
        //	c = infile.getchar();
        //	++pos;
      }
      else
      {
        // Something went wrong -> the last reading event from the file failed
        // -> we set the fail flag, this informs the user that something went wrong
        // We do not count this - just to be consistent before we bail out
        infile.clear(CFile::__fail_reason1 | CFile::__fail_flag  | infile.rdstate()); // Set both fail flags. || 20.7.2021: add rdstate, so that eof flags is not overwritten
        break;
      }
      c = infile.getchar();
      ++pos;
    }
    pos_in_seq = pos;
  }


  void getLineOfProtein_toupper_removegaps (CFile& infile, faststring::size_type &pos_in_seq)
  {
    char                     c;
    faststring::size_type       pos = pos_in_seq;

    //     ++pos;
    c = infile.getchar();
    char local_general_ambig = toupper_char(general_ambig);

    while (!infile.fail() && c != '\n')
    {
       c = toupper_char(c);

      if (c == 'X' || c=='?')
      {
        addCharToData(local_general_ambig);
      }
      else if (is_aa_or_aa_ambig_extended(c) )
      {
        addCharToData(c);
      }
      else if ( isspace(c) )// Spaces are removed
      {
      }
      else if ( c=='-' || c=='*')     // What do we do with gaps and *s - Here they are removed, but their positions are stored in order to reconstruct the original coordinates.
      {
        if (c=='-')
        {
        }
        else
        {
        }
      }
      else
      {
        infile.clear(CFile::__fail_reason2 | CFile::__fail_flag  | infile.rdstate()); // Set both fail flags. || 20.7.2021: add rdstate, so that eof flags is not overwritten
        break;
      }
      c = infile.getchar();
      ++pos;
    }
    pos_in_seq = pos;
  }

  void getLineOfDNA_toupper_ambig2N_gaps2ambig (CFile& infile)
  {
    char    c;

    c = infile.getchar();
    while (!infile.fail() && c != '\n')
    {
       c = toupper_char(c);

      if ( is_DNA_base(c))
      {
        addCharToData(c);
      }
      else if ( c == 'N' || c == '?' || is_DNA_or_DNA_ambig(c) || c=='-')
      {
        addCharToData(general_ambig);
      }
      else if ( c == 'U')
      {
        addCharToData('U');
      }
      else if ( isspace(c) )   // We remove white spaces.
      {                        // I don't really know why we tolerate '*' ?
      }
      else
      {
        infile.clear(CFile::__fail_reason1 | CFile::__fail_flag  | infile.rdstate()); // Set both fail flags. || 20.7.2021: add rdstate, so that eof flags is not overwritten
        return;
      }
      c = infile.getchar();
    }
  }


public:
  CSequence_Mol(SequenceDataTypesEnum dt, char p_general_ambig = '\0'):
  type(dt), general_ambig(p_general_ambig),
  bitcode_recoded(false)
  {
    data.reserve(160);
    if (p_general_ambig == '\0') // Use the typical ambig character
    {
      if      (dt == SeqType_dna)     general_ambig = 'N';
      else if (dt == SeqType_protein) general_ambig = 'X';
    }
  }


  CSequence_Mol(SequenceDataTypesEnum dt, faststring name, faststring::size_type len, char p_general_ambig = '\0'):
  full_name(name),type(dt), general_ambig(p_general_ambig),
  bitcode_recoded(false)
  {
    if (len < faststring::npos)
      ++len;
    data.reserve(len);

    find_shortname_and_description();

    if (p_general_ambig == '\0') // Use the typical ambig character
    {
      if      (dt == SeqType_dna)     general_ambig = 'N';
      else if (dt == SeqType_protein) general_ambig = 'X';
    }
  }


  //**************************************************************************
  // pos1 must be the first column starting with 0.
  // pos2 must be the index after the last column.
  // pos2-pos1 must be the number of bases that are copied to this sequence.
  //**************************************************************************
  CSequence_Mol(const CSequence_Mol& s, faststring::size_type pos1=0, faststring::size_type pos2 = faststring::npos):
  full_name(s.full_name),
  short_name(s.short_name),
  description(s.description),
  type(s.type),
  general_ambig(s.general_ambig),
  bitcode_recoded(s.bitcode_recoded)
  {
    if (pos1 == 0 && pos2 >= s.length() )
      data = s.data;
    else
    {
      if (pos2 > s.length() )
        pos2 = s.length();
      if (pos1 >= pos2)
      {
        return;
      }
      data.assign(s.data, pos1, pos2-pos1);
    }
  }
  
  size_t memory_usage()
  {
    //    std::cerr << "memory this:        " << sizeof(*this) << std::endl;
    //    std::cerr << "full name capacity: " << full_name.capacity() << std::endl;
    //    std::cerr << " data capacity:     " << data.capacity() << std::endl;

    return (unsigned long) (sizeof(*this) + full_name.capacity() + data.capacity());
  }

  size_t memory_usage_detailed(unsigned long &mem_this, unsigned long &mem_name, unsigned long &mem_data)
  {
    mem_this = sizeof(*this);
    mem_name = full_name.capacity();
    mem_data = data.capacity();
    return mem_this + mem_name + mem_data;
  }

  // Should only be applied to dna sequences
  void convert_ambig2N_this_sequence()
  {
    if ( !(type == SeqType_dna || type == SeqType_rna) )
      return;

    // This is a DNA sequence:

    char * it     = data.begin();
    char * it_end = data.end();

    while (it != it_end)
    {
      char c       = *it;

      if (is_DNA_ambig(c))
      {
        // Convert all ambiguity codes to N's
        *it = general_ambig;
      }
      ++it;
    }
  }

  bool has_only_gaps()
  {
    char *it     = data.begin();
    char *it_end = data.end();

    while (it != it_end)
    {
      if (*it != '-')
        return false;
      ++it;
    }
    return true;
  }

  bool has_only_ambigs()
  {
    char *it     = data.begin();
    char *it_end = data.end();

    if (type == SeqType_dna || type == SeqType_rna)
    {
      while (it != it_end)
      {
        if (!is_DNA_ambig(*it))
          return false;
        ++it;
      }
      return true;
    }
    else if (type == SeqType_protein)
    {
      while (it != it_end)
      {
        if (!is_aa_ambig(*it))
          return false;
        ++it;
      }
      return true;
    }
    else // This result is undefined, since the data type is not specified.
    {
      return false;
    }
  }

  bool has_only_ambigs_or_gap()
  {
    //    std::cout << "Call to: has_only_ambigs_or_gap" << std::endl;

    char *it     = data.begin();
    char *it_end = data.end();

    if (type == SeqType_dna || type == SeqType_rna)
    {
      //      std::cout << "Datatype is dna or rna" << std::endl;

      while (it != it_end)
      {
        //		std::cout << "Checking " << *it << std::endl;
        if ( !(is_DNA_ambig(*it) || *it == '-') )
        {
          //	  	  std::cout << "Returning false." << std::endl;
          return false;
        }
        ++it;
      }
      //            std::cout << "Returning true." << std::endl;
      return true;
    }
    else if (type == SeqType_protein)
    {
      while (it != it_end)
      {
        if ( !(is_aa_ambig(*it) || *it == '-') )
          return false;
        ++it;
      }
      return true;
    }
    else // This result is undefined, since the data type is not specified.
    {
      std::cout << "Returning default false for unspecified type." << std::endl;
      return false;
    }
  }


  faststring::size_type length() const
  {
    return data.size();
  }

  SequenceDataTypesEnum get_datatype()
  {
    return type;
  }

  faststring type_as_string()
  {
    if (type == SeqType_dna)
      return faststring("DNA");
    else if (type == SeqType_rna)
      return faststring("RNA");
    else if (type == SeqType_protein)
      return faststring("Protein");
    return "Unknown";
  }

  void set_taxon_and_sequence(const faststring& taxonName,
                              const faststring& seqData, SequenceDataTypesEnum datatype=SeqType_auto_detect_type)
  {
    full_name = taxonName;
    
    find_shortname_and_description();
    /*     if ( (pos = taxonName.find(' '))!=faststring::npos) */
    /*     { */
    /*       short_name  = taxonName.substr(0,pos+1); */
    /*       description = taxonName.substr(pos+1); */
    /*     } */
    /*     else */
    /*     { */
    /*       short_name  = taxonName; */
    /*       description.clear(); */
    /*     } */
    data = seqData;
    
    //
    //    toupper_this_sequence();
    
    // Should not be necessary:
    //    compute_numbers_of_residues();
    if (datatype == SeqType_auto_detect_type || datatype == SeqType_unknown)
    {
//      auto_detect_datatype();
      //    std::cerr << "Data type: " << type_as_string()  << std::endl;
      //    std::cerr << "General ambig: " << general_ambig  << std::endl;
//      if (auto_detect_ambig_detect_conflict() )
//      {
//        std::cerr << "Conflict in ambiguity code used in imported sequence." << std::endl;
//      }
    }
    //    std::cerr << "General ambig: " << general_ambig  << std::endl;
    
    change_ambigs();
    //    remove_gaps_stars_this_sequence(true, true);
    
    //    general_ambig = ambig;
    //    type = t;
  }


  void set_taxon_and_random_DNA_sequence(const faststring& taxonName,
                                         faststring::size_type len,
                                         double piA, double piC, double piG)
  {
    size_t pos;

    type = SeqType_dna;
    general_ambig = 'N';

    data.clear();
    data.reserve(len+1);

    full_name = taxonName;
    find_shortname_and_description();
    if ( (pos = taxonName.find(' '))!=faststring::npos)
    {
      short_name  = taxonName.substr(0,pos+1);
      description = taxonName.substr(pos+1);
    }
    else
    {
      short_name  = taxonName;
      description.clear();
    }

    // Compute random sequence:
    faststring::size_type i;
    double             z;

    double bA = piA;
    double bC = bA+piC;
    double bG = bC+piG;

    for (i=0; i<len; ++i)
    {
      z = rand()/(RAND_MAX+1.0);
      if ( z < bA )
        data.push_back('A');
      else if ( z < bC )
        data.push_back('C');
      else if ( z < bG )
        data.push_back('G');
      else
        data.push_back('T');
    }

    //    std::cerr << "Data type: " << type_as_string()  << std::endl;
    //    std::cerr << "General ambig: " << general_ambig  << std::endl;
    //    std::cerr << "General ambig: " << general_ambig  << std::endl;
  }


  char& operator[](faststring::size_type pos)
  {
    return data[pos];
  }

  char get_pos(faststring::size_type pos)
  {
    return data[pos];
  }

  void trim_seq_name(faststring::size_type max_len)
  {
    full_name.shorten(max_len);
    find_shortname_and_description();
  }

  void trim_seq_name(const char *trim_at_these_symbols)
  {
    full_name.shorten_to_first_occurrence_of(trim_at_these_symbols);
    find_shortname_and_description();
  }
  
  /*   void set_pos(faststring::size_type pos, char c) */
  /*   { */
  /*     data[pos] = c; */
  /*   } */

  bool is_bitcoderecoded()
  {
    return bitcode_recoded;
  }

/*
  faststring::size_type full_original_length() const
  {
    // Remember that some numbers include others:
    // basic_seq_stats.number_of_DNARNA_amgibs includes 'Ns' and '?', so they are not added again:
    // basic_seq_stats.number_of_DNARNA_bases  includes 'Us'
    if (type == SeqType_dna || type == SeqType_rna )
    {
      return basic_seq_stats.number_of_DNARNA_bases + basic_seq_stats.number_of_DNARNA_amgibs + basic_seq_stats.number_of_gaps + basic_seq_stats.number_of_stars;
    }
    else if (type == SeqType_protein)
    {
      return basic_seq_stats.number_of_AAs + basic_seq_stats.number_of_gaps + basic_seq_stats.number_of_stars+ basic_seq_stats.number_of_Xs;
    }
    else
    {
      return basic_seq_stats.number_raw_original_length;
    }
  }
*/

   /*
  faststring::size_type length_DNARNA_bases() const
  {
    return basic_seq_stats.number_of_DNARNA_bases;
  }
*/

  const char* getName() // const
  {
    return short_name.c_str();
  }

  const faststring& getName_faststring() const
  {
    return short_name;
  }

  const faststring& getSeq_faststring() const
  {
    return data;
  }

  faststring getPhylipName(int mlen)
  {
    faststring str(short_name);

    str.resize(mlen, ' ');
    return str;
  }

  const char* getFullName() // const
  {
    return full_name.c_str();
  }

  const faststring getFullName_faststring() const
  {
    return full_name;
  }

  const char* getDescription() // const
  {
    return description.c_str();
  }

  const char* getSeqStr()
  {
    return data.c_str();
  }

  faststring getPartialSeq(faststring::size_type pos, faststring::size_type n)
  {
    return data.substr(pos, n);
  }

  const char* getSeqBegin()
  {
    return data.begin();
  }

  const char* getSeqEnd()
  {
    return data.end();
  }

  const char* getSeqRBegin()
  {
    return data.rbegin();
  }

  const char*getSeqREnd()
  {
    return data.rend();
  }

/*
  void compute_numbers_of_residues()
  {
    const char *it=data.begin(), *it_end=data.end();

    if (numbers_of_residues_computed)
    {
      std::cerr << "Request for computing number of residues - even though they have already been determined" << std::endl;
      return;
    }

    reset_numbers_of_residues();

    for (;it<it_end; ++it)
    {
      char c = *it;

      if (is_DNA_or_RNA_base(c))
      {
         ++basic_seq_stats.number_of_DNARNA_bases;
      }

      if (is_DNA_ambig(c))
        ++basic_seq_stats.number_of_DNARNA_amgibs;

      if (is_aa_or_aa_ambig_extended(c) && c!= 'X' && c != '?')
        ++basic_seq_stats.number_of_AAs;

      if (c == 'N')
      {
        ++basic_seq_stats.number_of_Ns;
      }

      if (c == '-')
      {
        ++basic_seq_stats.number_of_gaps;
      }

      if (c == '*')
      {
        ++basic_seq_stats.number_of_stars;
      }

      if (c == 'X')
      {
        ++basic_seq_stats.number_of_Xs;
      }

      numbers_of_residues_computed = true;
    }
  }

  void auto_detect_datatype()
  {
    if (!numbers_of_residues_computed)
      compute_numbers_of_residues();

    faststring::size_type aa_stuff     = basic_seq_stats.number_of_AAs + basic_seq_stats.number_of_Xs + basic_seq_stats.number_of_stars;
    faststring::size_type DNARNA_stuff = basic_seq_stats.number_of_DNARNA_bases + basic_seq_stats.number_of_DNARNA_amgibs;

    // reine DNA: 100 > 0.25 * 100 - false - DNA

    if (aa_stuff < 0.1*data.length() && DNARNA_stuff < 0.1*data.length() )
      type = SeqType_unknown;
    else if (aa_stuff > 0.25*DNARNA_stuff)
    {
//      if (number_of_Us > 0)
//        type = rna;
//      else
        type = SeqType_dna;
    }
    else
    {
      type = SeqType_protein;
    }
  }

  bool auto_detect_ambig_detect_conflict()  // Return true in case of ambig conflict, else false
  {
    if (!numbers_of_residues_computed)
      compute_numbers_of_residues();

    if (type == SeqType_protein)
      general_ambig = 'X';
    else general_ambig = 'N';

    // Determine conlict:
    return ((basic_seq_stats.number_of_Xs > 0 || ((type==SeqType_dna || type==SeqType_rna) && basic_seq_stats.number_of_Ns > 0)));
  }
*/

  void assignSeq(faststring &newseq, processing_flag flag=(processing_flag)0)
  {
    data.assign(newseq);
  }

  void maskSeq(const char *start, const char *end, char c)
  {
    char *pos = const_cast<char*>(start);
    while (pos != end)
    {
      if (*pos != '-')
        *pos = c;
      ++pos;
    }
  }

  void maskSeq_also_mask_gaps(const char *start, const char *end, char c)
  {
    char *pos = const_cast<char*>(start);
    while (pos != end)
    {
      *pos = c;
      ++pos;
    }
  }

  void maskSeq_tolower(const char *start, const char *end)
  {
    char *pos = const_cast<char*>(start);
    while (pos != end)
    {
      *pos = tolower(*pos);
      ++pos;
    }
  }

  faststring::size_type find_pos_first_symbol_not_in_alphabet(const char * alphabet=NULL)
  {
    if (alphabet == NULL)
    {
      if (type == SeqType_dna || SeqType_rna)
      {
        alphabet = DNARNA_symbols;
      }
      else if (type == SeqType_protein)
      {
        alphabet = aa_symbols;
      }
      else
      {
        return faststring::npos;
      }
    }

    return data.find_first_not_of(alphabet);
  }

  
  void determine_bounds_with_respect_to_gaps(faststring::size_type &b, faststring::size_type &e)
  {
    //    std::cerr << "This name :" << getName() << std::endl;
    //    std::cerr << data << std::endl;
    b = data.find_first_not_of('-'); // 0 based
    e = data.find_last_not_of('-');  // 0 based
  }

  // Fills str with the sequence data, that contains gaps and stars, that had been removed
  // temporarily.

  

  // Replace part of the sequence. pos1 and pos2 are in the destination sequence.
  // - erases region from pos1 to pos2
  //   pos1 is the first column pos2-pos1 is the number of bases that are removed,
  //   i.e. pos2 is the column after the last column that is to be removed.
  // - insert s with full length.
  void replace_part_of_sequence(CSequence_Mol &s, faststring::size_type pos1, faststring::size_type pos2)
  {
    if (type != s.type)
    {
      std::cerr << "Critical warning: Replacing part of sequence with different data type." << std::endl;
    }
    data.replace(pos1, pos2-pos1, s.data);
  }

  // Replace first with second, third with fourth char, etc replace string.
  void recode_with_replacement_string(faststring replace)
  {
    data.replace_chars(replace);
  }

  void degap()
  {
    data.remove_symbol('-');
  }

  void remove_with_vector(const std::vector<bool>& mask)
  {
    data.remove_with_vector(mask);
  }

  void append_sequence(CSequence_Mol &s)
  {
    if (type != s.type)
    {
      std::cerr << "Critical warning: Appending sequence with different data type." << std::endl;
    }
    data.append(s.data);
  }

  void append_residue_unchecked(char c)
  {
    data.push_back(c);
  }

/*
  void writeSequence_fill_in_gaps_and_stars(FILE *of, unsigned char_per_line=50)
  {
    if (bitcode_recoded)
      backrecode_bitcode_sequence();

    faststring::size_type upos = 0;
    const char *pos = data.begin();
    //    const char *end = data.end();

    faststring::size_type left_this_line = char_per_line;

    faststring::size_type N_all_char = data.size()+N;

    fprintf(of, ">%s\n", full_name.c_str());

    while (upos < N_all_char)
    {
      if (pos_in_gap_and_star_koords == N)
        break;

      {
        putc(*pos, of);
        ++pos;
      }
      ++upos;
      --left_this_line;

      if (left_this_line == 0)
      {
        putc('\n', of);

        left_this_line = char_per_line;
      }
    }

    while (upos < N_all_char)
    {
      putc(*pos, of);
      ++pos;
      ++upos;
      --left_this_line;

      if (left_this_line == 0)
      {
        putc('\n', of);

        left_this_line = char_per_line;
      }
    }
    if (left_this_line != char_per_line)    // We have already written a partial line
      putc('\n', of);
  }
*/

  void writeSequence(FILE *of, unsigned char_per_line=50) //const
  {
    writeSequence_fasta(of, 0, faststring::npos, full_name.c_str(), char_per_line);
  }

  void writeSequence_fasta(FILE *of, unsigned char_per_line=50) //const
  {
    writeSequence_fasta(of, 0, faststring::npos, full_name.c_str(), char_per_line);
  }

  void writeSequence_fasta_revComp(FILE *of, unsigned char_per_line=50) //const
  {
    writeSequence_fasta_revComp(of, 0, faststring::npos, full_name.c_str(), char_per_line);
  }

  void writeSequence(FILE *of,
                     faststring::size_type  pos_beg,
                     faststring::size_type  pos_end,
                     const char  *s_name,
                     unsigned    char_per_line=50) const
  {
    writeSequence_fasta(of, pos_beg, pos_end, s_name, char_per_line);
  }


  void writeSequence_fasta(FILE *of,
                           faststring::size_type    pos_beg, // 0 based
                           faststring::size_type    pos_end, // 0 based, stop index
                           const char               *s_name,
                           faststring::size_type    char_per_line=50) const
  {
    const char *beg = data.begin() + pos_beg;
    const char *end = data.begin() + pos_end;
    const char *tmp_end;

    /*     if (bitcode_recoded) */
    /*       backrecode_bitcode_sequence(); */

    if (pos_end == faststring::npos || end > data.end() )
      end = data.end();

    if (char_per_line == faststring::npos)
      char_per_line = data.length();

    fprintf(of, ">%s\n", s_name);
    
    tmp_end = beg + char_per_line;
    if (tmp_end > end)
      tmp_end = end;

    while (beg < end)
    {
      fprintf(of, "%.*s\n", int(tmp_end - beg), beg);

      beg = tmp_end;
      tmp_end = beg + char_per_line;
      if (tmp_end > end)
        tmp_end = end;
    }
  }

  // Some programs do not accept too many Ns in a row.
  // This function export sequences so that N-regions are collapsed
  // to a maximum number of successive Ns.
  // USE WITH CAUTION: This function changes your sequence data!!!
  // The output is not identical to your original sequence data.


  void writeSequence_fasta_collaps_NX_regions(FILE *of,
                                              faststring::size_type char_per_line=50,
                                              faststring::size_type maximum_number_of_successive_NXs=100)
  {
    writeSequence_fasta_collaps_NX_regions(of,
                                           0,
                                           faststring::npos,
                                           full_name.c_str(),
                                           char_per_line,
                                           maximum_number_of_successive_NXs);
  }

  void writeSequence_fasta_collaps_NX_regions(FILE *of,
                                              faststring::size_type  pos_beg, // 0 based
                                              faststring::size_type  pos_end, // 0 based , stop index
                                              const char          *s_name,
                                              faststring::size_type  char_per_line=50,
                                              faststring::size_type  maximum_number_of_successive_NXs=100) const
  {
    const char *beg = data.begin() + pos_beg;
    const char *end = data.begin() + pos_end;
    //    const char *tmp_end;

    /*     if (bitcode_recoded) */
    /*       backrecode_bitcode_sequence(); */

    if (pos_end == faststring::npos || end > data.end() )
      end = data.end();

    fprintf(of, ">%s\n", s_name);

    faststring::size_type chars_in_this_line=0;
    faststring::size_type tmp_succNs=0;

    while (beg != end)
    {
      char c = *beg;
      char c_toupper = toupper_char(c);

      if (c_toupper == 'N' || c_toupper == 'X')
      {
        ++tmp_succNs;
        if (tmp_succNs <= maximum_number_of_successive_NXs) // print the N or X
        {
          putc(c, of);
          ++chars_in_this_line;

          if (chars_in_this_line == char_per_line)
          {
            putc('\n', of);
            chars_in_this_line = 0;
          }
        }
        // unwritten else: do not print the N or X
      }
      else
      {
        putc(c, of);
        tmp_succNs = 0;
        ++chars_in_this_line;

        if (chars_in_this_line == char_per_line)
        {
          putc('\n', of);
          chars_in_this_line = 0;
        }
      }
      ++beg;
    }
    if (chars_in_this_line != 0)
      putc('\n', of);
  }

  void writeSequence_partial_fasta(FILE *of,
                                   faststring::size_type pos_beg, // 0 based
                                   faststring::size_type pos_end, // 0 based , stop index
                                   const char         *s_name,
                                   faststring::size_type           char_per_line=50,
                                   bool               degap=true) const
  {
    const char *beg = data.begin() + pos_beg;
    const char *end = data.begin() + pos_end;

    /*     if (bitcode_recoded) */
    /*       backrecode_bitcode_sequence(); */

    if (pos_end == faststring::npos || end > data.end() )
      end = data.end();

    fprintf(of, ">%s\n", s_name);

    faststring::size_type chars_in_this_line=0;
    
    if (degap)
    {
      while (beg != end)
      {
        char c = *beg;

        if (c != '-')
        {
          putc(c, of);
          ++chars_in_this_line;

          if (chars_in_this_line == char_per_line)
          {
            putc('\n', of);
            chars_in_this_line = 0;
          }
        }
        // unwritten else: do not print -
        ++beg;
      }
    }
    else
    {
      while (beg != end)
      {
        putc(*beg, of);
        ++chars_in_this_line;

        if (chars_in_this_line == char_per_line)
        {
          putc('\n', of);
          chars_in_this_line = 0;
        }
        ++beg;
      }
    }
  }

  void writeSequence_fasta_filter_sites(FILE *of,
                                        std::vector<faststring::size_type> &site_filter,
                                        int         offset,   // 0 if site_filter is 0 based, 1 if site_filter is 1 based
                                        const char  *s_name,
                                        faststring::size_type    char_per_line=50) const
  {
    const char *beg = data.begin();
    const char *end = data.end();

    /*     if (bitcode_recoded) */
    /*       backrecode_bitcode_sequence(); */

    fprintf(of, ">%s\n", s_name);

    faststring::size_type chars_in_this_line=0;
    int i=offset;

    {
      while (beg != end)
      {
        if (site_filter[i] != 0)
        {
          putc(*beg, of);
          ++chars_in_this_line;

          if (chars_in_this_line == char_per_line)
          {
            putc('\n', of);
            chars_in_this_line = 0;
          }
        }
        // Move on in sequence and site_filter vector
        ++beg;
        ++i;
      }
    }
  }

  void writeSequence_phylip(FILE *of,
                            faststring::size_type pos_beg, // 0 based
                            faststring::size_type pos_end, // 0 based, stop index.
                            const char  *s_name) const
  {
    /*     if (bitcode_recoded) */
    /*       backrecode_bitcode_sequence(); */

    const char *beg = data.begin() + pos_beg;
    // const char *end = data.begin() + pos_end;

    /*     if (pos_end == faststring::npos || end > data.end() ) */
    /*       end = data.end(); */

    faststring::size_type name_len = strlen(s_name);
    if (name_len < 11)
      fprintf(of, "%-10.10s", s_name);
    /*     else */
    /*     { */
    /*       char str[11]; */
    /*       std::strncpy(str, s_name, 10); */
    /*       str[10] = '\0'; */
    /*       fprintf(of, "%s", str); */
    /*     } */
    fprintf(of, "%.*s\n", (int)(pos_end - pos_beg), beg);
  }

  // Expects 0 bases coordinates. The range starts at pos_beg and pos_end
  // points to the base after the last position in the range.
  // Clearly pos_end > pos_beg. In the empty range, pos_end = pos_beg+1.

  // Deprecated name
  void writeSequence_revComp(FILE *of,
                             faststring::size_type pos_beg,  // 0 based index
                             faststring::size_type pos_end,  // 0 based index after the last pos.
                             const char *s_name,
                             faststring::size_type    char_per_line=50) const
  {
    writeSequence_fasta_revComp(of, pos_beg, pos_end, s_name, char_per_line);
  }

  void writeSequence_fasta_revComp(FILE *of,
                                   faststring::size_type pos_beg,  // 0 based index
                                   faststring::size_type pos_end,  // 0 based index after the last pos.
                                   const char *s_name,
                                   faststring::size_type    char_per_line=50) const
  {
    /*     if (bitcode_recoded) */
    /*       backrecode_bitcode_sequence(); */

    const char *beg = data.begin() + pos_beg;
    const char *end = data.begin() + pos_end;

    faststring::size_type interleaved = char_per_line;

    fprintf(of, ">%s\n", s_name);

    
    // Walking backwards, we shift the range:
    --beg;   // 0 based indices. beg is the rend index. Move beg to pos before "first pos". This is the stop index.
    --end;   // 0 based index. Move end to index at which we start to move.

    // Example: Print first 3 bases.
    // Then: pos_beg=0, pos_end=3, 0 based, pos_end is index behind last pos in range.
    // Thus: beg="-1+seq_beg" and end="2+seq_beg", so print 2, 1, 0, and stop at -1

    for (; beg < end; --end)
    {
      if (interleaved == 0)
      {
        fputc('\n', of);
        interleaved = char_per_line;
      }
      fputc( complementBase(*end), of ); // From CDnaString.h
      --interleaved;
    }
    fputc('\n', of);
  }

  faststring::size_type baseComposition(std::vector<unsigned> &v) const
  {
    if ( v.size() < 256 )
    {
      std::cerr << "Internal error. -24" << std::endl;
      exit(0);
    }

    char *pos = data.begin();
    char *end = data.end();

    faststring::size_type len = end - pos;

    while (pos != end)
    {
      ++v[*pos];
      ++pos;
    }
    return len;
  }

  void sliding_window_CG_AT_content(std::vector<double> &v, unsigned window_size)
  {
    char *end = data.end();

    const char *range_beg = data.begin();
    const char *range_end;

    double   CG_prop;
    unsigned CG_num;
    unsigned AT_num;

    while (1)
    {
      range_end = range_beg + window_size;

      if (range_end > end)
        range_end = end;

      if (range_end == range_beg)
        break;

      CG_AT_content_in_region(range_beg, range_end, AT_num, CG_num);

      //      std::cerr << AT_num << " " << CG_num << std::endl;

      if (CG_num == 0 && AT_num == 0) // e.g. if all chars are ambigs
        CG_prop = -1;
      else
        CG_prop = (double)CG_num/(CG_num+AT_num);
      
      v.push_back(CG_prop);

      range_beg = range_end;
    }
  }


  unsigned long  AsciiComposition_CpG_CpNpG(unsigned long  v[], unsigned long  &CpG, unsigned long  &CpNpG,
                                            unsigned long  &Nislands) const
  {
    char *pos = data.begin();
    char *end = data.end();
    char last = 'N';  // An island is not counted if located at the beginning or end of the sequence.
    
    faststring::size_type len = end - pos;
    
    
    end = end - 2;
    while (pos < end)
    {
      char c1 = *pos;
      
      if ( c1 == 'C' || c1 == 'c' )
      {
        char c2 = *(pos+1);
        char c3 = *(pos+2);
        if ( c2 == 'G' || c2 == 'g') ++CpG;
        if ( c3 == 'G' || c3 == 'g') ++CpNpG;
      }
      if (c1 == 'N' && last != 'N')
        ++Nislands;
      
      ++v[(int)c1];
      ++pos;
      last = c1;
    }
    
    // Now deal with the second last position:
    if (len > 1 )
    {
      char c1 = *pos;
      
      if ( c1 == 'C' || c1 == 'c' )
      {
        char c2 = *(pos+1);
        if ( c2 == 'G' || c2 == 'g' ) ++CpG;
      }
      if (c1 == 'N' && last != 'N')
        ++Nislands;
      ++v[(int)c1];
      ++pos;
      last = c1;
    }
    
    // Now deal with the last position:
    if (len > 0)
    {
      ++v[(int)*pos];
      if (*pos == 'N') // We don't count the island if at the end of the sequence
      {
        if  (last == 'N') // it has already been counted so undo this
        {
          --Nislands;
        }
        // else
        //    we don't count it
      }
    }
    return len;
  }

/*
  void debug_print_gap_coords(std::ostream &os)
  {
    faststring::size_type *git = gap_and_star_koords.begin();

    os << "Gap positions:" << std::endl;
    for (; git < gap_and_star_koords.end(); ++git)
      os << *git << std::endl;
  }
*/

  faststring::size_type getSeqCoord(const char *seq_pos)
  {
    return seq_pos - data.begin() + 1;
  }

  char get_ambiguity_character()
  {
    return general_ambig;
  }

/*
  faststring::size_type getSeqCoord_gap_correted(const char *seq_pos, faststring::size_type &seq_gap_offset)
  {
    faststring::size_type corrected_pos = seq_pos - data.begin()+1+seq_gap_offset;
    faststring::size_type i = seq_gap_offset;
    faststring::size_type N = gap_and_star_koords.size();

    //   std::cerr << "Entering getSeqCoord_gap_correted i: " << i << " N: " << N << " corr coord: " << corrected_pos << std::endl;

    while (i < N && gap_and_star_koords.get_unckecked((unsigned) i) < corrected_pos)
    {
      ++i;
      ++corrected_pos;
    }

    seq_gap_offset = i;

    //    std::cerr << "Exiting getSeqCoord_gap_correted i: " << i << " : " << N << " corr coord: " << corrected_pos << std::endl;

    return corrected_pos;
  }
*/

  void writeSequence(std::ostream &os, faststring::size_type char_per_line)
  {
    char *pos = data.begin();
    char *end = data.end();
    char *tmp_end;

    if (bitcode_recoded)
      backrecode_bitcode_sequence();

    os << ">" << full_name.c_str() << std::endl;
    //    os << ">" << full_name.c_str() << std::endl;

    tmp_end = pos + char_per_line;
    if (tmp_end > end)
      tmp_end = end;

    while (pos != end)
    {
      while (pos != tmp_end)
      {
        os << *pos;
        ++pos;
      }
      os << std::endl;
      tmp_end = pos + char_per_line;
      if (tmp_end > end)
        tmp_end = end;
    }
  }

  /*   void set_to_using_phylip_sequential_sequence() */
  /*   { */
  /*   } */

  /* REVISE THIS FUNCTION. DEFINE WHAT PARAMETER num DOES.
   void readRawFastaSequence(CFile& infile, faststring::size_type &num)
   {
   char     c   = '\0';
   bool     first_line = true;0

   full_name = "";
   data.clear();

   reset_numbers_of_residues();

   c = infile.getchar();
   while (!infile.fail() && isspace(c)) // Ignore spaces
   c = infile.getchar();
   while (!infile.fail() && c == ';')   // Ignore lines starting with ;
   {
   infile.ignore('\n');
   c = infile.getchar();
   while (!infile.fail() && isspace(c) )   // Ignore spaces
   c = infile.getchar();
   }

   if (c != '>')          // Should be new sequence
   {
   infile.clear(CFile::__fail_flag); // Sets the fail bit of the input stream to indicate that reading failed.
   return;
   }
   infile.getline(full_name);
   if (infile.fail()) // getline only fails, if 0 characters could be read.
   {
   return;
   }
   find_shortname_and_description();

   c = infile.peekchar();
   while ( !infile.fail() && c != '>' )
   {
   //      cerr << data.capacity() << endl;
   if (first_line)
   getRawLine(infile, num);
   else
   getRawLine(infile);
   first_line = false;
   c = infile.peekchar();
   }
   if (data.length() > 0)
   infile.clear(infile.rdstate() & ~CFile::__fail_flag); // Unset the fail flag.
   }
   */

  void readRawFastaSequence(CFile& infile)
  {
    char     c   = '\0';

    full_name = "";
    data.clear();

    c = infile.getchar();
    while (!infile.fail() && isspace(c)) // Ignore leading spaces
      c = infile.getchar();
    while (!infile.fail() && c == ';')   // Ignore lines starting with ;
    {
      infile.ignore('\n');
      c = infile.getchar();
      while (!infile.fail() && isspace(c) )   // Ignore spaces
        c = infile.getchar();
    }

    if (c != '>')          // Should be new sequence
    {
      infile.clear(CFile::__fail_flag | infile.rdstate()); // 20.7.2021: add rdstate, so that eof flags is not overwritten
      return;
    }
    infile.getline(full_name);
    if (infile.fail())
    {
      return;
    }
    find_shortname_and_description();

    c = infile.peekchar();
    while ( !infile.fail() && c != '>' )
    {
      //      cerr << data.capacity() << endl;
      getRawLine(infile);
      c = infile.peekchar();
    }
    if (data.length() > 0)
      infile.clear(infile.rdstate() & ~CFile::__fail_flag); // Unset the fail flag.
  }

  void readSeqName_ignore_Sequence_data(CFile& infile)
  {
    char     c   = '\0';

    full_name = "";
    data.clear();

    c = infile.getchar();
    while (!infile.fail() && isspace(c)) // Ignore spaces
      c = infile.getchar();
    while (!infile.fail() && c == ';')   // Ignore lines starting with ;
    {
      infile.ignore('\n');
      c = infile.getchar();
      while (!infile.fail() && isspace(c) )   // Ignore spaces
        c = infile.getchar();
    }

    if (c != '>')          // Should be new sequence
    {
      infile.clear(CFile::__fail_flag | infile.rdstate()); // 20.7.2021: add rdstate, so that eof flags is not overwritten
      return;
    }
    infile.getline(full_name);
    if (infile.fail())
    {
      return; // getline only fails, if 0 characters could be read
    }
    find_shortname_and_description();

    // Ignore the sequence data of this sequence:
    c = infile.peekchar();
    while ( !infile.fail() && c != '>' )
    {
      //      std::cerr << data.capacity() << std::endl;
      infile.ignore('\n');
      c = infile.peekchar();
    }
    // Should not fail, if we get here, since we successfully read a sequence name.
    infile.clear(infile.rdstate() & ~CFile::__fail_flag); // Unset the fail flag.
  }


  void readRawFastaSequence_toupper(CFile& infile)
  {
    char     c   = '\0';

    full_name = "";
    data.clear();

    c = infile.getchar();
    while (!infile.fail() && isspace(c))
      c = infile.getchar();
    while (!infile.fail() && c == ';')
    {
      infile.ignore('\n');
      c = infile.getchar();
      while (!infile.fail() && isspace(c) )
        c = infile.getchar();
    }

    if (c != '>')
    {
      infile.clear(CFile::__fail_flag | infile.rdstate()); // 20.7.2021: add rdstate, so that eof flags is not overwritten
      return;
    }
    infile.getline(full_name);
    if (infile.fail())
    {
      return;
    }
    find_shortname_and_description();

    c = infile.peekchar();
    while ( !infile.fail() && c != '>' )
    {
      //      std::cerr << data.capacity() << std::endl;
      getRawLine(infile);
      c = infile.peekchar();
    }
    toupper_this_sequence();
    if (data.length() > 0)
      infile.clear(infile.rdstate() & ~CFile::__fail_flag); // Unset the fail flag.
  }

  // it1 is the first position to fill. This cooridnate must be 0 based.
  // it2 is the position behind the last position to fill. It must be 0 based.

  void fill_range_with(faststring::size_type pos1, faststring::size_type pos2, char c)
  {
    char     *it1     = data.begin() + pos1;
    char     *it2     = data.begin() + pos2;
    char     *it_end  = data.end();

    if (it2 > it_end)
      it2 = it_end;

    while (it1  < it2)
    {
      *it1 = c;
      ++it1;
    }
  }

  //////*********************** Continue here:

  // Called from Phobos:
  // Stars should only be tolerated in protein sequences.
  void readFastaSequence_toupper_ambig2N_removegaps(CFile& infile)
  {
    char                c   = '\0';
    faststring::size_type  pos = 0;

    full_name = "";
    data.clear();

    c = infile.getchar();
    while (!infile.fail() && isspace(c))
      c = infile.getchar();
    while (!infile.fail() && c == ';')
    {
      infile.ignore('\n');
      c = infile.getchar();
      while (!infile.fail() && isspace(c) )
        c = infile.getchar();
    }

    if (c != '>')
    {
      infile.clear(CFile::__fail_flag | infile.rdstate()); // 20.7.2021: add rdstate, so that eof flags is not overwritten
      return;
    }
    infile.getline(full_name);
    if (infile.fail())
    {
      return;
    }
    find_shortname_and_description();

    c = infile.peekchar();
    while ( !infile.fail() && c != '>' )
    {
      //      std::cerr << data.capacity() << std::endl;
      if (type == SeqType_dna)
        getLineOfDNA_toupper_ambig2N_removegaps(infile, pos);
      else
        std::cout << "Not implemented. Use the readFastaSequence_toupper_removegaps(CFile& infile) function instead."
        << std::endl;
      c = infile.peekchar();
    }
    if (data.length() > 0)
      if ( !(infile.fail_reason1() || infile.fail_reason2() ) )
        infile.clear(infile.rdstate() & ~CFile::__fail_flag); // Unset the fail flag.
  }


  void readFastaSequence_toupper_removegaps(CFile& infile)
  {
    char                c   = '\0';
    faststring::size_type  pos = 0;

    full_name = "";
    data.clear();

    c = infile.getchar();
    while (!infile.fail() && isspace(c))
      c = infile.getchar();
    while (!infile.fail() && c == ';')
    {
      infile.ignore('\n');
      c = infile.getchar();
      while (!infile.fail() && isspace(c) )
        c = infile.getchar();
    }

    if (c != '>')
    {
      infile.clear(CFile::__fail_flag | infile.rdstate()); // 20.7.2021: add rdstate, so that eof flags is not overwritten
      return;
    }
    infile.getline(full_name);
    if (infile.fail())
    {
      return;
    }
    find_shortname_and_description();

    c = infile.peekchar();
    while ( !infile.fail() && c != '>' )
    {
      //      std::cerr << data.capacity() << std::endl;
      if (type == SeqType_dna)
        getLineOfDNA_toupper_removegaps(infile, pos);
      else
        getLineOfProtein_toupper_removegaps(infile, pos);
      c = infile.peekchar();
    }
    if (data.length() > 0)
      if ( !(infile.fail_reason1() || infile.fail_reason2() ) )
        infile.clear(infile.rdstate() & ~CFile::__fail_flag); // Unset the fail flag.

  }

  // Backward compatible function:
  // Called e.g. from Phobos
  void readFastaSequence_toupper_ambig2N_gaps2N(CFile& infile)
  {
    readFastaSequence_toupper_ambig2N_gaps2ambig(infile);
  }


  // Has recently been renamed: gaps2N -> gaps2ambig
  // since we want to prepare for reading proteins.
  void readFastaSequence_toupper_ambig2N_gaps2ambig(CFile& infile)
  {
    char     c   = '\0';

    full_name = "";
    data.clear();

    c = infile.getchar();
    while (!infile.fail() && isspace(c))
      c = infile.getchar();
    while (!infile.fail() && c == '#')
    {
      infile.ignore('\n');
      c = infile.getchar();
      while (!infile.fail() && isspace(c) )
        c = infile.getchar();
    }

    if (c != '>')
    {
      infile.clear(CFile::__fail_flag | infile.rdstate()); // 20.7.2021: add rdstate, so that eof flags is not overwritten
      return;
    }
    infile.getline(full_name);
    if (infile.fail())
    {
      return;
    }
    find_shortname_and_description();

    c = infile.peekchar();
    while ( !infile.fail() && c != '>' )
    {
      //      std::cerr << data.capacity() << std::endl;
      if (type == SeqType_dna)
        getLineOfDNA_toupper_ambig2N_gaps2ambig(infile);
      else
        std::cout << "Not implemented." << std::endl;
      c = infile.peekchar();
    }
    if (data.length() > 0)
      if ( !(infile.fail_reason1() || infile.fail_reason2() ) )
        infile.clear(infile.rdstate() & ~CFile::__fail_flag); // Unset the fail flag.
  }

  void recode_sequence_to_bitcode()
  {
    if (type == SeqType_dna)
    {
      char *it     = data.begin();
      char *it_end = data.end();

      while (it != it_end)
      {
        *it = recode_DNA(*it);
        ++it;
      }
      bitcode_recoded = true;
    }
    else if (type == SeqType_protein)
    {}
  }

  void backrecode_bitcode_sequence()
  {
    if (type == SeqType_dna)
    {
      char *it     = data.begin();
      char *it_end = data.end();

      while (it != it_end)
      {
        *it = backrecode_DNA(*it);
        ++it;
      }
      bitcode_recoded = false;
    }
    else if (type == SeqType_protein)
    {}
  }

  // Crashed in this function on linux - do not know why yet 64 bit problem??
  void check_allowed_symbols_in_sequence(const char *symbols)
  {
    faststring::size_type pos=0;

    pos = data.find_first_not_of(symbols, pos);
    while (pos != faststring::npos)
    {
      std::cerr << "Found illegal symbol in sequence: " << short_name
      << " Position " << pos << ". Found char " << data[pos] << std::endl;

      ++pos;
      pos = data.find_first_not_of(symbols, pos);
    }

  }


  void export_fasta()
  {
    

  }

  void export_phylip()
  {


  }

  void export_nexus()
  {


  }


  void readFastaSequence_generic(CFile& infile, processing_flag pflag)
  {
    char     c   = '\0';
    faststring::size_type  pos = 0;

    // TODO: Prepare to remove gaps and/or stars:

    // Increase the readability:
    bool   convert_toupper_bool  = (pflag & convert_toupper);   //
    bool   convert_ambig2N_bool  = (pflag & convert_ambig2N);   //
    bool   removegaps_bool       = (pflag & removegaps);        //
    bool   removestars_bool      = (pflag & removestars);       //
    bool   gaps_to_ambig_bool    = (pflag & gaps_to_ambig);     //

    full_name = "";
    data.clear();

    c = infile.getchar();
    while (!infile.fail() && isspace(c))
      c = infile.getchar();
    while (!infile.fail() && c == ';')
    {
      infile.ignore('\n');
      c = infile.getchar();
      while (!infile.fail() && isspace(c) )
        c = infile.getchar();
    }

    if (c != '>')
    {
      infile.clear(CFile::__fail_flag | infile.rdstate()); // 20.7.2021: add rdstate, so that eof flags is not overwritten
      return;
    }
    infile.getline(full_name);
    if (infile.fail()) // getline only fails, if 0 characters could be read.
    {
      return;
    }
    find_shortname_and_description();

    c = infile.peekchar();
    while ( !infile.fail() && c != '>' )
    {
      //      std::cerr << data.capacity() << std::endl;
      if (type == SeqType_dna)
        getLineOfDNA_generic(infile, pos, pflag);
      else if (type == SeqType_protein)
        getLineOfProtein_generic(infile, pos, pflag);
      else if (type == SeqType_molecular)
      {
        getRawLine(infile);
        if (convert_toupper_bool)
          toupper_this_sequence();

        if (convert_ambig2N_bool)
          convert_ambig2N_this_sequence();
        if (gaps_to_ambig_bool)
          gaps_to_ambig_this_sequence();
        if (removegaps_bool || removestars_bool)
          remove_gaps_stars_this_sequence(removegaps_bool, removestars_bool);
      }
      else if (type == SeqType_auto_detect_type)
      {
        getRawLine(infile);
        if (convert_toupper_bool)
          toupper_this_sequence();
//        auto_detect_datatype();

        if (convert_ambig2N_bool)
          convert_ambig2N_this_sequence();
        if (gaps_to_ambig_bool)
          gaps_to_ambig_this_sequence();
        if (removegaps_bool || removestars_bool)
          remove_gaps_stars_this_sequence(removegaps_bool, removestars_bool);
      }
      c = infile.peekchar();
    }
    if (data.length() > 0)
      if ( !(infile.fail_reason1() || infile.fail_reason2() ) )
        infile.clear(infile.rdstate() & ~CFile::__fail_flag); // Unset the fail flag.
  }

/*
  const basic_sequence_stats &get_basic_seq_stats()
  {
    return basic_seq_stats;
  }
*/

  bool range_contains_gaps_or_Ns(faststring::size_type pos1, faststring::size_type pos2)
  {
    faststring::size_type len = length();

    if (pos1 >= pos2 || pos1 > len) // Empty range
      return false;
    if (pos2 > len)
      pos2 = len;

    const char *p1, *p2;
    p1 = getSeqBegin() + pos1;
    p2 = getSeqBegin() + pos2;

    char c='\0';

    while (p1 != p2)
    {
      c = *p1;
      if (c == 'N' || c == 'n' || c == '-')
        break;
      ++p1;
    }
    if (c == 'N' || c == 'n' || c == '-')
      return true;
    else
      return false;
  }

  // Requirement: Amino acids have to encoded as  upper case characters.
  // Returns -1 in case the sequence is not proper blosum code.
  // Returns 0 or positive value indicating the number of Js that have been replaced.
  faststring::size_type check_protein_sequence_for_blosum_compatibility(faststring::size_type &refpos)
  {
    unsigned pos;  // We would need a new version of basic-DNA-RNA-AA-routines.h if we change this type.

    faststring::size_type  count = data.replace_char_count('J','X');
    bool OK    = is_aa_blosum_code_or_gap_sequence(data.begin(), data.end(), &pos);

    if (!OK)
    {
      refpos = pos;
      return -1;
    }

    return count;
  }

  void determine_indices_for_pattern(faststring pattern, std::vector<size_t> &indices)
  {
    size_t pos = data.find(pattern, 0);
    while(pos != faststring::npos)
    {
      indices.push_back(pos);
      pos = data.find(pattern,pos+1);
    }
  }

  // void readFastaSequence(CFile& infile) is deprecated:
  // Use one of the other read functions below.
  // This function used to by of type
  // toupper_ambig2N_removegaps

  /*   void readFastaSequence(CFile& infile) */
  /*   { */
  /*     #warning "The readFastaSequence function is deprecated." */

  /*     char     c   = '\0'; */

  /*     full_name = ""; */
  /*     data.clear(); */

  /*     c = infile.getchar(); */
  /*     while (!infile.fail() && isspace(c)) */
  /*       c = infile.getchar(); */
  /*     while (!infile.fail() && c == ';') */
  /*     { */
  /*       infile.ignore('\n'); */
  /*       c = infile.getchar(); */
  /*       while (!infile.fail() && isspace(c) ) */
  /* 	c = infile.getchar(); */
  /*     } */

  /*     if (c != '>') */
  /*     { */
  /*       infile.clear(CFile::__fail_flag); */
  /*       return; */
  /*     } */
  /*     infile.getline(full_name); */
  /*     if (infile.fail()) */
  /*     { */
  /*       return; */
  /*     } */
  /*     find_shortname_and_description(); */

  /*     c = infile.peekchar(); */
  /*     while ( !infile.fail() && c != '>' ) */
  /*     { */
  /*       //      cerr << data.capacity() << endl; */
  /*       if (type == dna) */
  /* 	getLineOfDNA_toupper_ambig2N_removegaps(infile); */
  /*       else */
  /* 	std::cout << "Not implemented." << std::endl; */
  /*       c = infile.peekchar(); */
  /*     } */
  /*   } */


  //bool less_than_full_seqname_using_pointer(CSequence_Mol *, CSequence_Mol *);
  //bool less_than_full_seqname_case_insensitive_using_pointer(CSequence_Mol *, CSequence_Mol *);
  //bool 


  friend bool less_than_full_seqname_using_pointer(CSequence_Mol *ps1, CSequence_Mol *ps2)
  {
    return ps1->full_name < ps2->full_name;
  }

  friend bool less_than_full_seqname_caseinsensitive_using_pointer(CSequence_Mol *ps1, CSequence_Mol *ps2)
  {
    return less_than_faststring_case_insensitive(ps1->full_name, ps2->full_name);
  }

  friend bool less_than_full_seqname_first_unsigned_using_pointer(CSequence_Mol *ps1, CSequence_Mol *ps2)
  {
    return less_than_faststring_first_unsigned(ps1->full_name, ps2->full_name);
  }

};

#endif
