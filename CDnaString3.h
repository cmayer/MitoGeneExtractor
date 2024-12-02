#ifndef DNASTRING_H
#define DNASTRING_H

#include "primefactors.h"
#include "faststring3.h"
#include <cctype>
#include "basic-DNA-RNA-AA-routines.h"


// Todo:
//       - Check for correct nucl. in push_back, set routine

inline bool isReducibleUnit_global(const char *beg, const char* end)
{
  unsigned       len = end-beg;
  const char     *test_pos;

  int            factor_index;
  unsigned       factor;
  unsigned       factor_len;
  unsigned       i;
  bool           composite;

  factor_index = 0;
  while( (factor = primeFactors[len][factor_index]) != 0 )
  {
    composite = true;
    factor_len = len/factor;
    test_pos = beg;
    for (i=1; i<factor; ++i, test_pos += factor_len)
      if ( strncmp(test_pos, (test_pos + factor_len), factor_len) != 0)
	composite = false;

    if (composite)
      return true;
    ++factor_index;
  }
  return false;
}



class CDnaString : public faststring
{
 public:
/*   using faststring::c_str; */
/*   using faststring::reserve; */
/*   using faststring::size; */
/*   using faststring::capacity; */
/*   using faststring::clear; */
/*   using faststring::check_pos; */
/*   using faststring::set_unckecked; */
/*   using faststring::assign; */
/*   using faststring::push_back; */
/*   using faststring::set_unchecked; */
/*   using faststring::begin; */
/*   using faststring::end; */
/*   using faststring::tolower; */
/*   using faststring::toupper; */


  CDnaString():faststring()
  {}

  CDnaString(const char *s):faststring(s)
  {}

  CDnaString(const faststring &s):faststring(s)
  {}

  void toupper_assign(const char *str)
  {
    faststring::assign(str);
    toupper();
  }

  void toupper_assign(const char *str_begin, const char *str_end)
  {
    faststring::assign(str_begin, str_end);
    toupper();
  }

  void toupper_push_back(char c)
  {
    faststring::push_back(std::toupper(c));
  }

  void toupper_set_unckecked(unsigned pos, char c)
  {
    faststring::set_unckecked(pos, std::toupper(c));
  }

  bool toupper_set(unsigned pos, char c)
  {
    return faststring::set(pos, std::toupper(c));
  }

  unsigned countUnknownBases()
  {
    char *b = begin();
    char *e = end();
    char c;
    unsigned errors = 0;

    while (b != e)
    {
      c = *b;
      if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
	  c != 'R' && c != 'Y' &&
	  c != 'M' && c != 'K' &&
	  c != 'B' && c != 'V' &&
	  c != 'D' && c != 'H' &&
	  c != 'N' )
	++errors;
      ++b;
    }
    return errors;
  }


  unsigned setToComplementOf(const CDnaString &a)
  {
    unsigned   len      = a.size();
    const char *pos     = a.begin();
    const char *pos_end = a.end();
    unsigned   error    = 0;

    reserve(len);
    clear();

    while( pos != pos_end )
    {
      switch (*pos)
      {
        case 'A':  faststring::push_back('T'); break;
        case 'T':  faststring::push_back('A'); break;
        case 'G':  faststring::push_back('C'); break;
        case 'C':  faststring::push_back('G'); break;

        case 'R':  faststring::push_back('Y'); break;
        case 'Y':  faststring::push_back('R'); break;

        case 'M':  faststring::push_back('K'); break;
        case 'K':  faststring::push_back('M'); break;

        case 'B':  faststring::push_back('V'); break;
        case 'V':  faststring::push_back('B'); break;

        case 'D':  faststring::push_back('H'); break;
        case 'H':  faststring::push_back('D'); break;

        case 'N':  faststring::push_back('N'); break;
        case '-':  faststring::push_back('-'); break;

        case 'a':  faststring::push_back('t'); break;
        case 't':  faststring::push_back('a'); break;
        case 'g':  faststring::push_back('c'); break;
        case 'c':  faststring::push_back('g'); break;

        case 'r':  faststring::push_back('y'); break;
        case 'y':  faststring::push_back('r'); break;

        case 'm':  faststring::push_back('k'); break;
        case 'k':  faststring::push_back('m'); break;

        case 'b':  faststring::push_back('v'); break;
        case 'v':  faststring::push_back('b'); break;

        case 'd':  faststring::push_back('h'); break;
        case 'h':  faststring::push_back('d'); break;

        case 'n':  faststring::push_back('n'); break;

        default:
	  ++error;
	  push_back('?');
        break;
      };
      ++pos;
    }
    return error;
  }

  unsigned setToReverseComplementOf(const CDnaString &a)
  {
    unsigned   len      = a.size();
    const char *pos     = a.rbegin();
    const char *pos_end = a.rend();
    unsigned   error    = 0;

    reserve(len);
    clear();

    while( pos != pos_end )
    {
      switch (*pos)
      {
        case 'A':  faststring::push_back('T'); break;
        case 'T':  faststring::push_back('A'); break;
        case 'G':  faststring::push_back('C'); break;
        case 'C':  faststring::push_back('G'); break;

        case 'R':  faststring::push_back('Y'); break;
        case 'Y':  faststring::push_back('R'); break;

        case 'M':  faststring::push_back('K'); break;
        case 'K':  faststring::push_back('M'); break;

        case 'B':  faststring::push_back('V'); break;
        case 'V':  faststring::push_back('B'); break;

        case 'D':  faststring::push_back('H'); break;
        case 'H':  faststring::push_back('D'); break;

        case 'N':  faststring::push_back('N'); break;
        case '-':  faststring::push_back('-'); break;
        default:
	  ++error;
	  push_back('?');
        break;
      };
      --pos;
    }
    return error;
  }


  void setToReverseOf(const CDnaString &a)
  {
    size_type   len      = a.size();
    const char *pos     = a.rbegin();
    const char *pos_end = a.rend();

    reserve(len);
    clear();

    while( pos != pos_end )
    {
      faststring::push_back(*pos);
      --pos;
    }
  }

  // Returns true if a subpattern has been found,
  // false if no subpattern has been found. 
  bool setToMinAlphaRepModCyclicPermOf(const CDnaString &s)
  {
    size_type      len   = s.size();
    char           *tmp  = new char [2*len+1];
    const char     *orig = const_cast<CDnaString&>(s).c_str();

    char           *pos  = tmp;
    char           *pos_end = pos + len;
    //unsigned       i;
    char           *min_str;
    int             strncmp_result;

    bool            found_subpattern = false;

    strcpy(tmp,     orig);
    strcpy(tmp+len, orig);

    min_str = tmp;
    ++pos;
    for (; pos != pos_end; ++pos)
    {
      strncmp_result = strncmp(pos, min_str, len);
      if ( strncmp_result < 0 )
	min_str = pos;
      else if ( strncmp_result == 0 )
      {
	found_subpattern = true;
      }
    }
    *(min_str+len) = '\0';     /* At this position within tmp we have the end of the unit */
    faststring::assign(min_str);
    delete [] tmp;
    return found_subpattern;
  }

  // Returns true if a subpattern has been found,
  // false if no subpattern has been found. 
  bool setToMinAlphaRepModRevCompAndCyclicPermOf(const CDnaString &s)
  {
    CDnaString tmp;
    bool       found_subpattern;

    tmp.setToReverseComplementOf(s);
    found_subpattern = setToMinAlphaRepModCyclicPermOf(tmp);   // *this contains minimum representation of rev comp of s

    tmp.setToMinAlphaRepModCyclicPermOf(s); // tmp contains minimum representation of s itself

    if ( strcmp(c_str(), tmp.c_str() ) > 0 )
    {
      assign( tmp.c_str() );
    }
    return found_subpattern;
  }


  bool isCompositeUnitLength(int x)
  {
    return primeFactors[x][0];
  }

  // Old and inefficient:
/*   bool isReducibleUnit_old() */
/*   { */
/*     unsigned       len   = size(); */
/*     char           *tmp  = new char [2*len+1]; */
/*     const char     *orig = c_str(); */

/*     int            factor_index; */
/*     unsigned       shift; */

/*     strcpy(tmp,     orig); */
/*     strcpy(tmp+len, orig); */

/*     factor_index = 0; */
/*     while( (shift = primeFactors[len][factor_index]) != 0 ) */
/*     { */
/*       if ( strncmp(tmp, tmp+shift, len) == 0) */
/* 	return true; */
/*       ++factor_index; */
/*     } */
/*     return false; */
/*   } */

  // New and works but is redundent - equivalent to current version
/*   bool isReducibleUnit_ok() */
/*   { */
/*     unsigned       len   = size(); */
/*     const char     *orig = c_str(); */

/*     int            factor_index; */
/*     unsigned       factor; */
/*     unsigned       factor_len; */
/*     unsigned       i; */
/*     bool           composite; */

/*     factor_index = 0; */
/*     while( (factor = primeFactors[len][factor_index]) != 0 ) */
/*     { */
/*       composite = true; */
/*       factor_len = len/factor; */
/*       orig = c_str(); */
/*       for (i=1; i<factor; ++i, orig += factor_len) */
/* 	if ( strncmp(orig, orig+factor_len, factor_len) != 0) */
/* 	  composite = false; */

/*       if (composite) */
/* 	return true; */
/*       ++factor_index; */
/*     } */
/*     return false; */
/*   } */

  bool isReducibleUnit()
  {
    return isReducibleUnit_global(c_str(), c_str()+size());
  }


  friend bool operator<(const CDnaString &a, const CDnaString &b)
  {
    // Possibly more efficient would be a strncmp. If equally long, take the shorter one.
    // Not tested yet - thats why its not used!

    //    unsigned minlen = (a.length() < b.length() ? a.length() : b.length());
    //   short    test   =  strncmp(a.begin(), b.begin(), minlen);

    //    if (test == 0)
    //     return a.length() < b.length();
    //   else
    //      return (test < 0);

    // Hi compiler - how do you like this: Haahaahaahaaaaaaaa
    return (strcmp(const_cast<CDnaString&>(a).c_str(), const_cast<CDnaString&>(b).c_str() ) < 0);
  }

};




#endif
