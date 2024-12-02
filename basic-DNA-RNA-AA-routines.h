#ifndef DNA_ROUTINES_H
#define DNA_ROUTINES_H

/* Todo:  count DNA bases, guess data type, is AA, etc.  */

// Changes 10/11.11.2016:
// * Removed a bug: (c == 'O' && c == 'U') should read (c == 'O' || c == 'U').
// The impact of this error is low. Basically, we overlook O or U as not allowed amino acids.
// Since they should not occur in normal amino acid sequences, it does not make much of a difference.
// * added the blosum symbols, lists and functions
// * added some upper lists and functions.
// * added a look-up table for the blosum symbols. This is much more efficient.
// * removed z from aa_symbols, aa_symbols_ext and _upper, since it as an ambiguity character.


#include <cstring>
#include <cctype>

static const char DNARNA_symbols[]         = "aAcCgGtTuU";
static const char DNARNA_symbols_upper[]   = "ACGTU";
static const char aa_symbols[]             = "aArRnNdDcCeEqQgGhHiIlLkKmMfFpPsStTwWyYvVxX";
static const char aa_symbols_upper[]       = "ARNDCEQGHILKMFPSTWYVX";

// Extended also includes: O and U:
static const char aa_symbols_ext[]         = "aArRnNdDcCeEqQgGhHiIlLkKmMfFpPsStTwWyYvVxXuUoO*";
static const char aa_symbols_ext_upper []  = "ARNDCEQGHILKMFPSTWYVXUO*";

static const char non_aa_in_ABC [] = "BJOUZbjouz";

// AA ambiguity characters: J is a very rare ambiguity character also included in this set.
// J is not allowed in many programs.
static const char aa_ambig [] = "bBzZjJxX";

// Usually a safe list of allowed characters are those in the blosum matrix:
// Lists the symbols that are allowed in the blosum matrix: Contains B,Z but not J as ambiguity character.
// See also the look up tables:  is_aa_blosum_code[] and is_aa_blosum_code_or_gap[]

static const char aa_symbols_blosum[]         = "aArRnNdDcCeEqQzZgGhHiIlLkKmMfFpPsStTwWyYvVbBxX*";
static const char aa_symbols_blosum_upper []  = "ARNDCEQZGHILKMFPSTWYVBX*";

static const char non_aa_symbols_blosum[]         = "jJoOuU";
static const char non_aa_symbols_blosum_upper []  = "JOU";

static const char DNARNA_iupac_symbols []          = "aAcCgGtTuURYSWryswKMBDkmbdHVNhvn";
static const char DNARNA_iupac_symbols_upper []    = "ACGTURYSWKMBDHVN";

static const char phylip_DNARNA_allowed_symbols [] = "?-nNaAcCgGtTuURYSWryswKMBDkmbdHVhvoOxX";

static const char phylip_neither_aa_dna_symbol [] = "jJzZ";
static const char allowed_non_ABC [] = "-?*";

// O and U are special amino acids that are usually not in the models.
// Chars that are not an aa symbol: B, J, O, U, Z
// X is the amibuity code

// Includes the symbols which occur in the blosum matrices: 
const int is_aa_blosum_code_lookup[] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,1,0,0,0,0,0,0,0, // 40-49  45: - is not included, 42: * included 
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,0,0,1,1,1,1,1, // 60-69  63: ? is not included
  1,1,1,1,0,1,1,1,1,0, // 70-79
  1,1,1,1,1,0,1,1,1,1, // 80-89
  1,0,0,0,0,0,0,1,1,1, // 90-99
  1,1,1,1,1,1,0,1,1,1, //100-109
  1,0,1,1,1,1,1,0,1,1, //110-119
  1,1,1,0,0,0,0,0 //120-127
};

// Includes the symbols which occur in the blosum matrices and the gap character:
const int is_aa_blosum_code_or_gap_lookup[] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,1,0,0,1,0,0,0,0, // 40-49  45: - is included, 42: * included 
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,0,0,1,1,1,1,1, // 60-69  63: ? is not included
  1,1,1,1,0,1,1,1,1,0, // 70-79
  1,1,1,1,1,0,1,1,1,1, // 80-89
  1,0,0,0,0,0,0,1,1,1, // 90-99
  1,1,1,1,1,1,0,1,1,1, //100-109
  1,0,1,1,1,1,1,0,1,1, //110-119
  1,1,1,0,0,0,0,0 //120-127
};

// Inlcudes all DNA symbols, no ambiguity characters
const char is_DNA_lookup[128] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,0,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,0,0,1,0,1,0,0, // 60-69     63:'?', 65:'A' to  69:'E'
  0,1,0,0,0,0,0,0,0,0, // 70-79
  0,0,0,0,1,0,0,0,0,0, // 80-89
  0,0,0,0,0,0,0,1,0,1, // 90-99
  0,0,0,1,0,0,0,0,0,0, //100-109
  0,0,0,0,0,0,1,0,0,0, //110-119
  0,0,0,0,0,0,0,0      //120-127
};

// Inlcudes all DNA symbols and gap '-', no ambiguity characters
const char is_DNA_or_GAP_lookup[128] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,1,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,0,0,1,0,1,0,0, // 60-69     63:'?', 65:'A' to  69:'E'
  0,1,0,0,0,0,0,0,0,0, // 70-79
  0,0,0,0,1,0,0,0,0,0, // 80-89
  0,0,0,0,0,0,0,1,0,1, // 90-99
  0,0,0,1,0,0,0,0,0,0, //100-109
  0,0,0,0,0,0,1,0,0,0, //110-119
  0,0,0,0,0,0,0,0      //120-127
};

// Includes all ambiguity characters as well as ?
const char is_DNA_or_AMBIG_lookup[128] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,0,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,1,0,1,1,1,1,0, // 60-69     63:'?', 65:'A' to  69:'E' .,.,.,?,.,A,B,C,D,E,
  0,1,1,0,0,1,0,1,1,0, // 70-79                               F,G,H,I,J,K,L,M,N,O,
  0,0,1,1,1,0,1,1,0,1, // 80-89     U is not included         P,Q,R,S,T,U,V,W,X,Y,
  0,0,0,0,0,0,0,1,1,1, // 90-99                               Z,.,.,.,.,.,.,a,b,c,
  1,0,0,1,1,0,0,1,0,1, //100-109                              d,e,f,g,h,i,j,k,l,m,
  1,0,0,0,1,1,1,0,1,1, //110-119                              n,o,p,q,r,s,t,u,b,w,
  0,1,0,0,0,0,0,0      //120-127                              x,y,z,
};

// Inlcudes all ambiguity characters as well as ? and -
const char is_DNA_or_AMBIG_or_GAP_lookup[128] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,1,0,0,0,0, // 40-49     45: - is included
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,1,0,1,1,1,1,0, // 60-69     63:'?', 65:'A' to  69:'E' .,.,.,?,.,A,B,C,D,E,
  0,1,1,0,0,1,0,1,1,0, // 70-79                               F,G,H,I,J,K,L,M,N,O,
  0,0,1,1,1,0,1,1,0,1, // 80-89     U is not included         P,Q,R,S,T,U,V,W,X,Y,
  0,0,0,0,0,0,0,1,1,1, // 90-99                               Z,.,.,.,.,.,.,a,b,c,
  1,0,0,1,1,0,0,1,0,1, //100-109                              d,e,f,g,h,i,j,k,l,m,
  1,0,0,0,1,1,1,0,1,1, //110-119                              n,o,p,q,r,s,t,u,b,w,
  0,1,0,0,0,0,0,0      //120-127                              x,y,z,
};


// Inlcudes all ambiguity characters as well as ?. Does not include DNA/RNA characters.
const char is_DNA_AMBIG_lookup[128] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,0,0,0,0,0, // 40-49     45: - is not included
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,1,0,0,1,0,1,0, // 60-69     63:'?', 65:'A' to  69:'E' .,.,.,?,.,A,B,C,D,E,
  0,0,1,0,0,1,0,1,1,0, // 70-79                               F,G,H,I,J,K,L,M,N,O,
  0,0,1,1,0,0,1,1,0,1, // 80-89     U is not included         P,Q,R,S,T,U,V,W,X,Y,
  0,0,0,0,0,0,0,0,1,0, // 90-99                               Z,.,.,.,.,.,.,a,b,c,
  1,0,0,0,1,0,0,1,0,1, //100-109                              d,e,f,g,h,i,j,k,l,m,
  1,0,0,0,1,1,0,0,1,1, //110-119                              n,o,p,q,r,s,t,u,v,w,
  0,1,0,0,0,0,0,0      //120-127                              x,y,z,
};

// Inlcudes all ambiguity characters as well as '?' and '-'. Does not include DNA/RNA characters.
const char is_DNA_AMBIG_OR_gap_lookup[128] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,1,0,0,0,0, // 40-49     45: - is not included
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,1,0,0,1,0,1,0, // 60-69     63:'?', 65:'A' to  69:'E' .,.,.,?,.,A,B,C,D,E,
  0,0,1,0,0,1,0,1,1,0, // 70-79                               F,G,H,I,J,K,L,M,N,O,
  0,0,1,1,0,0,1,1,0,1, // 80-89     U is not included         P,Q,R,S,T,U,V,W,X,Y,
  0,0,0,0,0,0,0,0,1,0, // 90-99                               Z,.,.,.,.,.,.,a,b,c,
  1,0,0,0,1,0,0,1,0,1, //100-109                              d,e,f,g,h,i,j,k,l,m,
  1,0,0,0,1,1,0,0,1,1, //110-119                              n,o,p,q,r,s,t,u,v,w,
  0,1,0,0,0,0,0,0      //120-127                              x,y,z,
};


// Inlcudes all non R/Y ambiguity characters as well as ?. Does not include DNA/RNA characters.
const char is_DNA_AMBIG_non_R_Y_lookup[128] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,0,0,0,0,0, // 40-49     45: - is not included
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,1,0,0,1,0,1,0, // 60-69     63:'?', 65:'A' to  69:'E' .,.,.,?,.,A,B,C,D,E,
  0,0,1,0,0,1,0,1,1,0, // 70-79                               F,G,H,I,J,K,L,M,N,O,
  0,0,0,1,0,0,1,1,0,0, // 80-89     U is not included         P,Q,R,S,T,U,V,W,X,Y,
  0,0,0,0,0,0,0,0,1,0, // 90-99                               Z,.,.,.,.,.,.,a,b,c,
  1,0,0,0,1,0,0,1,0,1, //100-109                              d,e,f,g,h,i,j,k,l,m,
  1,0,0,0,0,1,0,0,1,1, //110-119                              n,o,p,q,r,s,t,u,v,w,
  0,0,0,0,0,0,0,0      //120-127                              x,y,z,
};


/* inline bool is_DNA_base_upper(char c) */
/* {  */
/*   if (c == 'A' || c == 'C' || c == 'G' || c == 'T' ) */
/*     return true; */
/*   return false; */
/* } */

/* inline bool is_DNA_base_lower(char c) */
/* {  */
/*   if (c == 'a' || c == 'c' || c == 'g' || c == 't' ) */
/*     return true; */
/*   return false; */
/* } */

// Only defined if c <= 127
inline bool is_DNA_base(char c)
{
  return is_DNA_lookup[(int)c];
}

inline bool is_DNA_base_or_GAP(char c)
{
  return is_DNA_or_GAP_lookup[(int)c];
}

inline bool is_DNA_or_DNA_ambig(char c)
{
  return is_DNA_or_AMBIG_lookup[(int)c];
}

inline bool is_DNA_or_DNA_ambig_or_GAP(char c)
{
  return is_DNA_or_AMBIG_or_GAP_lookup[(int)c];
}

inline bool is_DNA_or_RNA_base(char c)
{
  if (c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
      c == 'a' || c == 'c' || c == 'g' || c == 't' || c == 'U' || c == 'u')
    return true;
  return false;
}

inline bool is_allowed_non_ABC(char c)
{
  if (c == '-' || c == '?' || c == '*')
    return true;
  else
    return false;
}

inline bool is_phylip_aa_dna_symbol(char c) // Does not include -? or similar
{
  // Only one character is not allowed:
  if (c=='o' || c == 'O')
    return false;
  if (c >='a' && c <= 'z')
  {
      return true;
  }
  else if (c >='A' && c <= 'Z')
  {
      return true;
  }
  return false;
}

inline bool is_gap(char c)
{
  return c == '-';
}


// Only defined if c <= 127
// Includes only DNA ambig symbols. Does not include DNA symbols. See also is_DNA_or_AMBIG below.
inline bool is_DNA_ambig(char c)
{
  return is_DNA_AMBIG_lookup[(int)c];
/*   if (c == 'N'  ||  c == '?'  || c == 'n' || */
/*       c == 'R'  ||  c == 'Y'  || c == 'S' || c == 'W' || */
/*       c == 'r'  ||  c == 'y'  || c == 's' || c == 'w' || */
/*       c == 'K'  ||  c == 'M'  || c == 'B' || c == 'D' || */
/*       c == 'k'  ||  c == 'm'  || c == 'b' || c == 'd' || */
/*       c == 'H'  ||  c == 'V'  || */
/*       c == 'h'  ||  c == 'v') */
/*     return true; */
/*   else */
/*     return false; */
}

// Only defined if c <= 127
// Includes only DNA ambig symbols and '-'. Does not include DNA symbols. 
// See also is_DNA_or_AMBIG below.
inline bool is_DNA_ambig_OR_gap(char c)
{
  return is_DNA_AMBIG_OR_gap_lookup[(int)c];
}



/* inline bool is_DNA_iupac_ambig_upper(char c) */
/* { */
/*   if (c == 'N'  ||  c == '?'  || */
/*       c == 'R'  ||  c == 'Y'  || c == 'S' || c == 'W' || */
/*       c == 'K'  ||  c == 'M'  || c == 'B' || c == 'D' || */
/*       c == 'H'  ||  c == 'V') */
/*     return true; */
/*   else */
/*     return false; */
/* } */

inline bool is_DNA_ambig_non_RY(char c)
{
  return is_DNA_AMBIG_non_R_Y_lookup[(int)c];
}








// A list can be found here: https://en.wikipedia.org/wiki/Amino_acid
// The acceptance of J as an ambiguous amino acid symbol is low.
// J does not occur in the blosm matrix and is not allowed in some programs.
inline bool is_aa_ambig(char c)
{
  if (c == 'b' || c == 'B' || c == '?' || c == 'x' || c == 'X' || 
      c == 'j' || c == 'J' || c == 'z' || c == 'Z')
    return true;
  else
    return false;
}

inline bool is_aa_blosum_code(unsigned char c)
{
  return (is_aa_blosum_code_lookup[c] != 0);
}

inline bool is_aa_blosum_code_or_gap(unsigned char c)
{
  return (is_aa_blosum_code_or_gap_lookup[c] != 0);
}



inline int is_aa_upper(char c)
{
  if (c >='A' && c <= 'Z')
  {
    // AA ambigs:
    if (c == 'B' || c == 'X' || c == 'J' || c == 'Z')
      return false;

    // The few symbols that represent nothing.
    // Indeed they are listed as the 21st and 22nd amino acid in https://en.wikipedia.org/wiki/Amino_acid
    // But they are usually not handled by programs, so we return false.
    if (c == 'O' || c == 'U') // OU
      return false;
    return true;
  }
  return false;
}


inline int is_aa_or_aa_ambig(char c) // includes '?' and  B, Z, J, X
{
  if (c >='a' && c <= 'z')
  {
    if (c != 'o' && c != 'u') // Only not allowed characters in a..z are: ou
      return true;
  }
  else if (c >='A' && c <= 'Z')
  {
    if (c != 'O' && c != 'U') // Only not allowed characters in a..z are: OU
      return true;
  }

  if (c == '?' || c == '*')
    return true;

  return false;
  //  return (strchr(aa_symbols, c) != NULL);
}

inline int is_aa_or_aa_ambig_extended(char c) // O and U are allowed  // includes '?' and  B, Z, J, X
{
  if (c >='a' && c <= 'z')
  {
      return true;
  }
  else if (c >='A' && c <= 'Z')
  {
      return true;
  }

  if (c == '?' || c == '*')
    return true;

  return false;
  //  return (strchr(aa_symbols, c) != NULL);
}


inline char complementBase_of_upper(char c)
{
  switch (c) {
  case 'A':  c = 'T'; break;
  case 'T':  c = 'A'; break;
  case 'G':  c = 'C'; break;
  case 'C':  c = 'G'; break;
    
  case 'R':  c = 'Y'; break;
  case 'Y':  c = 'R'; break;

  case 'M':  c = 'K'; break;
  case 'K':  c = 'M'; break;

  case 'B':  c = 'V'; break;
  case 'V':  c = 'B'; break;

  case 'D':  c = 'H'; break;
  case 'H':  c = 'D'; break;

  case 'N':  c = 'N'; break;
  case '?':  c = '?'; break;
  case '-':  c = '-'; break;

  case 'S':  /*c = 'S';*/  break;
  case 'W':  /*c = 'W';*/  break;

  default:
    c = '*';
    break;
  }
  return c;
}

inline char complementBase_of_lower(char c)
{
  switch (c) {
  case 'a':  c = 't'; break;
  case 't':  c = 'a'; break;
  case 'g':  c = 'c'; break;
  case 'c':  c = 'g'; break;
    
  case 'r':  c = 'y'; break;
  case 'y':  c = 'r'; break;

  case 'm':  c = 'k'; break;
  case 'k':  c = 'm'; break;

  case 'b':  c = 'v'; break;
  case 'v':  c = 'b'; break;

  case 'd':  c = 'h'; break;
  case 'h':  c = 'd'; break;

  case 'n':  c = 'n'; break;
  case '?':  c = '?'; break;
  case '-':  c = '-'; break;

  case 's':  /*c = 's';*/  break;
  case 'w':  /*c = 'w';*/  break;

  default:
    c = '*';
    break;
  }
  return c;
}

inline char complementBase(char c)
{
  if ( isupper(c) )
    return complementBase_of_upper(c);
  else
    return complementBase_of_lower(c);
}


inline void CG_AT_content_in_region(const char *beg, const char *end, unsigned &AT, unsigned &CG)
{
  AT = 0;
  CG = 0;

  char c;

  while (beg != end)
  {
    c = *beg;
    if (c == 'C' || c == 'c' || c == 'G' || c == 'g')
      ++CG;
    else if (c == 'A' || c == 'a' || c == 'T' || c == 't')
      ++AT;

    ++beg;
  }
}

inline void content_in_region_DNA(const char *beg, const char *end, unsigned &AT, unsigned &CG, unsigned &gap, unsigned &ambig, unsigned &unknown)
{
  AT      = 0;
  CG      = 0;
  gap     = 0;
  ambig   = 0;
  unknown = 0;

  char c;

  while (beg != end)
  {
    c = *beg;
    if (c == 'C' || c == 'c' || c == 'G' || c == 'g')
      ++CG;
    else if (c == 'A' || c == 'a' || c == 'T' || c == 't')
      ++AT;
    else if (c == '-')
      ++gap;
    else if (is_DNA_ambig(c) )
      ++ambig;
    else
      ++unknown;

    ++beg;
  }
}

inline bool is_aa_blosum_code_or_gap_sequence(const char *p1, const char *p2, unsigned *offending_pos)
{
  const char *startp = p1;

  while (p1 != p2)
  {
    if (!is_aa_blosum_code_or_gap(*p1) )
    {
      *offending_pos = p1 - startp;
      return false;
    }
    ++p1;
  }
  return true;
}



// Recoding symbols:
// A char can hold 8 bits, enough for 8 characters to be "in or out",
// to be switched on or off. Thus a char can be used to hold
// all bases including ambiguities as follows:
// Say A, C, G, T are the bits 1,2,3,4,
// the ambiguity {A,C} is 1+2 = 3

// A:       65, 97  -> 1
// C:       67, 99  -> 2
// G:       71,103  -> 4
// T:       84,116  -> 8
// R(AG):   82,114  -> 5
// Y(CT):   89,121  -> 10
// S(CG):   83,115  -> 6
// W(TA):   87,119  -> 9
// K(TG):   75,107  -> 12
// M(CA):   77,109  -> 3
// B(CTG):  66,98   -> 14
// D(ATG):  68,100  -> 13
// H(ATC):  72,104  -> 11
// V(ACG):  86,118  -> 7
// N(ACGT): 78,110  -> 15
// -        45      -> 127 !!!!
// ?        63      -> 15

const char ascii2recode_DNA[128] = {
  0,0,0,0,0,0,0,0,0,0,    //  0- 9
  0,0,0,0,0,0,0,0,0,0,    // 10-19
  0,0,0,0,0,0,0,0,0,0,    // 20-29
  0,0,0,0,0,0,0,0,0,0,    // 30-39
  0,0,0,0,0,127,0,0,0,0,  // 40-49
  0,0,0,0,0,0,0,0,0,0,    // 50-59
  0,0,0,15,0,1,14,2,13,0, // 60-69
  0,4,11,0,0,12,0,3,15,0, // 70-79
  0,0,5,6,8,0,7,9,0,10,   // 80-89
  0,0,0,0,0,0,0,1,14,2,   // 90-99
  13,0,0,4,11,0,0,0,0,0,  //100-109
  0,0,0,0,5,0,8,0,7,9,    //110-119
  0,10,0,0,0,0,0,0        //120-127
};

const char recode2ascii_DNA[128] = {
  0,65,67,77,71,82,83,86,84,87, //  0- 9
  89,72,75,68,66,78,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0,    // 20-29
  0,0,0,0,0,0,0,0,0,0,    // 30-39
  0,0,0,0,0,0,0,0,0,0,    // 40-49
  0,0,0,0,0,0,0,0,0,0,    // 50-59
  0,0,0,0,0,0,0,0,0,0,    // 60-69
  0,0,0,0,0,0,0,0,0,0,    // 70-79
  0,0,0,0,0,0,0,0,0,0,    // 80-89
  0,0,0,0,0,0,0,0,0,0,    // 90-99
  0,0,0,0,0,0,0,0,0,0,    //100-109
  0,0,0,0,0,0,0,0,0,0,    //110-119
  0,0,0,0,0,0,0,45        //120-127
};

const char recodeISDNAambig[128] = {
  0,0,0,1,0,1,1,1,0,1,    //  0- 9
  1,1,1,1,1,1,0,0,0,0,    // 10-19
  0,0,0,0,0,0,0,0,0,0,    // 20-29
  0,0,0,0,0,0,0,0,0,0,    // 30-39
  0,0,0,0,0,0,0,0,0,0,    // 40-49
  0,0,0,0,0,0,0,0,0,0,    // 50-59
  0,0,0,0,0,0,0,0,0,0,    // 60-69
  0,0,0,0,0,0,0,0,0,0,    // 70-79
  0,0,0,0,0,0,0,0,0,0,    // 80-89
  0,0,0,0,0,0,0,0,0,0,    // 90-99
  0,0,0,0,0,0,0,0,0,0,    //100-109
  0,0,0,0,0,0,0,0,0,0,    //110-119
  0,0,0,0,0,0,0,0         //120-127
};

const char recodeISDNA[128] = {
  0,1,1,0,1,0,0,0,1,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,0,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,0,0,0,0,0,0,0, // 60-69
  0,0,0,0,0,0,0,0,0,0, // 70-79
  0,0,0,0,0,0,0,0,0,0, // 80-89
  0,0,0,0,0,0,0,0,0,0, // 90-99
  0,0,0,0,0,0,0,0,0,0, //100-109
  0,0,0,0,0,0,0,0,0,0, //110-119
  0,0,0,0,0,0,0,0 //120-127
};


// A:       65, 97  -> 0
// C:       67, 99  -> 1
// G:       71,103  -> 2
// T:       84,116  -> 3
// R(AG):   82,114  -> 4
// Y(CT):   89,121  -> 5
// S(CG):   83,115  -> 6
// W(TA):   87,119  -> 7
// K(TG):   75,107  -> 8
// M(CA):   77,109  -> 9
// B(CTG):  66,98   -> 10
// D(ATG):  68,100  -> 11
// H(ATC):  72,104  -> 12
// V(ACG):  86,118  -> 13
// N(ACGT): 78,110  -> 14
// -        45      -> 127
// ?        63      -> 14
// all others:      -> 16

const char ascii2simpleIndex_DNA[128] = {
  16,16,16,16,16, 16,16,16,16,16,     //  0- 9
  16,16,16,16,16, 16,16,16,16,16,     // 10-19
  16,16,16,16,16, 16,16,16,16,16,     // 20-29
  16,16,16,16,16, 16,16,16,16,16,     // 30-39
  16,16,16,16,16,127,16,16,16,16,     // 40-49
  16,16,16,16,16, 16,16,16,16,16,     // 50-59
  16,16,16,14,16,  0,10, 1,11,16,     // 60-69
  16, 2,12,16,16,  8,16, 9,14,16,     // 70-79
  16,16, 4, 6, 3, 16,13, 7,16, 5,     // 80-89
  16,16,16,16,16, 16,16, 0,10, 1,     // 90-99
  11,16,16, 2,12, 16,16, 8,16, 9,     //100-109
  14,16,16,16, 4,  6, 3,16,13, 7,     //110-119
  16, 5,16,16,16, 16,16,16            //120-127
};

const char simpleIndex2ascii_DNA[128] = {
  65,67,71,84,82,89,83,87,75,77,     //  0- 9
  66,68,72,86,78, 0, 0, 0, 0, 0,     // 10-19
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 20-29
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 30-39
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 40-49
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 50-59
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 60-69
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 70-79
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 80-89
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     // 90-99
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     //100-109
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     //110-119
   0, 0, 0, 0, 0, 0, 0,45            //120-127
};



inline unsigned char recode_DNA(unsigned char c)
{
  return ascii2recode_DNA[c];
}

inline unsigned char backrecode_DNA(unsigned char c)
{
  return recode2ascii_DNA[c];
}

// Only true of character is an ambiguity character for DNA in the iupaccode
inline unsigned char recode_is_DNA_ambig(unsigned char c)
{
  return recodeISDNAambig[c];
}

inline unsigned char recode_is_DNA(unsigned char c)
{
  return recodeISDNA[c];
}

// Recodes can store multiple DNA bases in one char. This function determines whether the base is present in the reode.
inline bool is_base_in_recode(unsigned char recode, unsigned char base)
{
  return recode & recode_DNA(base);
}

inline bool is_valid_recode_DNA(unsigned char c)
{
  return recodeISDNA[c];
}

inline bool is_valid_recode_DNA_or_DNAambig(unsigned char c)
{
  return recodeISDNA[c] || recodeISDNAambig[c];
}

inline bool is_recode_gap(unsigned char c)
{
  return (c == 127);
}


inline unsigned char simpleIndex2DNA(unsigned char c)
{
  return simpleIndex2ascii_DNA[c];
}

inline unsigned char DNA2simpleIndex(unsigned char c)
{
  return ascii2simpleIndex_DNA[c];
}

inline bool is_valid_simple_Index_DNA(unsigned char c)
{
  return simpleIndex2ascii_DNA[c]!=0;
}

// A	Ala	Alanine         01   65
// R	Arg	Arginine        02   82
// N	Asn	Asparagine      04   78
// D	Asp	Aspartic acid   08   68
// C	Cys	Cysteine        16   67
// Q	Gln	Glutamine       32   81
// E	Glu	Glutamic acid   64   69
// G	Gly	Glycine        128   71
// H	His	Histidine      256   72
// I	Ile	Isoleucine     512   73
// L	Leu	Leucine       1024   76
// K	Lys	Lysine        2048   75 
// M	Met	Methionine    4096   77
// F	Phe	Phenylalanine 8192   70
// P	Pro	Proline      16384   80
// S	Ser	Serine       32768   83
// T	Thr	Threonine    65536   84
// W	Trp	Tryptophan  131072   87
// Y	Tyr	Tyrosine    262144   89
// V	Val	Valine      524288   86
// -        45      -> 0
// ?        63      -> 1048575
// X        88      -> 1048575

/* const int ascii2recode_aa[128] = { */
/*   0,0,0,0,0,0,0,0,0,0, //  0- 9 */
/*   0,0,0,0,0,0,0,0,0,0, // 10-19 */
/*   0,0,0,0,0,0,0,0,0,0, // 20-29 */
/*   0,0,0,0,0,0,0,0,0,0, // 30-39 */
/*   0,0,0,0,0,64,0,0,0,0, // 40-49 */
/*   0,0,0,0,0,0,0,0,0,0, // 50-59 */
/*   0,0,0,15,0,1,14,2,13,0, // 60-69 */
/*   0,4,11,0,0,12,0,3,15,0, // 70-79 */
/*   0,0,5,6,8,0,7,9,0,10, // 80-89 */
/*   0,0,0,0,0,0,0,1,14,2, // 90-99 */
/*   13,0,0,4,11,0,0,0,0,0, //100-109 */
/*   0,0,0,0,5,0,8,0,7,9, //110-119 */
/*   0,10,0,0,0,0,0,0 //120-127 */
/*}; */

/* Not yet implemented:
const int32_t ascii2recode_AA[128] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,127,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,15,0,1,14,2,13,0, // 60-69
  0,4,11,0,0,12,0,3,15,0, // 70-79
  0,0,5,6,8,0,7,9,0,10, // 80-89
  0,0,0,0,0,0,0,1,14,2, // 90-99
  13,0,0,4,11,0,0,0,0,0, //100-109
  0,0,0,0,5,0,8,0,7,9, //110-119
  0,10,0,0,0,0,0,0 //120-127
};

const char recode2ascii_AA[128] = {
  0,65,67,77,71,82,83,86,84,87, //  0- 9
  89,72,75,68,66,78,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,0,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,0,0,0,0,0,0,0, // 60-69
  0,0,0,0,0,0,0,0,0,0, // 70-79
  0,0,0,0,0,0,0,0,0,0, // 80-89
  0,0,0,0,0,0,0,0,0,0, // 90-99
  0,0,0,0,0,0,0,0,0,0, //100-109
  0,0,0,0,0,0,0,0,0,0, //110-119
  0,0,0,0,0,0,0,45 //120-127
};

const char recodeISAAambig[128] = {
  0,0,0,1,0,1,1,1,0,1, //  0- 9
  1,1,1,1,1,1,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,0,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,0,0,0,0,0,0,0, // 60-69
  0,0,0,0,0,0,0,0,0,0, // 70-79
  0,0,0,0,0,0,0,0,0,0, // 80-89
  0,0,0,0,0,0,0,0,0,0, // 90-99
  0,0,0,0,0,0,0,0,0,0, //100-109
  0,0,0,0,0,0,0,0,0,0, //110-119
  0,0,0,0,0,0,0,0 //120-127
};

const char recodeISAA[128] = {
  0,1,1,0,1,0,0,0,1,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,0,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,0,0,0,0,0,0,0, // 60-69
  0,0,0,0,0,0,0,0,0,0, // 70-79
  0,0,0,0,0,0,0,0,0,0, // 80-89
  0,0,0,0,0,0,0,0,0,0, // 90-99
  0,0,0,0,0,0,0,0,0,0, //100-109
  0,0,0,0,0,0,0,0,0,0, //110-119
  0,0,0,0,0,0,0,0 //120-127
};

*/

//********************************
// AA to simple index:
//********************************

// A	Ala	Alanine         00   65,  97
// R	Arg	Arginine        01   82, 114
// N	Asn	Asparagine      02   78, 110
// D	Asp	Aspartic acid   03   68, 100
// C	Cys	Cysteine        04   67,  99
// Q	Gln	Glutamine       05   81, 113
// E	Glu	Glutamic acid   06   69, 101
// G	Gly	Glycine         07   71, 103
// H	His	Histidine       08   72, 104
// I	Ile	Isoleucine      09   73, 105
// L	Leu	Leucine         10   76, 108
// K	Lys	Lysine          11   75, 107
// M	Met	Methionine      12   77, 109
// F	Phe	Phenylalanine   13   70, 102
// P	Pro	Proline         14   80, 112
// S	Ser	Serine          15   83, 115
// T	Thr	Threonine       16   84, 116
// W	Trp	Tryptophan      17   87, 119
// Y	Tyr	Tyrosine        18   89, 121
// V	Val	Valine          19   86, 118
// -        45      -> 32
// ?        63      -> 64
// X        88      -> 64, 120

const char ascii2simpleIndex_AA[128] = {
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127,  //  0-  9
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127,  // 10- 19
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127,  // 20- 29
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127,  // 30- 39
    127, 127, 127, 127, 127,  32, 127, 127, 127, 127,  // 40- 49
    127, 127, 127, 127, 127, 127, 127, 127, 127, 127,  // 50- 59
    127, 127, 127,  64, 127,   0, 127,   4,   3,   6,  // 60- 69
     13,   7,   8,   9, 127,  11,  10,  12,   2, 127,  // 70- 79
     14,   5,   1,  15,  16, 127,  19,  17,  64,  18,  // 80- 89
    127, 127, 127, 127, 127, 127, 127,   0, 127,   4,  // 90- 99
      3,   6,  13,   7,   8,   9, 127,  11,  10,  12,  //100-109
      2, 127,  14,   5,   1,  15,  16, 127,  19,  17,  //110-119
     64,  18, 127, 127, 127, 127, 127, 127,            //120-129
};


const char simpleIndex2ascii_AA[128] = {
   'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',     //  0- 9
   'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',     // 10-19
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     // 20-29
     0,   0, '-',   0,   0,   0,   0,   0,   0,   0,     // 30-39
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     // 40-49
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     // 50-59
     0,   0,   0,   0, 'X',   0,   0,   0,   0,   0,     // 60-69
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     // 70-79
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     // 80-89
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     // 90-99
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     //100-109
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     //110-119
     0,   0,   0,   0,   0,   0,   0,   0                //120-127
};

inline unsigned char simpleIndex2AA(unsigned char c)
{
  return simpleIndex2ascii_AA[c];
}

inline unsigned char AA2simpleIndex(unsigned char c)
{
  return ascii2simpleIndex_AA[c];
}

inline bool is_valid_simple_Index_AA(unsigned char c)
{
  return simpleIndex2ascii_AA[c]!=0;
}

#endif

