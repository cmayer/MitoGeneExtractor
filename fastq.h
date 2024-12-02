#ifndef FASTQ_H
#define FASTQ_H

#include <iostream>
#include "faststring3.h"
#include "CFile/CFile2_3.h"

class fastq_record
{
  faststring identifier;
  faststring seq;
  faststring second_identifier;
  faststring qual;

  public:
  // Sets state of infile to __fail_reason1 if the content is not valid
  void read_next_record(CFile &infile, unsigned processing_flag)
  {
    infile.getline(identifier);
    while (identifier.empty() )
    {
      infile.getline(identifier);
      if (infile.eof())
	return;
    }

    infile.getline(seq);
    infile.getline(second_identifier);
    infile.getline(qual);
    if (is_invalid() )
    {
      infile.clear(infile.rdstate() | CFile::__fail_reason1);
    }
    identifier.erase_front(1); // Remove the '@' from the front of the sequence name.

    if (processing_flag == 1)
    {
      seq.ToUpper();
    }
  }

  int is_invalid()
  {
    if (identifier.length() == 0)
      return 1;
    if (identifier[0] != '@')
      return 1;
    if (identifier.length() <= 1)
      return 2;
    if (seq.length() != qual.length())
      return 3;
    return 0;
  }

  void trim(unsigned beg, unsigned end)
  {
    // The indices are the indices of the good range.
    // end is the length of the new sequence after trimming only the back end. 
    seq.shorten(end);
    qual.shorten(end);
    // beg is the number of nucleotides trimmed from the front
    seq.erase_front(beg);
    qual.erase_front(beg);
  }

  void clear_seq()
  {
    seq.clear();
    qual.clear();
  }
  
  unsigned length()
  {
    return seq.length();
  }

  const faststring &get_identifier()
  {
    return identifier;
  }

  const faststring &get_seq()
  {
    return seq;
  }

  const faststring &get_second_identifier()
  {
    return second_identifier;
  }

  const faststring &get_qual()
  {
    return qual;
  }

  const char *get_identifier_cstr()
  {
    return identifier.c_str();
  }

  const char *get_seq_cstr()
  {
    return seq.c_str();
  }

  const char *get_second_identifier_cstr()
  {
    return second_identifier.c_str();
  }

  const char *get_qual_cstr()
  {
    return qual.c_str();
  }

  void print(FILE *out)
  {
    fputs(get_identifier_cstr(), out);
    fputc('\n', out);
    fputs(get_seq_cstr(), out);
    fputc('\n', out);
    fputs(get_second_identifier_cstr(), out);
    fputc('\n', out);
    fputs(get_qual_cstr(), out);
    fputc('\n', out);
  }

  void print(std::ostream &os)
  {
    os << identifier << '\n'
       << seq  << '\n'
       << second_identifier << '\n'
       << qual << '\n';
  }

};


#endif
