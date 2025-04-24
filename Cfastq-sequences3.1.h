#ifndef  CFASTQ_SEQUENCES_H
#define  CFASTQ_SEQUENCES_H

#include <iostream>
#include <vector>
#include "CSequence_Mol3.1.h"
#include "CSequences3.1.h"
#include "CFile/CFile2_3.h"
#include <list>
#include "fastq.h"
#include "faststring3.h"


inline short load_fastq_save_in_CSequences(CFile &infile, const char *filename,  CSequences3 *pseqs, unsigned processing_flag = 1)
{
  short status_error_flag = 0;
  fastq_record *tmp;

  while(status_error_flag == 0)
  {
    tmp = new fastq_record;
    tmp->read_next_record(infile, processing_flag);

    if (infile.fail_reason1())
    {
      std::cerr << "ERROR: Found malformed fastq record starting in line "
      << infile.line()-4 << " of the input file "
      << filename  << std::endl;
      tmp->print(std::cerr);
      return -1;
    }
    else if (infile.fail() ) // This includes eof.
      status_error_flag = 1;
    else
    {
      pseqs->add_seq_to_dataset(SeqType_dna, tmp->get_identifier(), tmp->get_seq(), 'N');
    }
    delete tmp;
  }
  return 0;
}

class fastq_sequences
{
  std::list<fastq_record *> list_of_sequences;

public:
  fastq_sequences()
  {}

  // This function wants to have the opened infile as well as the file name.
  // The file name is used to print an error message. In the future this could be handled nicer
  // e.g. by using the return value and let the caller print the error message.

  short read_fastq(CFile &infile, const char *filename,  unsigned processing_flag = 1)
  {
    short status_error_flag = 0;
    fastq_record *tmp;

    while(status_error_flag == 0)
    {
      tmp = new fastq_record;
      tmp->read_next_record(infile, processing_flag);

      if (infile.fail_reason1())
      {
        std::cerr << "ERROR: Found malformed fastq record starting in line "
        << infile.line()-4 << " of the input file "
        << filename  << std::endl;
        tmp->print(std::cerr);
        return -1;
      }
      else if (infile.fail() ) // This includes eof.
        status_error_flag = 1;
      else
        list_of_sequences.push_back(tmp);
    }
    return 0;
  }


  void print(std::ostream &os)
  {
    std::list<fastq_record *>::iterator it, it_end;
    it     = list_of_sequences.begin();
    it_end = list_of_sequences.end();
    //    int count = 0;
    
    while(it != it_end)
    {
      //      os << "Number " << count << std::endl;
      (**it).print(os);
      ++it;
      //      ++count;
    }
  }

  void add_sequences_to_CSequences_object(CSequences3 *pseqs)
  {
    std::list<fastq_record *>::iterator it, it_end;
    it     = list_of_sequences.begin();
    it_end = list_of_sequences.end();
    
    while(it != it_end)
    {
      fastq_record *tmp = *it;
      pseqs->add_seq_to_dataset(SeqType_dna, tmp->get_identifier(), tmp->get_seq(), 'N');
      delete tmp;
      ++it;
      //      ++count;
    }
  }

  unsigned size()
  {
    return list_of_sequences.size();
  }
  
};





#endif
