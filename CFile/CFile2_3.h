#ifndef CFILE_H
#define CFILE_H

#include <fstream>
#include <string>
#include "../faststring3.h"

// Should not exceed 1 000 000 without making sure that the stack is large enough.
#ifndef BUFFERSIZE
#define BUFFERSIZE 1000000
#endif
#define UNDOs      1


// Changes:
// 22.08.2009: Added getline function for faststring. Disadvantage: additional dependence on faststring.h
// 07.01.2011: Changed Version to CFile2_1.h
// 07.01.2011: Now uses the faststring2.h
// 26.12.2020: Corrected fail_flag usage:
//             (1) Added rdstate() function to make it compatible with ifstream objects.
//             (2) Reading file to string. Do not fail if >0 chars read.
//             (3) getline functions: Do fail if thing read, except delim == first char. 

// Idea: ios::   set Buffer size
//               Maybe use fstream fail flags instead of own flags. Not sure this is a good idea. Consequences are not thought out.

// TODO: Old mac format not supported.
//       This requires to allow two successive calls to the internal ungetchar command.


class CFile : private std::ifstream
{
 private:
  // There is a problem here: The stack size is limited. One should consider to move this to the heap.
  // allocate it with malloc, new or declare it as static - the latter is not good in a class since a
  // static class member is the same for all instances of the class, which we do not want here.

  char      buffer[BUFFERSIZE];
  char*     buffer_end;
  char*     buffer_pos;
  
  char      __status;
  unsigned  __line;

  //  bool      __eof;
  //  bool      __openFailed;

  unsigned fill_buffer(unsigned overlap)
  {
    unsigned          good_overlap = buffer_end - buffer;
    std::streamsize   n;

    if (good_overlap < overlap)
    {
      overlap = good_overlap;
    }

    if (overlap > 0)
      std::memmove(buffer, buffer_end - overlap, overlap);

    std::ifstream::read(buffer + overlap, BUFFERSIZE - overlap);
    n = std::ifstream::gcount();

    if ( n == 0 )
    {
      __status |=  __eof_flag;     // Set eof flag
      __status &= ~__good_flag;    // Unset good flag
      __status |=  __fail_flag;    // Setting the fail flag is not always correct. Needs to be unset in routines that read more than one char. 
    }

    buffer_pos = buffer+overlap;
    buffer_end = buffer_pos + n;

    return overlap;
  }

  char getchar_intern()  // should only be done if we did not fail yet!!!! This would allow us to recover!!!
  {
    if ( buffer_pos == buffer_end )
    {
      fill_buffer(UNDOs);
      if (__status & __fail_flag)
	return '\0';
    }
    return *buffer_pos++;
  }

  void ungetchar_intern()  // should only be done if we did not fail yet!!!! This would allow us to recover!!!
  {
    if (buffer_pos != buffer)
      --buffer_pos;
    else
    {
      __status |=  __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
    }
  }

 public:

  enum {__eof_flag = 1, __good_flag = 2,  __fail_flag = 4, __bad_flag = 8, __fail_reason1 = 16, __fail_reason2 = 32, __fail_reason3 = 64,
	__fail_reason4 = 128};

  void open(const char *name)
  {
    ffopen(name);
  }

  void open(std::string name)
  {
    ffopen(name);
  }

  void ffopen(const char *name)
  {
    __status = __good_flag;

    std::ifstream::open(name);

    if ( std::ifstream::fail() )
    {
      __status |=  __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
    }

    __line        = 1;
    buffer_end    = buffer;
    buffer_pos    = buffer;
  }

  void ffopen(std::string name)
  {
    ffopen(name.c_str());
  }

  void ffclose()
  {
    std::ifstream::close();
  }

  void close()
  {
    ffclose();
  }

  bool exists()  // -- deprecated -- do not use this any more - check fail() instead !!!!!!!!
  {
    return !fail(); 
  }

  bool fail()
  {
    return (__status & __fail_flag);
  }

  bool good()
  {
    return (__status & __good_flag);
  }

  unsigned line()
  {
    return __line;
  }

  bool eof()
  {
    return (__status & __eof_flag);
  }

  // Deprecated
  char status()
  {
    return __status;
  }

  char rdstate()
  {
    return __status;
  }
  
  bool fail_reason1()
  {
    return (__status & __fail_reason1);
  }

  bool fail_reason2()
  {
    return (__status & __fail_reason2);
  }

  bool fail_reason3()
  {
    return (__status & __fail_reason3);
  }

  bool fail_reason4()
  {
    return (__status & __fail_reason4);
  }

  void rewind()
  {
    clear();

    __line        = 1;
    buffer_end    = buffer;
    buffer_pos    = buffer;

    std::ifstream::seekg (0, std::ios::beg);
    //    clear();
  }

  void clear(char s = __good_flag)
  {
    __status = s;
    std::ifstream::clear();
  }

  void ungetchar()
  {
    if (*(buffer_pos-1) == '\n' && !(__status & __fail_flag))
      --__line;
    ungetchar_intern();
  }

  void ignore(int delim = -1)
  {
    while (__status == __good_flag && getchar() != delim ){}
  }

  char peekchar()
  {
    char c = getchar();
    ungetchar();
    return c;
  }

  char peek()
  {
    return peekchar();
  }

  // Will set the fail flag if no char can be returned.
  // Functions such as those in CSequence_Mol2_1.h should unset the fail flag if reading the sequence was successful in the end.
  char getchar()
  {
    char c;

    c = getchar_intern();

    if ( c < 14 && !(__status & __fail_flag) )
    {
      if (  c == '\r' )
      {
	// Overwrite the last reading position that contained the \r with \n
	*(buffer_pos-1) = '\n';           // Should always be valid! Works for 1 Undo
	c = getchar_intern();
	if ( c != '\n' && !(__status & __fail_flag) )
	{
	  // 	std::cerr << "Old mac file format currently not supported." << '\n'
	  ungetchar();    /* old mac format, else dos format     */
	}
	c = '\n';
	++__line;
      }
      else if ( c == '\n')
	++__line;
    }
    return c;  
  }


  char getrawchar()
  {
    return getchar_intern();
  }

  void getline(faststring& str, char delim='\n')
  {
    char c = getchar();
    str.clear();
    while ( c != delim && !(__status & __fail_flag) )
    {
      str.push_back(c);
      c = getchar();
    }
    if ((__status & __fail_flag) && (str.size() > 0 || c == delim) )
      __status &= ~__fail_flag; // Unset fail flag by using & on the complement of the fail flag;
  }

  void getline(std::string& str, char delim='\n')
  {
    char c = getchar();
    str="";
    while ( c != delim && !(__status & __fail_flag) )
    {
      str.push_back(c);
      c = getchar();
    }
    if ((__status & __fail_flag) && (str.size() > 0 || c == delim) )
      __status &= ~__fail_flag; // Unset fail flag by using & on the complement of the fail flag;
  }

  void getline(char* cstr, unsigned n, char delim='\n')
  {
    char     c = getchar();
    unsigned i = 0;

    while ( !(__status & __fail_flag ) && i < n-1 &&  c != delim)
    {
      cstr[i] = c;
      ++i;
      c = getchar(); 
    }

    if ((__status & __fail_flag) && (i > 0 || c == delim) )
      __status &= ~__fail_flag; // Unset fail flag by using & on the complement of the fail flag;
    cstr[i] = '\0';
  }

  void readFileIntoString(faststring &str)
  {
    char c;
    
    str.clear();
    
    c = getchar();
    while (!(__status & __fail_flag))
    {
      str.push_back(c);
      c = getchar();
    }
    if (str.length() > 0)
      __status &= ~__fail_flag;  // Unset fail flag by using & on the complement of the fail flag;
  }
  
  // Double check that the fail_flag is set correctly.
  char lastchar()
  {
    char c;

    if ( buffer_pos != buffer_end )
      c = *(buffer_pos-1);
    else
    {
      __status |= __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
      return 0;
    }

    /* It could be a '\r' in mac format */
    if (c == '\r')
      c = '\n';

    return c;
  }

  char lastrawchar()
  {
    char c;

    if ( buffer_pos != buffer_end )
      c = *(buffer_pos-1);
    else
    {
      __status |= __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
      return 0;
    }
    return c;
  }


};




#endif
