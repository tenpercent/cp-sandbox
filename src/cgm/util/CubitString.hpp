//- Class: CubitString
//- 
//- Description: This file defines the CubitString class which is a 
//- simple implementation of a character string. Basic functionality 
//- is provided as well as equality/inequality tests, subscripting into 
//- a string, searching a string for another string, and stream I/O.
//- This class was written by someone else, but I forget who and cannot
//- find where I got it from. I am still looking.
//-
//- Owner: Greg Sjaardema
//- Author: <unknown at this time, I am looking>
//- Checked by: Jim Hipp, August 31, 1994
//- Version: $Id: 

#if !defined(STRING_HPP)
#define STRING_HPP

#include <assert.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "CubitDefines.h"

#include <iostream>

#include "CubitUtilConfigure.h"

class CubitStringRep;

const size_t MAX_POS = static_cast<size_t>(-1);    // "SIZE_T_MAX"

class CUBIT_UTIL_EXPORT CubitString
{
public:
  //- Heading:  Constructors / Destructor

  CubitString();
  //- Default constructor

  CubitString(const CubitString& s);
  //- Copy Constructor

  CubitString(const char *s);
  //- Create a string from a char*

  CubitString(const int i);
  //- Create a string from a integer

  CubitString(const double f, const unsigned int s_length = 0, const unsigned int sig_digits = 0);
  //- Create a string from a double.  
  //- Use either fixed point or scientific notation, whichever is shorter.
  //- s_length is the maximum string length: If s_length > 0, then
  //- string will contain no spaces and be close to s_length long
  //- without going over. Hence precision is variable.

  ~CubitString();
  //- Default destructor

  CubitString& operator=(const CubitString& s);
  //- Assignment
  
  CubitString& operator+=(const CubitString& s);
  CubitString& operator+=(const char *c);
  CUBIT_UTIL_EXPORT friend CubitString operator+(const CubitString& s1, const CubitString& s2);
  CUBIT_UTIL_EXPORT friend CubitString operator+(const CubitString& s1, const char *c2);
  CUBIT_UTIL_EXPORT friend CubitString operator+(const char *c1, const CubitString& s2);
  //- Concatenation

  CUBIT_UTIL_EXPORT friend bool operator<=(const CubitString&, const CubitString&);
  CUBIT_UTIL_EXPORT friend bool operator>=(const CubitString&, const CubitString&);
  CUBIT_UTIL_EXPORT friend bool operator<(const CubitString&, const CubitString&);
  CUBIT_UTIL_EXPORT friend bool operator>(const CubitString&, const CubitString&);
  
  bool operator==(const CubitString &s) const;
  bool operator!=(const CubitString &s) const;
  bool operator==(const char *s) const;
  bool operator!=(const char *s) const;
  //- Predicates: Compare equality/inequality with other CubitStrings or 
  //- with pointers to a list of characters.
  
  char get_at(size_t pos) const;
  void put_at(size_t pos, char c);
  CubitString substr(size_t first, size_t count = MAX_POS) const;
  CubitString segmentstr(size_t first, size_t first_not = MAX_POS) const;
	//- Added by RY of CAT 11/10/00
	//- substr returns CubitString from first to the next count characters in the string
	//-       i.e "abcdefghijk", pos = 2, n = 6, returns "cdefgh"
	//- segmentstr returns CubitString from first to the character just before first_not.
	//-       i.e "abcdefghijk", pos = 2, n = 6, returns "cdef"
  
  void to_lower();
  static void to_lower(char *string);
  void to_upper();
  static void to_upper(char *string);
  static void trim(char* &string);
  //- Subscripting
  
  size_t find(const CubitString& s, size_t pos = 0) const;
  size_t find_first_of(const CubitString& s, size_t pos = 0) const;
  size_t find_first(char c, size_t pos = 0) const;
  size_t find_last (char c, size_t pos = MAX_POS) const;
#ifdef BOYD15
  size_t find_first_not_of(const CubitString& s, size_t pos = 0) const;
  //- Searching
#endif
  CUBIT_UTIL_EXPORT friend std::ostream & operator<<(std::ostream &, const CubitString&);
  CUBIT_UTIL_EXPORT friend std::istream & operator>>(std::istream & is, CubitString&);
  //- I/O

  size_t length() const;
  const char *c_str() const;
  //- Miscellaneous
  
private:
  CubitStringRep *rep;
};

inline void CubitString::to_lower(char *string) 
{
    // convert this string to lower case
  char *p = string;
  while (*p != '\0')
  {
    *p = tolower (*p);
    p++;
  }
}

inline void CubitString::to_upper(char *string) 
{
    // convert this string to upper case
  char *p = string;
  while (*p != '\0')
  {
    *p = toupper (*p);
    p++;
  }
}

inline void CubitString::trim(char* &string)
{
	while(isspace(*string))
		string++;
	char *p = (string + strlen(string)) - 1;
	while(isspace(*p))
		p--;
	p++;
	*p = '\0';
}

#endif

