//- Class: DIntArray
//- Description: DIntArray is an array of ints making use of the vector
//-              container from the STL.  The only modifications 
//-              that have been made are the addition of an append,
//-              size, and overloaded operator[] functions to match
//-              the current use of the class.  
//- Owner: Kirk Walton
//- Checked by:

#ifndef DINTARRAY_HPP
#define DINTARRAY_HPP

#include <vector>
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT DIntArray 
{
public:
  DIntArray(){};
  
  ~DIntArray(){};
  
  void append(int item);
  int size();
  int operator[](int index);
  
private:

#ifndef HP
  std::vector<int> intArray;
#else
  vector<int> intArray;
#endif
  
};

#endif

