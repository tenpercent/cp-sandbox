//-------------------------------------------------------------------------
// Filename      : RandomMersenne.hpp 
//
// Purpose       : Main include for random number generator
//                 This random number generator came from the authors of the Mersenne
//                 random number generator, Makoto Matsumoto and Takuji Nishimura.
//                 See the definition file RandomMersenne.cpp for the full copyright
//                 statement.
//
//                 I have made the following changes, starting from the (faster version
//                 of the) source code provided by the authors.  I have made the following
//                 changes:
//                 - put the generator into a C++ class, to simplify use of inline functions
//                   and member instead of global variables
//                 - made the period parameters enum variables instead of #defines,
//                   formerly static state variables member variables, and formerly
//                   global functions member functions.
//                 - removed check for initialization from next_state, since I require
//                   initialization when the class object is created
//
//                 Other Notes:
//                 - I made the no-argument constructor private, to force input of
//                   a seed value.  The random sequence gets initialized upon construction.
//                 - function names in this class don't follow our naming convention, because
//                   I wanted to stay as close to the original source code as possible
//                 - most of the public functions in this class are inlined, in the hopes that
//                   at least some compilers will honor that; these functions all call a
//                   non-inlined function next_state() though.
//
//
// Creator       : Tim Tautges (mostly code from aforementioned site)
//
// Date          : 08/02
//
// Owner         : Tim Tautges
//-------------------------------------------------------------------------

#ifndef RANDOM_MERSENNE_HPP
#define RANDOM_MERSENNE_HPP

#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT RandomMersenne {                  // encapsulate random number generator

public:

    //! constructor initialized random number generator, so it takes a seed
  RandomMersenne(unsigned long seed);

    //! constructor initialized random number generator, so it takes a seed array
  RandomMersenne(unsigned long init_key[], unsigned long key_length);

    //! initialize by an array with array-length
    //! init_key is the array for initializing keys
    //! key_length is its length
  void init_by_array(unsigned long init_key[], unsigned long key_length);

    //! generates a random number on [0,0xffffffff]-interval 
  unsigned long genrand_int32();
  
    //! generates a random number on [0,0x7fffffff]-interval 
  long genrand_int31();

    //! generates a random number on [0,1]-real-interval 
  double genrand_real1();

    //! generates a random number on [0,1)-real-interval 
  double genrand_real2();

    //! generates a random number on (0,1)-real-interval 
  double genrand_real3();

    //! generates a random number on [0,1) with 53-bit resolution
  double genrand_res53();

private:

    //! private constructor to prevent constructing without seed
  RandomMersenne() {}

    //! initializes state[N] with a seed
  void init_genrand(unsigned long s);

  enum {
    N = 624, M = 397,
    MATRIX_A = 0x9908b0dfU, //! constant vector a    
    UPPER_MASK = 0x80000000U, /* most significant w-r bits */
    LOWER_MASK = 0x7fffffffU /* least significant r bits */
  };

  unsigned long mt[N]; //! the array for the state vector  
  unsigned mti;
  int left;
  int initf;
  unsigned long *next;

};    

inline RandomMersenne::RandomMersenne(unsigned long seed) 
{
  init_genrand(seed);
}

inline RandomMersenne::RandomMersenne(unsigned long init_key[], 
                                      unsigned long key_length)
{
  init_by_array(init_key, key_length);
}

/* generates a random number on [0,0x7fffffff]-interval */
inline long RandomMersenne::genrand_int31()
{
  return (long) (genrand_int32() >> 1);
}

#endif  

