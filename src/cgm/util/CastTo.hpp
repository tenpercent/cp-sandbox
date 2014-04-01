//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename      : CastTo.hpp
//
// Purpose       : Provides a macro to convert an object from one type to
//                 any other type, provided the cast is deemed safe.
//
// Special Notes : Uses the get_address(NewType) function to determine whether
//                 a cast is safe or not and then to return the appropriate
//                 address. 
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/09/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef CAST_TO_HPP
#define CAST_TO_HPP

// Macro to cast a pointer to one class to a pointer to another class.
// Constraint: The first class must support a function get_address() and the
// second class must have an enum entry corresponding to it who's syntax
// is "classTwo"_TYPE.
//
// NOTE: Multiple inheritance and virtual inheritance are accounted for
//       appropriately in the get_address functions.
//
// NOTE: If a NULL pointer is given to CAST_TO, it always returns NULL.
//
// NOTE: If the get_address function returns a NULL, then a NULL is returned.

#define CAST_TO(classOnePtr, classTwo) \
   dynamic_cast<classTwo*>(classOnePtr)

#define STATIC_CAST_TO(classOnePtr, classTwo) \
   static_cast<classTwo*>(classOnePtr)

#define CAST_LIST_TO_PARENT(classOneList, classTwoList) \
   { \
      (classOneList).reset(); \
      (classTwoList).clean_out(); \
      for (int CAST_LIST_IT = (classOneList).size(); CAST_LIST_IT > 0; CAST_LIST_IT--) \
        (classTwoList).append((classOneList).get_and_step()); \
   }
 
#define CAST_LIST(classOneList, classTwoList, classTwo) \
   { \
      classTwo *CAST_LIST_PTR_TWO; \
      (classOneList).reset(); \
      (classTwoList).clean_out(); \
      for (int CAST_LIST_IT = (classOneList).size(); CAST_LIST_IT > 0; CAST_LIST_IT--) { \
        CAST_LIST_PTR_TWO = dynamic_cast<classTwo*>( (classOneList).get_and_step() ); \
        if (CAST_LIST_PTR_TWO != NULL) (classTwoList).append(CAST_LIST_PTR_TWO); \
      } \
   }

#define CONST_CAST_LIST(classOneList, classTwoList, classTwo) \
   { \
      classTwo *CAST_LIST_PTR_TWO; \
      (classTwoList).clean_out(); \
      for (int CAST_LIST_IT = 0; CAST_LIST_IT < (classOneList).size(); CAST_LIST_IT++) { \
        CAST_LIST_PTR_TWO = dynamic_cast<classTwo*>( (classOneList).next(CAST_LIST_IT) ); \
        if (CAST_LIST_PTR_TWO != NULL) (classTwoList).append(CAST_LIST_PTR_TWO); \
      } \
   }
#endif

