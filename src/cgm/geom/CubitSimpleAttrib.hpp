//- Class: CubitSimpleAttrib
//- Owner: Greg Nielson
//- Description: This class is used to carry information from
//-              the cubit attribute classes to the solid model
//-               attributes.
//- Checked By:
//- Version: $Id:

#ifndef CUBIT_SIMPLE_ATTRIB_HPP
#define CUBIT_SIMPLE_ATTRIB_HPP

#if defined(WIN32)
#pragma warning(push)
#pragma warning(disable : 4251)  // hide warnings about template dll exports
#endif

#include "CubitString.hpp"
#include "CubitDefines.h"
#include "DLIList.hpp"
#include "CubitGeomConfigure.h"

class CUBIT_GEOM_EXPORT CubitSimpleAttrib {
private:
  CubitString characterType;  // defines the attribute type 

// This is the old data before the lists below were added.
// This data will be removed and will be added into the lists
// rather than keeping it separate.  The characterType above will
// allow to distinguish between "NEW_SIMPLE_ATTRIB" and older types
//   Byron 1/16/03
//
//  CubitString stringData;  // attribute data of type string
//  CubitString moreStringData;  //  attribute data of type string
//  double doubleData;  // attribute data of type double
//  int integerData;  // attribute data of type int

// THE FOLLOWING MEMBERS ARE USED IF 
// characterType IS EQUIVALENT TO "NEW_SIMPLE_ATTRIB".
// THE ATTRIBUTE TYPE IS DEFINED BY THE FIRST
// MEMBER OF THE stringDataList LIST.

  DLIList<CubitString*> stringDataList;  // List of CubitString attribute data

  DLIList<double*> doubleDataList;  // List of double attribute data

  DLIList<int*> intDataList;  // List of int attribute data

  static CubitBoolean pushAttribs;
    //- if true, when saving attributes, push them to non-primary bridges
    //- on merged ref entities
  
  void set( const DLIList<int>& int_data );
  void set( const DLIList<double>& real_data );
  void set( const DLIList<CubitString*>& string_data );
    //- Keep these private.  When re-using an existing CSA,
    //- it is too likely that the caller will forget to call
    //- all three methods, some possibly with empty lists, 
    //- ensuring all old data is cleared out.  Instead, force
    //- caller to use the public 'set' method and pass all
    //- data at once.  Also, these methods treat the passed
    //- lists as circular: copying them beginning at the current
    //- index NOT the reset position.
  
public:
    void character_type(CubitString new_type);
  //- sets (or resets) the attribute type

  CubitSimpleAttrib();
  //- default constructor

  CubitSimpleAttrib(CubitString new_type);
  //- constructor use in save/restore
  
  CubitSimpleAttrib(const CubitString new_character_type,
                    const CubitString new_string_data,
                    const CubitString new_more_string_data,
                    const int new_integer_data = 0,
                    const double new_double_data = 0.);
  //- constructor

  CubitSimpleAttrib(const char* new_character_type,
                    const char* new_string_data = NULL,
                    const char* new_more_string_data = NULL,
                    const int *new_integer_data = NULL,
                    const double *new_double_data = NULL);
  //- constructor  

  CubitSimpleAttrib(DLIList<CubitString*> *string_list,
                    DLIList<double> *double_list = NULL,
                    DLIList<int> *int_list = NULL);
  //- constructor  

  CubitSimpleAttrib(CubitSimpleAttrib *csa_ptr);
    //- copy constructor, sort of
  
  ~CubitSimpleAttrib();
   //- destructor

  CubitSimpleAttrib& operator=(CubitSimpleAttrib &other);


  void initialize_from_lists_of_ptrs(DLIList<CubitString*> *string_list,
                                     DLIList<double*> *double_list,
                                     DLIList<int*> *int_list);

  void set( const CubitString&            type,
            const DLIList<int>*           int_data    = 0,
            const DLIList<double>*        real_data   = 0,
            const DLIList<CubitString*>*  string_data = 0 );
    //- Unlike the 'initialize_from_lists' method below,
    //- this method will clear existing data from the
    //- lists before adding the new (and reuse the 
    //- allocated heap space.)  
    //-
    //- NOTE:  If one of passed list pointers is NULL, all
    //-        data of the corresponding type will be 
    //-        removed from this attribute.
    //- NOTE:  Lists are treated as circular.  The copying
    //-        of data into this attribute begins at the
    //-        CURRENT position, NOT THE RESET position of
    //-        the list.
  
  void initialize_from_lists(DLIList<CubitString*>* string_list,
                             DLIList<double*>* double_list,
                             DLIList<int*>* int_list);
    //- used to copy one csa to another
    //- This really does "append_from_lists". (i.e. It appends
    //- to the data already in this attribute, if any.)
 
  void initialize_from_lists(DLIList<CubitString*> *string_list,
                             DLIList<double> *double_list,
                             DLIList<int> *int_list);
    //- used to copy one csa to another

  CubitString character_type(); // {return characterType;}
  //- returns the type of attribure

  
  static void initialize_settings();

  DLIList<CubitString*>* string_data_list() {return &stringDataList;} 

  void string_data_list( DLIList<CubitString*>* new_list); 

  DLIList<double*>* double_data_list() {return &doubleDataList;}

  void double_data_list(const DLIList<double*>* new_list);

  DLIList<int*>* int_data_list() {return &intDataList;}

  void int_data_list(const DLIList<int*>* new_list);

  
  static CubitBoolean equivalent(CubitSimpleAttrib *first_attrib_ptr,
                                 CubitSimpleAttrib *second_attrib_ptr);
    //- return true if the two csa's are equivalent

  void print();
    //- print the contents of this simpleattrib

#ifdef BOYD14
  CubitStatus save(FILE *save_file);
  CubitStatus restore(FILE *restore_file);
#endif

  static CubitBoolean get_push_attribs() {return pushAttribs;};
  static void set_push_attribs(CubitBoolean flag) {pushAttribs = flag;};
    //- get/set pushAttribs
  
};


#if defined(WIN32)
#pragma warning(pop)
#endif

#endif








