//- Class:          CAEntityName
//- Owner:          Greg Nielson
//- Description:    Cubit attribute for entity names.
//- Checked by:
//- Version:

#ifndef CA_ENTITY_NAME_HPP
#define CA_ENTITY_NAME_HPP

#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CubitDefines.h"
#include "CADefines.hpp"

class RefEntity;

class CUBIT_GEOM_EXPORT CAEntityName: public CubitAttrib
{
private:
 
  DLIList<CubitString*> entityNames;

#ifdef BOYD17
  DLIList<CubitString*> entityOptions;
#endif

public:

  virtual ~CAEntityName();

  CAEntityName(RefEntity*);

  CAEntityName(RefEntity*, CubitSimpleAttrib *);
    //- create a CAEN from a simple attribute

  //HEADER- RTTI and safe casting functions.
  virtual const type_info& entity_type_info() const
     { return typeid(CAEntityName);}
  //R- The geometric modeler type
  //- This function returns the type of the geometric modeler.

  CubitStatus actuate();

  CubitStatus update();

  CubitStatus reset();
    //- reset function, cleans out name lists

  CubitSimpleAttrib *split_owner();

  void merge_owner(CubitAttrib *deletable_attrib);

  CubitSimpleAttrib* cubit_simple_attrib();

  int int_attrib_type() {return CA_ENTITY_NAME;}

#ifdef BOYD14
  void remove_caentity_name(CubitString old_entity_name);

  void remove_caentity_option(CubitString old_entity_option);

  void add_csa_data(CubitSimpleAttrib *csa_ptr);
    //- add data in csa_ptr to this attrib; for compatibility with old
    //- sat files, which have multiple name attributes
#endif

  void print();
  
};

CubitAttrib* CAEntityName_creator(RefEntity* entity, CubitSimpleAttrib *p_csa);

#endif

