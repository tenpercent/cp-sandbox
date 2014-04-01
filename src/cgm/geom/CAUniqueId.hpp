//- Class:          CAUniqueId
//- Owner:          Greg Nielson
//- Description:    Cubit attribute for mesh interval.
//- Checked by:
//- Version:

#ifndef CA_UNIQUE_ID_HPP
#define CA_UNIQUE_ID_HPP

#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CADefines.hpp"

class CUBIT_GEOM_EXPORT CAUniqueId: public CubitAttrib
{
private:

  int uniqueId;
    //- the unique id tag for an entity

  static DLIList<CAUniqueId *> allCAUniqueIds;
    //- list of all CAUI's; used in actuate_all function

  static bool autoUniqueId;
    //- flag controlling whether uids are automatically created (even when no other 
    //- CA's request them)

public:

  virtual ~CAUniqueId();

  CAUniqueId(RefEntity*);

  CAUniqueId(RefEntity*, CubitSimpleAttrib*);
    //- make a CAMI from a simple attribute

  //HEADER- RTTI and safe casting functions.
  virtual const type_info& entity_type_info() const
     { return typeid(CAUniqueId);}
  //R- The geometric modeler type
  //- This function returns the type of the geometric modeler.


  CubitStatus actuate();

  CubitStatus update();

  CubitStatus reset() {return CUBIT_SUCCESS;};
    //- don't need to do anything, as all the data gets assigned
    //- and not appended

  CubitSimpleAttrib* cubit_simple_attrib();

  int unique_id() { return uniqueId;}

  void unique_id (int id) {uniqueId = id;};

  int int_attrib_type() {return CA_UNIQUE_ID;};

  static CubitStatus actuate_all();
    //- actuate all the CAUI's on the list, then empty the list

  static bool auto_unique_id();
  static void auto_unique_id(const bool flag);
    //- get/set autoUniqueId

  virtual void print();
    //- print the value of this attribute
};

inline bool CAUniqueId::auto_unique_id()
{
  return autoUniqueId;
}

inline void CAUniqueId::auto_unique_id(const bool flag) 
{
  autoUniqueId = flag;
}

CubitAttrib* CAUniqueId_creator(RefEntity* entity, CubitSimpleAttrib *p_csa);

#endif

