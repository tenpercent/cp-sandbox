//- Class:          CAMergePartner
//- Owner:          Greg Nielson
//- Description:    Cubit attribute for merge partners.
//- Checked by:
//- Version:

#ifndef CA_MERGE_PARTNER_HPP
#define CA_MERGE_PARTNER_HPP

#include "CubitAttrib.hpp"
#include "CubitString.hpp"
#include "CubitDefines.h"
#include "CADefines.hpp"

class RefEntity;
class RefEntity;
template <class X> class DLIList;

class CUBIT_GEOM_EXPORT CAMergePartner: public CubitAttrib
{

private:
  int mergeID;
    //- a unique id to verify proper merging

  int isSurvivor;
  CubitSense bridge_sense_;

public:
  CAMergePartner(RefEntity* );

  CAMergePartner(RefEntity*, CubitSimpleAttrib*);
    //- construct a CAMP from a simple attribute

  void initialize();
    //- initialize the data for this CA

  virtual ~CAMergePartner();

  //HEADER- RTTI and safe casting functions.
  virtual const type_info& entity_type_info() const
     { return typeid(CAMergePartner);}
  //R- The geometric modeler type
  //- This function returns the type of the geometric modeler.

  CubitStatus actuate();
  
  CubitStatus actuate_list(DLIList<RefEntity*>);

  CubitStatus update();

  CubitStatus reset();
    //- just reset the merge id

  CubitSimpleAttrib* cubit_simple_attrib();
  
  CubitSimpleAttrib* cubit_simple_attrib(CubitString);
  
  int merge_id(){return mergeID;}

  int int_attrib_type() {return CA_MERGE_PARTNER;}
    //- returns the enumerated attribute type
    
  CubitSense bridge_sense()
    { return bridge_sense_; }
  
  static void set_survivor( CubitSimpleAttrib* csa, int is_survivor );
  static CubitBoolean is_survivor( CubitSimpleAttrib* csa );

  static void set_bridge_sense(CubitSimpleAttrib* csa, CubitSense sense);
  static CubitSense get_bridge_sense( CubitSimpleAttrib* csa_ptr );
  
  static void set_saved_id( CubitSimpleAttrib* csa, int id );
  static int get_saved_id( CubitSimpleAttrib* csa );
    // returns zero if not set.
  
  void merge_prepare(DLIList<RefEntity*> &merge_list);
    //- pass back a list of mergable entities with same unique id

  virtual void print();
    //- print the info in this attribute

};

inline CubitStatus CAMergePartner::reset() 
{
  mergeID = -1;
  hasUpdated = CUBIT_FALSE;
  hasActuated = CUBIT_FALSE;
  hasWritten = CUBIT_FALSE;
  deleteAttrib = CUBIT_FALSE;
  bridge_sense_ = CUBIT_UNKNOWN;
  return CUBIT_SUCCESS;
}

CubitAttrib* CAMergePartner_creator(RefEntity* entity, CubitSimpleAttrib *p_csa);

#endif

