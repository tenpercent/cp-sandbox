//- Class:       CubitAttribUser
//- Owner:       Greg Nielson
//- Description: Base class inherited by the RefEntity class that
//-              gives attribute functionality to the RefEntitys.
//- Checked By:
//- Version: 

#ifndef CUBIT_ATTRIB_USER_HPP
#define CUBIT_ATTRIB_USER_HPP

// #include "CubitString.hpp"

class RefEntity;

#include <typeinfo>
#if !defined(NT)
using std::type_info;
#endif

#include "CubitDefines.h"
#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CubitSimpleAttrib.hpp"
#include "CubitGeomConfigure.h"

class TopologyBridge;

//template <class X> class DLIList;
//class CubitAttrib;
//class CubitSimpleAttrib;
//class CAEntityName;


class CUBIT_GEOM_EXPORT CubitAttribUser
{
private:
//  CubitStatus remove_attrib_geometry_entity (CubitSimpleAttrib*);
    CubitAttribUser( const CubitAttribUser& );
    void operator=( const CubitAttribUser&);

  
protected:
  CubitAttrib* headAttrib;

public:
// bad function CubitAttribUser() {headAttrib = NULL;}

  virtual ~CubitAttribUser();

  CubitAttribUser(CubitAttrib* = NULL);

  CubitAttrib* get_cubit_attrib (int attrib_type,
                                 CubitBoolean create_if_missing = CUBIT_TRUE);

  CubitStatus add_cubit_attrib (CubitAttrib* cubit_attrib_ptr);
    //- appends this attrib to the list

private:
//________  Change Code by DZ of Cat,  3/11/99 12:19:10 PM  ________
  CubitStatus put_simple_attrib (CubitSimpleAttrib* cubit_simple_attrib_ptr,
     CubitBoolean append_it = CUBIT_TRUE);
//________  Change End by DZ of Cat,  3/11/99 12:19:10 PM  ________

public:
  void append_simple_attribute(TopologyBridge* bridge, CubitSimpleAttrib* attrib_ptr);
  void append_attrib_internal( TopologyBridge* bridge, CubitSimpleAttrib* attrib_ptr );
  void remove_all_simple_attribute(TopologyBridge* bridge);
  void remove_simple_attribute(TopologyBridge* bridge, CubitSimpleAttrib* attrib_ptr);

  CubitStatus clear_simple_attribs();
    // Remove all CubitSimpleAttrib from TopologyBridges.
  
public:
  CubitStatus auto_read_cubit_attrib ();
    //moves all attributes from the SME to CUBIT

  CubitStatus read_cubit_attrib (int attrib_type);
    //moves a specific cubit attribute type from the SME to CUBIT
  
  //CubitStatus read_cubit_attrib(CubitBoolean read_children = CUBIT_FALSE);
    // reads all attribs from SME to CUBIT, and optionally for children
  
  CubitStatus write_specific_cubit_attrib (CubitAttrib* cubit_attrib_ptr);
    //moves a specific cubit attribute from CUBIT to the SME
  
  CubitStatus write_cubit_attrib_by_type (int attrib_type);
    //moves a specific cubit attribute type from CUBIT to the SME

  CubitStatus write_cubit_attrib_list (DLIList<CubitAttrib*> attrib_list);
    //moves a list of cubit attributes from CUBIT to the SME

  CubitStatus write_cubit_attribs ();
    //moves all cubit attributes from CUBIT to the SME

public:
  void split_owner(DLIList<CubitSimpleAttrib*> &csa_list);
    //- if owner is to be split, get simple attribs for new entity

  void merge_owner(RefEntity *deletable_entity);
    //- if owner is to be merged, combine attribs from deletable_entity with this

  void transf_owner(const CubitVector &matrow1,
                    const CubitVector &matrow2,
                    const CubitVector &matrow3,
                    const CubitVector &translate_vec,
                    const double scale_factor);
    //- called if owner is to be transformed

private:
  void auto_create_for_merge(RefEntity *deletable_entity);
    //- create any attribs on deletable_entity and not on this

public:
  //static CubitBoolean cubit_simple_attrib_equivalent(CubitSimpleAttrib*,
  //                                                   CubitSimpleAttrib*);

  //int cubit_attrib_equivalent(CubitAttrib*,CubitAttrib*);

  //CubitStatus make_simple_attrib (CubitSimpleAttrib**,
  //                                CubitAttrib*);

#ifdef BOYD14
  //int cubit_attrib_list_length ();
#endif
  

  CubitStatus actuate_cubit_attrib (int attrib_type);

private:
  CubitStatus actuate_cubit_attrib (DLIList<CubitAttrib*>);
  CubitStatus actuate_cubit_attrib (CubitAttrib*);

public:
  static CubitStatus actuate_cubit_attrib (DLIList<RefEntity*>,int);
    //- actuate attributes of the specified type on entities in the list

  //CubitStatus actuate_cubit_attrib ();

  CubitStatus auto_actuate_cubit_attrib (CubitBoolean from_constructor = CUBIT_TRUE,
                                         CubitBoolean after_geom_changes = CUBIT_FALSE);

  //CubitStatus update_cubit_attrib ();

  //CubitStatus update_cubit_attrib (CubitAttrib*);

  CubitStatus update_cubit_attrib (int);

  //CubitStatus update_cubit_attrib (DLIList<CubitAttrib*>);

  CubitStatus auto_update_cubit_attrib ();
    //- for this cau, automatically create and update ca's

public:
  static CubitStatus auto_update_cubit_attrib (DLIList<RefEntity*> &entity_list,
                                               CubitBoolean write_too = CUBIT_TRUE);
    //- for entity_list, auto create, update and write ca's

  static CubitStatus clear_all_simple_attrib( DLIList<RefEntity*>& entity_list );
    //- remove all CubitSimpleAttrib from TopologyBridges.

private:
  static void auto_reset_cubit_attrib(DLIList<RefEntity*> ref_ents);
    //- set the update flag on all attribs on these entities and their children to false

public:
  //void auto_reset_cubit_attrib();
    //- set the update flag on all my attribs to false

  void find_cubit_attrib_type (int, DLIList<CubitAttrib*>&) const;
  
#ifdef BOYD14
  //CubitBoolean cubit_attrib_exists (int type);
#endif

private:
  CubitStatus remove_cubit_attrib (DLIList<CubitAttrib*>);

public:
  CubitStatus remove_cubit_attrib (int attrib_type);

  CubitStatus remove_cubit_attrib (CubitAttrib*);

  //CubitStatus remove_cubit_attrib ();

  CubitStatus remove_attrib_geometry_entity (CubitAttrib*);

  //CubitStatus remove_attrib_geometry_entity ();

  void get_cubit_attrib_list (DLIList<CubitAttrib*>&);

  int num_cubit_attrib();
    // returns number of cubit attributes

  void set_updated_flag(CubitBoolean flag);
#ifdef BOYD14
  //void set_actuated_flag(CubitBoolean flag);
#endif
  void set_written_flag(CubitBoolean flag);
#ifdef BOYD14
  //void set_delete_flag(CubitBoolean flag);
    //- sets the given flags for all CAs on this CAU
#endif

  //void print_attribs();
    //- finds and prints all CA's on this entity
};

#endif
                     
  

  
  
