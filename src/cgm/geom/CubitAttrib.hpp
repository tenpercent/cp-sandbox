//- Class:          CubitAttrib
//- Owner:          Greg Nielson
//- Description:    pure virtual base class from which all the specific
//-                 cubit attributes are derived.
//- Checked By:
//- Version:

#ifndef CUBIT_ATTRIB_HPP
#define CUBIT_ATTRIB_HPP

#include <typeinfo>
#if !defined(NT)
using std::type_info;
#endif

#include "CubitString.hpp"
#include "CubitDefines.h"
#include "CubitVector.hpp"

#include "CubitAttribManager.hpp"
#include "CGMApp.hpp"
#include "CubitGeomConfigure.h"

class RefEntity;
template <class X> class DLIList;
class CubitSimpleAttrib;
class CubitAttribUser;

// this enum used to index static arrays; first attribute must start at 0,
// last attribute must be CA_LAST_CA (some static arrays dimensioned using
// this item)
/*
enum CubitAttributeType {
    CA_UNDEFINED = -1,
    CA_MERGE_PARTNER = 0, // must be 0!!!
    CA_ENTITY_NAME,
    CA_GROUP,
    CA_ENTITY_ID,
//    CA_PARTITION_VG,
//    CA_COMPOSITE_VG,
//    CA_VIRTUAL_VG,
    CA_UNIQUE_ID,
    CA_DEFERRED_ATTRIB,
    CA_MESH_INTERVAL,
    CA_MESH_RELATIVE_LENGTH,
    CA_GENESIS_ENTITY,
    CA_MESH_SCHEME,
    CA_SMOOTH_SCHEME,
    CA_VERTEX_TYPE,
    CA_MESH_CONTAINER,
    CA_BODIES,
    CA_ENTITY_COLOR,
#ifdef CAT
    CA_VERTEX_FORCE,
    CA_SURFACE_FORCE,
    CA_CURVE_FORCE,
    CA_VERTEX_DISPLACEMENT,
    CA_SURFACE_DISPLACEMENT,
    CA_CURVE_DISPLACEMENT,
    CA_VOLUME_DISPLACEMENT,
    CA_SURFACE_PRESSURE,
    CA_CURVE_PRESSURE,
    CA_SURFACE_TEMPERATURE,
    CA_CURVE_TEMPERATURE,
    CA_VERTEX_TEMPERATURE,
    CA_SURFACE_HEATFLUX,
    CA_CURVE_HEATFLUX,
    CA_SURFACE_CONVECTION,
    CA_CURVE_CONVECTION,
    CA_SURFACE_CONTACT,
    CA_CURVE_CONTACT,
#endif
    CA_MERGE_STATUS,
    CA_LAST_CA // must be last item!!!
};
*/
class CubitAttrib;

class CUBIT_GEOM_EXPORT CubitAttribFactory 
{
public:
  friend class CubitAttrib;
  
protected:
  CubitAttribFactory() {};
  
  virtual ~CubitAttribFactory() {};

private:
  virtual CubitAttrib *create_cubit_attrib(const int attrib_type,
                                           RefEntity *entity) = 0;
  
  virtual CubitAttrib *create_cubit_attrib(CubitSimpleAttrib *csa_ptr,
                                           RefEntity *entity) = 0;
};

class CUBIT_GEOM_EXPORT CubitAttrib
{
public:

  CubitAttrib(RefEntity *attrib_owner);
    //- bare constructor

  virtual ~CubitAttrib();
    //- destructor

//  static CubitAttrib *create_cubit_attrib(const int attrib_type,
//                                          RefEntity *entity);
    //- create a cubit attribute of the specified type and pass it back

//  static CubitAttrib *create_cubit_attrib(CubitSimpleAttrib *csa_ptr,
//                                          RefEntity *entity);
    //- create a cubit attribute from a CSA and pass it back

    //HEADER-: get/set functions

  RefEntity* attrib_owner() {return attribOwnerEntity;}
  void attrib_owner(RefEntity* new_attrib_owner)
    {attribOwnerEntity = new_attrib_owner;}
    //- get/set the attrib owner

  CubitBoolean has_actuated() {return hasActuated;}
  void has_actuated(CubitBoolean set_has_actuated)
    { hasActuated = set_has_actuated;}
    //- get/set hasActuated

  CubitBoolean has_updated() {return hasUpdated;}
  void has_updated(CubitBoolean set_has_updated)
    { hasUpdated = set_has_updated;}
    //- get/set hasUpdated

  virtual CubitBoolean has_written() const;
  virtual void has_written(CubitBoolean set_has_written);
    //- get/set hasWritten
  
  void delete_attrib(CubitBoolean set_remove) {deleteAttrib = set_remove;}
  CubitBoolean delete_attrib() {return deleteAttrib;}
    //- get/set the deleteAttrib flag

  CubitAttrib* next_attrib() {return nextAttrib;}
  CubitStatus set_next_attrib(CubitAttrib* next_attrib_ptr);
    //- get/set the next attribute in the chain

//  static CubitStatus auto_create_attribs(CubitAttribUser *cau);
    //- create attribs whose auto update flag is set
  
//  static CubitStatus auto_update_attribs(CubitAttribUser *cau);
    //- create attribs whose auto update flag is set
  
//  static void auto_actuate_flag(int attrib_type, CubitBoolean value);
//  static CubitBoolean auto_actuate_flag(int attrib_type);
//  static void auto_update_flag(int attrib_type, CubitBoolean value);
//  static CubitBoolean auto_update_flag(int attrib_type);
//  static void auto_write_flag(int attrib_type, CubitBoolean value);
//  static CubitBoolean auto_write_flag(int attrib_type);
//  static void auto_read_flag(int attrib_type, CubitBoolean value);
//  static CubitBoolean auto_read_flag(int attrib_type);
    //- get/set actuate, update, write, and read flags, static versions

//  static CubitBoolean* auto_read_flag_ptr(int attrib_type);
//  static CubitBoolean* auto_actuate_flag_ptr(int attrib_type);
//  static CubitBoolean* auto_update_flag_ptr(int attrib_type);
//  static CubitBoolean* auto_write_flag_ptr(int attrib_type);
      //- get array ptr for actuate, update, write, and read flags, static versions

//  static void auto_flag(int flag);
//  static int auto_flag();
    //- get/set all attribute flags
  //- returns -1, 0 or 1.

  void auto_actuate_flag(CubitBoolean value);
  CubitBoolean auto_actuate_flag();
  void auto_update_flag(CubitBoolean value);
  CubitBoolean auto_update_flag();
  void auto_write_flag(CubitBoolean value);
  CubitBoolean auto_write_flag();
  void auto_read_flag(CubitBoolean value);
  CubitBoolean auto_read_flag();
    //- get/set actuate, update, write, and read flags, non-static versions

  CubitBoolean actuate_in_constructor();
    //- return value of actuateInConstructor for this attribute type
  
  CubitBoolean actuate_after_geom_changes();
    //- return value of actuateAfterGeomChanges for this attribute type

//  static CubitBoolean overwrite_flag() {return overwriteFlag;}
//  static void overwrite_flag(CubitBoolean flag) {overwriteFlag = flag;}
    //- get/set overwriteFlag
  
  virtual CubitSimpleAttrib* cubit_simple_attrib() = 0;
    //- return a cubitSimpleAttrib for this CA

  virtual CubitStatus actuate() = 0;
    //- actuate this attrib
  
  virtual CubitStatus actuate_list(DLIList<RefEntity*>);
    //- actuate this attrib on each of the entities passed in

  virtual CubitStatus update() = 0;
    //- update this attrib

  virtual CubitStatus reset() = 0;
    //- reset any lists or other info in this attribute
  
  virtual CubitSimpleAttrib *split_owner();
    //- split this attrib; pass back a new simple attrib if desired

  virtual void merge_owner(CubitAttrib *deletable_ca_ptr);
    //- merge this attrib with deletable_ca_ptr (keep this)

  virtual void transf_owner(const CubitVector &matrow1,
                            const CubitVector &matrow2,
                            const CubitVector &matrow3,
                            const CubitVector &translate_vec,
                            const double scale_factor);
    //- transform this attrib with the data passed in

    virtual int int_attrib_type() = 0;
    //- return the enumerated type of attribute

  const char *att_internal_name() {
    return CGMApp::instance()->attrib_manager()->att_internal_name( int_attrib_type() );
  }
    //- return the internal name of this attribute
  
//  static const char *att_internal_name(CubitAttributeType attrib_type) {
//    return (attrib_type < 0 || attrib_type > CA_LAST_CA)
//          ? NULL : attInternalNames[attrib_type]; 
//  } 
    //- return the internal name of this CA given the enumerated attribute type
  
//  static const char *att_name(CubitAttributeType attrib_type) {
//    return (attrib_type < 0 || attrib_type > CA_LAST_CA)
//          ? NULL : attTypeNames[attrib_type]; 
//  } 
    //- return the internal name of this CA given the enumerated attribute type
  
//  static CubitAttributeType attrib_type(const char* name);
    //- Converts a string to the corresponding int

//  static CubitAttributeType attrib_type_from_internal_name(const char* name);
    //- Converts an "internal name" string to the corresponding int

//  static CubitAttributeType attrib_type(CubitSimpleAttrib *csa_ptr);
    //- returns the type given a simple attribute

//HEADER- RTTI and safe casting functions.
  virtual const type_info& entity_type_info() const
     { return typeid(CubitAttrib);}
  //R- The geometric modeler type
  //- This function returns the type of the geometric modeler.

//  static int equivalent(CubitAttrib* first_attrib_ptr,
//                        CubitAttrib* second_attrib_ptr);
    //- return true if the two ca's are equivalent
  
  int equivalent(CubitSimpleAttrib* csa_ptr);
    //- return true if the csa and this are equivalent
  
//  static CubitBoolean attrib_exists(int attrib_type, CubitBoolean print = CUBIT_TRUE);
    //- find if the attrib exists, and print some details  

  virtual void print();
    //- print some details about this attrib

//  static CubitStatus import_actuate(DLIList<RefEntity*> &entity_list);
    //- given a Body list, actuates first the merge partner attribute
    //- on all entities in the list, then actuates all other types of
    //- attributes

//  static void set_cubit_attrib_factory(CubitAttribFactory *factory)
//    {cubitAttribFactory = factory;}
    //- user-provided attribute factory

//  static int number_of_attributes();
    //- Returns number of defined attributes
  
//  static void clear_attrib_importeds();
    //- clear the attribImported array

//  static void report_attrib_importeds();
    //- print info on the attribImported array

protected:

  virtual void remove_attribute();
  //- Removes 'this' attribute
  virtual void add_attribute();
  //- and adds a new one with the changed data.

  CubitBoolean hasActuated;
    //- flag telling whether this CA has actuated

  CubitBoolean hasUpdated;
    //- flag telling whether this CA has updated

  CubitBoolean hasWritten;
    //- flag telling whether this CA has a corresponding geometry attribute
    //- (on the geometric entity assoc'd with the owner)
  
  CubitBoolean deleteAttrib;
    //- flag telling whether this CA should be deleted

  RefEntity* attribOwnerEntity;
    //- ref entity to which this CA is associated
  
  CubitAttrib* nextAttrib;
    //- next attribute in the chain
  
//   static const char* attTypeNames[];
    //- names of attribute types
  
//   static const char* attInternalNames[];
    //- names of attributes, used in acis files; can't be the same as
    //- attTypeNames because those have to be user-recognizable and can
    //- be multiple words (e.g. "entity name")

//  static CubitBoolean autoActuateFlag[CA_LAST_CA];
//  static CubitBoolean autoUpdateFlag[CA_LAST_CA];
//  static CubitBoolean autoWriteFlag[CA_LAST_CA];
//  static CubitBoolean autoReadFlag[CA_LAST_CA];
//  static CubitBoolean actuateInConstructor[CA_LAST_CA];
//  static CubitBoolean actuateAfterGeomChanges[CA_LAST_CA];
    //- auto flags for each type of CA

//  static bool attribImported[CA_LAST_CA];
    //- static flag used in import type reporting

//  static CubitBoolean overwriteFlag;
    //- if CUBIT_TRUE, any attribute can overwrite default files
    //- (used for mesh container attribute default filenames)

private:

//  static CubitAttribFactory *cubitAttribFactory;
    //- application-provided attribute factory
  
};

// static versions 
//inline CubitBoolean CubitAttrib::auto_actuate_flag(int attrib_type)
//{return autoActuateFlag[attrib_type];}

//inline CubitBoolean CubitAttrib::auto_update_flag(int attrib_type)
//{return autoUpdateFlag[attrib_type];}

//inline CubitBoolean CubitAttrib::auto_write_flag(int attrib_type)
//{return autoWriteFlag[attrib_type];}

//inline CubitBoolean CubitAttrib::auto_read_flag(int attrib_type)
//{return autoReadFlag[attrib_type];}

//inline CubitBoolean* CubitAttrib::auto_read_flag_ptr(int attrib_type)
//{return &autoReadFlag[attrib_type];}

//inline CubitBoolean* CubitAttrib::auto_actuate_flag_ptr(int attrib_type)
//{return &autoActuateFlag[attrib_type];}

//inline CubitBoolean* CubitAttrib::auto_update_flag_ptr(int attrib_type)
//{return &autoUpdateFlag[attrib_type];}

//inline CubitBoolean* CubitAttrib::auto_write_flag_ptr(int attrib_type)
//{return &autoWriteFlag[attrib_type];}

// non-static versions 
inline void CubitAttrib::auto_actuate_flag(CubitBoolean value)
{
  CGMApp::instance()->attrib_manager()->set_auto_actuate_flag(int_attrib_type(), value);
}

inline CubitBoolean CubitAttrib::auto_actuate_flag()
{
  return CGMApp::instance()->attrib_manager()->auto_actuate_flag(int_attrib_type());
}
  
inline CubitBoolean CubitAttrib::actuate_in_constructor()
{
  return CGMApp::instance()->attrib_manager()->actuate_in_constructor(int_attrib_type());
}
    //- get/set actuate, update, write, and read flags, non-static versions

inline CubitBoolean CubitAttrib::actuate_after_geom_changes()
{
  return CGMApp::instance()->attrib_manager()->actuate_after_geom_changes(int_attrib_type());
}
    //- get/set actuate, update, write, and read flags, non-static versions

inline CubitStatus CubitAttrib::set_next_attrib(CubitAttrib* next_attrib_ptr)
{
  nextAttrib = next_attrib_ptr;
  return CUBIT_SUCCESS;
}

inline CubitSimpleAttrib *CubitAttrib::split_owner()
{
    //- split this attrib; pass back a new simple attrib if desired

    // by default, get rid of this attribute and don't copy to new entity
  deleteAttrib = CUBIT_TRUE;
  return NULL;
}

inline void CubitAttrib::merge_owner(CubitAttrib *)
{
    //- merge this attrib with deletable_ca_ptr (keep this)

    // by default, get rid of this attribute
  deleteAttrib = CUBIT_TRUE;
  return;
}

inline void CubitAttrib::transf_owner(const CubitVector &,
                                      const CubitVector &,
                                      const CubitVector &,
                                      const CubitVector &,
                                      const double)
{
    //- by default, do nothing
}

//inline int CubitAttrib::number_of_attributes()
//{return CA_LAST_CA;}

#endif

