//- Class:          CADeferredAttrib
//- Owner:          Tim Tautges
//- Description:    Cubit attribute for storing attributes that will be used later,
//-       i.e. for entities which don't yet exist in the model
//- Checked by:
//- Version:
//- 
//- Deferred attributes are used to store data on one solid model entity
//- which actually apply to a different entity, one which doesn't exist
//- in the solid model as stored (usually something like virtual geometry)
//-
//- This attribute saves/updates something like this:
//- 1) Normal attribute (mesh interval, scheme, etc.) updates
//- 2) Normal attribute writes:
//-    - creates simple attribute, csa_ptr, and passes to 
//-      GeometryEntity::append_simple_attrib
//-    - For virtual entities, this is implemented in 
//-      VirtualEntity::append_simple_attrib_virt
//- 3) Inside VE::append_simple_attrib_virt, the simple attribute is combined
//-    with a uid for the owning ref entity and the two are used to make a
//-    new *simple* deferred attribute
//- 4) The simple deferred attribute is put onto the first underlying entity
//-    using TopologyBridge::append_simple_attrib
//- 
//- To restore:
//- 1) Take dcsa, create CADeferredAttrib from it, and put it on owning entity:
//-    - CADA stores uid from dcsa
//-    - Then the dcsa is converted to a csa (CADA stuff removed) and stored 
//-      in CADA on owning entity
//- 2) When virtual geometry is created from actuating an attribute, CADA is notified;
//-    CADA goes through the list of CADA's, looking for any that apply to the
//-    new entity (based on uid)
//- 3) If any CADA's are found for the new entity, their csa data are converted to
//-    a real CA, and that CA is put on the unactuatedCAs list in CADA.
//- 4) At the end of CubitAttribUser::auto_actuate_cubit_attrib, CADA::cleanup_cadas
//-    is called, which keeps going through the unactuatedCAs list until nothing
//-    new happens (no new unactuated cas)
//-

#ifndef CA_DEFERRED_ATTRIB_HPP
#define CA_DEFERRED_ATTRIB_HPP

#include "CubitAttrib.hpp"
#include "DLIList.hpp"

#include "CADefines.hpp"

class RefEntity;
class CubitSimpleAttrib;
class SDLCADeferredAttribList;

class CUBIT_GEOM_EXPORT CADeferredAttrib: public CubitAttrib
{
private:

  int uniqueId;
    //- the deferred id tag for an entity

  CubitSimpleAttrib *thisCSA;
    //- the deferred attribute information

  static SDLCADeferredAttribList allCADeferredAttribs;
    //- list of all CADA's; used in actuate_all function

  static DLIList<CubitAttrib*> unactuatedCAs;
    //- list of new CubitAttrib's that haven't yet been actuated

  static CubitStatus cleanup_cadas_private(const CubitBoolean from_constructor,
                                           const CubitBoolean after_geom_changes);
    //- moves between the global CADA list and the unactuated list:
    //- 
    //- 1. tries to actuate all CADAs on unactuated list
    //- 2. tries to assign_to_owner all CADAs on global list
public:
  virtual ~CADeferredAttrib();

//  CADeferredAttrib(RefEntity*);

  CADeferredAttrib(RefEntity*, CubitSimpleAttrib*);
    //- make a CADA from a simple attribute

  //HEADER- RTTI and safe casting functions.
  virtual const type_info& entity_type_info() const
     { return typeid(CADeferredAttrib);}
  //R- The geometric modeler type
  //- This function returns the type of the geometric modeler.

 CubitStatus actuate();

  CubitStatus update();

  CubitStatus reset();

  CubitSimpleAttrib* cubit_simple_attrib();

  int unique_id() { return uniqueId;}

  void unique_id (int id) {uniqueId = id;}
private:
//  CubitSimpleAttrib *this_csa() {return thisCSA;};

  CubitStatus init_csa(CubitSimpleAttrib *csa_ptr);

  int int_attrib_type() {return CA_DEFERRED_ATTRIB;};

  CubitStatus assign_to_owner(CubitAttribUser *owner = NULL);
    //- looks for an entity with the right uid, assigns itself to
    //- that entity if found

  static CubitStatus get_deferred_attribs(const int uid,
                                          SDLCADeferredAttribList &cada_list);
    //- given a uid, return a list of CADAs with corresponding uid

public:
  static DLIList<CubitAttrib*> get_unactuated_deferred_attribs() { return unactuatedCAs; }
    //- get all unactuated deferred attribs

  static CubitStatus cleanup_cadas(const CubitBoolean from_constructor,
                                   const CubitBoolean after_geom_changes);
    //- moves between the global CADA list and the unactuated list:
    //- 
    //- 1. tries to actuate all CADAs on unactuated list
    //- 2. tries to assign_to_owner all CADAs on global list

  static CubitStatus owner_created(RefEntity *new_owner, const int uid);
    //- for a newly created ref entity, assigns any CADA with the same uid to the
    //- new entity
private:
//  static CubitBoolean is_match(CubitSimpleAttrib *csa_ptr, const int uid);
    //- returns true if the simple attribute is deferred type and matches
    //- uid

#ifdef BOYD14
  static CubitSimpleAttrib *dcsa_from_csa(CubitSimpleAttrib *csa_ptr,
                                          const int uid);
    //- given a csa, make a deferred cubit simple attribute by changing
    //- the type (inserts type at front of string list) and inserting
    //- uid
#endif
  
  static CubitSimpleAttrib *csa_from_dcsa(CubitSimpleAttrib *csa_ptr,
                                          const int uid = 0);
    //- given a deferred csa, convert it to a normal csa by removing
    //- first type string and first int; if first int doesn't match
    //- uid passed in, NULL is returned

public:
  static CubitBoolean add_unactuated_ca(CubitAttrib *ca_ptr);
    //- adds an unactuated ca (usually one that didn't actuate because it
    //- depends on some piece of geometry that doesn't exist yet)

  static CubitBoolean remove_unactuated_ca(CubitAttrib* ca_ptr);
    //- If attribute gets destroted before it actuates, remove it
    //- from the list.
};

CubitAttrib* CADeferredAttrib_creator(RefEntity* entity, CubitSimpleAttrib *p_csa);

#endif

