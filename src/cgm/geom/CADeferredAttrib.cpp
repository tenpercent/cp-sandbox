#include "CubitSimpleAttrib.hpp"
#include "ToolDataUser.hpp"
#include "TDUniqueId.hpp"
#include "CastTo.hpp"
#include "RefEntity.hpp"
#include "SDLCADeferredAttribList.hpp"
#include "DLIList.hpp"

DLIList<CubitAttrib*> CADeferredAttrib::unactuatedCAs;
SDLCADeferredAttribList CADeferredAttrib::allCADeferredAttribs;

CubitAttrib* CADeferredAttrib_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CADeferredAttrib *new_attrib = NULL;

  // Deferred attributes are a special case -- they should only get created
  // when restoring from a file, in which case p_csa should be non-NULL
  if (NULL != p_csa)
  {
    new_attrib = new CADeferredAttrib(entity, p_csa);
  }
  return new_attrib;
}


CADeferredAttrib::~CADeferredAttrib()
{
  while (allCADeferredAttribs.move_to(this))
    allCADeferredAttribs.remove();

  delete thisCSA;
}

CADeferredAttrib::CADeferredAttrib(RefEntity *, CubitSimpleAttrib *csa_ptr) 
        : CubitAttrib(NULL)
{
  init_csa(csa_ptr);
  hasUpdated = CUBIT_TRUE;
}

CubitStatus CADeferredAttrib::init_csa(CubitSimpleAttrib *csa_ptr) 
{
  int csa_type = CGMApp::instance()->attrib_manager()->attrib_type(csa_ptr);
  if (CA_DEFERRED_ATTRIB != csa_type)
  {
    assert(false);
    return CUBIT_FAILURE;
  }

    // initialize this according to csa_ptr

    // first get the uniqueId off the csa
  csa_ptr->int_data_list()->reset();
  uniqueId = *csa_ptr->int_data_list()->get();
  assert(uniqueId > 0);

    // now check to see if we have this CADA alreedy; if so, set the
    // delete flag and exit
  CubitBoolean already_have = CUBIT_FALSE;
  if (allCADeferredAttribs.move_to(uniqueId)) {
    do {
      if (allCADeferredAttribs.get_and_step()->equivalent(csa_ptr)) {
        already_have = CUBIT_TRUE;
        break;
      }
    }
    while (!allCADeferredAttribs.is_at_beginning() &&
           allCADeferredAttribs.get()->unique_id() == uniqueId);
  }

  if (already_have) {
    deleteAttrib = CUBIT_TRUE;
    return CUBIT_SUCCESS;
  }

    // copy the info on the csa; need a new one, since
    // we don't own the original
  thisCSA = csa_from_dcsa(csa_ptr);
  assert(thisCSA != 0);

    // add this to the global list
  allCADeferredAttribs.append(this);

  return CUBIT_SUCCESS;
}

CubitStatus CADeferredAttrib::actuate()
{
    // test to see if we can assign this CADA to a new owner; set
    // flag accordingly
  if (assign_to_owner()) {
    hasActuated = CUBIT_TRUE;
    return CUBIT_SUCCESS;
  }
  else return CUBIT_FAILURE;
}

CubitStatus CADeferredAttrib::update()
{
    // the hasUpdated flag should always be true by the time
    // we get here, since the function that put this CA on the
    // owning entity set the flag
  assert(hasUpdated == CUBIT_TRUE);
  return CUBIT_SUCCESS;
}

CubitStatus CADeferredAttrib::reset()
{
    // do nothing; this CA manages its own duplicates, so no need
    // to worry about them in parent code
  return CUBIT_SUCCESS;
}

CubitSimpleAttrib* CADeferredAttrib::cubit_simple_attrib()
{
  return new CubitSimpleAttrib(thisCSA);
}

CubitStatus CADeferredAttrib::assign_to_owner(CubitAttribUser *owner)
{
    //- looks for an entity with the right uid, assigns itself to
    //- that entity if found

  if (uniqueId == 0) return CUBIT_FAILURE;
  
  if (owner == NULL) {
    ToolDataUser *tdu = TDUniqueId::find_td_unique_id(uniqueId);
    owner = CAST_TO(tdu, CubitAttribUser);
  }

  if (owner == NULL) return CUBIT_FAILURE;

    // ok, we have an owner; create a new CA using the csa, assigning
    // it to the ref entity
  RefEntity *ref_ent = CAST_TO(owner, RefEntity);
  assert (ref_ent != 0);

//  CubitAttrib *new_cubit_attrib =
//    CGMApp::instance()->attrib_manager()->create_cubit_attrib(thisCSA, ref_ent);
  int attrib_type = CGMApp::instance()->attrib_manager()->attrib_type(thisCSA);
  CubitAttrib *new_cubit_attrib =
    CGMApp::instance()->attrib_manager()->create_cubit_attrib(attrib_type, ref_ent, thisCSA);

    // now remove this CADA from the global list and add the new CA
    // to the unactuated list
  allCADeferredAttribs.move_to(this);
  allCADeferredAttribs.remove();

    // new attribute might be NULL, if there was already one there
  if (NULL != new_cubit_attrib) unactuatedCAs.append(new_cubit_attrib);

    // ok, all done
  
  deleteAttrib = CUBIT_TRUE;
  return CUBIT_SUCCESS;
}

CubitStatus CADeferredAttrib::get_deferred_attribs(const int uid,
                                                   SDLCADeferredAttribList &cada_list)
{
    // find a deferred attribute for the entity with the uid passed in
  cada_list.clean_out();
  if (allCADeferredAttribs.move_to(uid)) {
    do
      cada_list.append(allCADeferredAttribs.get_and_step());
    while (!allCADeferredAttribs.is_at_beginning() &&
           allCADeferredAttribs.get()->unique_id() == uid);
  }
  
  if (cada_list.size() > 0) return CUBIT_SUCCESS;
  else return CUBIT_FAILURE;
}

CubitStatus CADeferredAttrib::cleanup_cadas(const CubitBoolean from_constructor,
                                            const CubitBoolean after_geom_changes) 
{
  CubitStatus status = CUBIT_FAILURE;
  
  if (CUBIT_TRUE == from_constructor)
    status = CADeferredAttrib::cleanup_cadas_private(CUBIT_TRUE, CUBIT_FALSE);
  
    // exit if we're not to actuate for after_geom_changes
  if (CUBIT_TRUE == after_geom_changes) 
    status = CADeferredAttrib::cleanup_cadas_private(CUBIT_FALSE, CUBIT_TRUE);

  return status;  
}

CubitStatus CADeferredAttrib::cleanup_cadas_private(const CubitBoolean from_constructor,
                                                    const CubitBoolean after_geom_changes) 
{
    // moves between the global CADA list and the unactuated list:
    // 
    // 1. tries to actuate all CADAs on unactuated list
    // 2. tries to assign_to_owner all CADAs on global list
    // 
    // after each of these steps, if anything happened, the loop is
    // repeated
    // this function should be called as part of the auto actuate process

    // first call for from_constructor and !after_geom_changes attributes
  CubitBoolean done = CUBIT_FALSE;
  CubitStatus did_something = CUBIT_FAILURE;

  while (done == CUBIT_FALSE) {
    done = CUBIT_TRUE;
    int i;
    for (i = unactuatedCAs.size(); i > 0; i--) {
      CubitAttrib *attrib = unactuatedCAs.get();
      if( attrib == 0 )
        continue;

        // check the auto actuate flag for this CA; since this function
        // is only called from the auto actuating process, we don't need
        // to check whether the user requested that this attribute be actuated
      if (attrib->auto_actuate_flag() == CUBIT_TRUE &&
          (!from_constructor || attrib->actuate_in_constructor()) &&
          (after_geom_changes || !attrib->actuate_after_geom_changes()))
      {

          // if the attribute has already actuated, but is still in our list,
          // count it as doing something
        if ( /*attrib->has_actuated()  || */
            attrib->actuate() == CUBIT_SUCCESS) {
          PRINT_DEBUG_90("Actuated a CADA of type %s\n",
                         attrib->att_internal_name());
          did_something = CUBIT_SUCCESS;
          done = CUBIT_FALSE;
          unactuatedCAs.remove();
        }
        else unactuatedCAs.step();
      }
      else unactuatedCAs.step();
    }
    
    if (done == CUBIT_TRUE) break;

      // if we did something in the previous loop, some of our CADAs might
      // now have owners; check and see
    for (i = allCADeferredAttribs.size(); i > 0; i--) {
      CADeferredAttrib *cada = allCADeferredAttribs.get_and_step();
      if (cada->assign_to_owner() == CUBIT_SUCCESS) {
        PRINT_DEBUG_90("Assigned a CADA to a new owner in CADA::cleanup_cadas\n");
        did_something = CUBIT_SUCCESS;
        done = CUBIT_FALSE;
      }
    }
  }

    // if there wasn't anthing to do in the first place, we didn't
    // really fail
  if (did_something == CUBIT_FAILURE && unactuatedCAs.size() == 0)
    did_something = CUBIT_SUCCESS;
  
  return did_something;
}

CubitStatus CADeferredAttrib::owner_created(RefEntity *new_owner, const int uid) 
{
    // for a newly created ref entity, assigns any CADA with the same uid to the
    // new entity

    // get any CADAs with the same uid
  SDLCADeferredAttribList attrib_list;
  get_deferred_attribs(uid, attrib_list);

    // now assign them to the owner
  for (int i = attrib_list.size(); i > 0; i--)
    attrib_list.get_and_step()->assign_to_owner(new_owner);

  if (attrib_list.size() > 0) return CUBIT_SUCCESS;
  else return CUBIT_FAILURE;
}
/*
CubitBoolean CADeferredAttrib::is_match(CubitSimpleAttrib *csa_ptr,
                                        const int uid)
{
    //- returns true if the simple attribute is deferred type and matches
    //- uid
  assert(csa_ptr && uid > 0);
  csa_ptr->string_data_list()->reset();
  csa_ptr->string_data_list()->step();
  csa_ptr->int_data_list()->reset();
  if (!strcmp(csa_ptr->string_data_list()->get()->c_str(),
              CGMApp::instance()->attrib_manager()->att_internal_name(CA_DEFERRED_ATTRIB)) &&
      *csa_ptr->int_data_list()->get() == uid)
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}
*/


CubitSimpleAttrib *CADeferredAttrib::csa_from_dcsa(CubitSimpleAttrib *csa_ptr,
                                                   const int uid)
{
  
    //- given a deferred csa, convert it to a normal csa by removing
    //- first type string and first int; if first int doesn't match
    //- uid passed in, NULL is returned
  assert(csa_ptr != 0);

  csa_ptr->string_data_list()->reset();
  if (strcmp(csa_ptr->string_data_list()->get()->c_str(),
             CGMApp::instance()->attrib_manager()->att_internal_name(CA_DEFERRED_ATTRIB)))
      // csa isn't deferred type - return
    return NULL;

  csa_ptr->int_data_list()->reset();
  if (uid != 0 && *csa_ptr->int_data_list()->get() != uid)
      // csa uid doesn't match - return
    return NULL;

    // else we have a match - build new csa
  CubitSimpleAttrib *new_csa = new CubitSimpleAttrib(csa_ptr);
  
  new_csa->string_data_list()->reset();
  delete new_csa->string_data_list()->remove();
  
  new_csa->int_data_list()->reset();
  delete new_csa->int_data_list()->remove();

  return new_csa;
}

bool CADeferredAttrib::add_unactuated_ca(CubitAttrib *ca_ptr) 
{
  return (CUBIT_TRUE == unactuatedCAs.append_unique(ca_ptr) ? true : false);
}

bool CADeferredAttrib::remove_unactuated_ca( CubitAttrib* ca_ptr )
{
  int item_index = unactuatedCAs.where_is_item( ca_ptr );
  if( item_index >= 0 )
  {
      // This is often called when other code is stepping through
      // the list, so don't change the current list position.
    int curr_index = unactuatedCAs.get_index();
    int step = item_index - curr_index;
    unactuatedCAs.step( step );
    assert( unactuatedCAs.get() == ca_ptr );
    unactuatedCAs.change_to(0);
    unactuatedCAs.step( -step );
    assert( unactuatedCAs.get_index() == curr_index );
    return true;
  }
  else
    return false;
}
