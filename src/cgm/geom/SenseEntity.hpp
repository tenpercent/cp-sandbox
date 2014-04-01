//-------------------------------------------------------------------------
// Filename      : SenseEntity.hpp
//
// Purpose       : This class is the base class of all sense entities such,
//                 as CoVolume, CoFace, CoEdge, and CoVertex.
//
// Special Notes : This is a pure virtual class.
//
// Creator       : Xuechen Liu 
//
// Creation Date : 07/11/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef SENSE_ENTITY_HPP
#define SENSE_ENTITY_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "TopologyEntity.hpp"
#include "CubitObservable.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN MACROS DEFINITIONS     **********
// ********** END MACROS DEFINITIONS       **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class BasicTopologyEntity ;
class GroupingEntity;
template <class X> class DLIList ;
// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT  SenseEntity : public TopologyEntity,
                    public CubitObservable
{
  public :
  
  inline SenseEntity();
   
  virtual ~SenseEntity();

  inline BasicTopologyEntity* get_basic_topology_entity_ptr() ;
     //R BasicTopologyEntity*
     //R- A pointer to the BasicTopologyEntity which the current sense
     //R- entity is associated with.
     //- This function returns a pointer to the BasicTopologyEntity which
     //- the current sense entity is associated with.
		 
	inline GroupingEntity* get_grouping_entity_ptr();
	  //R GroupingEntity*
		//R- A pointer to the parent GroupingEntity.

  BasicTopologyEntity* get_parent_basic_topology_entity_ptr();
   
  CubitStatus attach_basic_topology_entity(BasicTopologyEntity* basicTopologyEntityPtr) ;
     //R CubitStatus
     //R- CUBIT_SUCCESS/FAILURE
     //I basicTopologyEntityPtr
     //I- The pointer to a BasicTopologyEntity which will be attached as
     //I- a BasicTopologyEntity of this SenseEntity.
     //- This function is used to attach a BasicTopologyEntity to this 
     //- SenseEntity. The operation can fail if the BTE is not the
     //- appropriate type. In that case the function returns
     //- CUBIT_FAILURE. If the attachment is done successfully, the function
     //- returns CUBIT_SUCCESS.
     
   CubitStatus switch_basic_topology_entity(BasicTopologyEntity* new_bte );
    //- Change the BasicTopologyEntity attached to this SenseEntity
    //- to the passed BTE.  This is used by unmerge.
   
   inline void set_sense(CubitSense sense) ;
     //R void 
     //I sense
     //I- The sense to be set.
     //- This function sets the sense of the SenseEntity.
     
   void reverse_sense();
     //- Reverse the sense of this SenseEntity.
   
   inline CubitSense get_sense() const ;
     //R CubitSense
     //R- The sense of the sense entity - CUBIT_FORWARD/REVERSED.
     //- This function returns the sense of the SenseEntity.
   
   inline SenseEntity* next();
     //R SenseEntity*
     //R- A SenseEntity pointer
     //- This function returns a pointer to the SenseEntity that is
     //- a child of this one in the DAG. If it is the end of the line,
     //- then a NULL pointer is returned.
   
   inline SenseEntity* previous();
     //R SenseEntity*
     //R- A SenseEntity pointer
     //- This function returns a pointer to the SenseEntity that is
     //- a "parent" of this one in the DAG. If this SenseEntity is
     //- already at the beginning of the list, then a NULL pointer is 
     //- returned.
   
   inline SenseEntity* next_on_bte();
    //R SenseEntity*
    //R- The next sense entity in the linked of list of 
    //R- sense entities associated with the child bte.
   
   virtual int get_parents( DLIList<ModelEntity*>* list = 0 ) const;
   virtual int get_children( DLIList<ModelEntity*>* list = 0 ) const;
   
   
  protected :
   
   virtual CubitBoolean query_append_parents( DLIList<ModelEntity*>& list );
   virtual CubitBoolean query_append_children(DLIList<ModelEntity*>& list );
   
   virtual CubitStatus remove_child_link( ModelEntity* );
   
   CubitStatus disconnect_all_children( DLIList<ModelEntity*>* children = 0 );
   CubitStatus disconnect_all_parents( DLIList<ModelEntity*>* parents = 0 );
   
   private :  
   
    // functions for use by GroupingEntity only
   friend class GroupingEntity;
   inline CubitStatus gpe_insert_after (SenseEntity* prev_ptr);
   inline CubitStatus gpe_insert_before(SenseEntity* next_ptr);
   inline CubitStatus gpe_remove();
   inline CubitStatus set_grouping_entity_ptr(GroupingEntity* gpe_ptr);
   inline void swap_gpe_list_ptrs();
   
    // functions for use by BasicTopologyEntity only
   friend class BasicTopologyEntity;
   inline SenseEntity* set_bte_next(SenseEntity* next_ptr);
   inline CubitStatus set_basic_topology_entity_ptr(BasicTopologyEntity* bte_ptr);
   
   CubitSense mySense;
   
   GroupingEntity* myParent;
   SenseEntity* nextInParent;
   SenseEntity* prevInParent;
   
   BasicTopologyEntity* myChild;
   SenseEntity* nextInChild;
   
   SenseEntity( const SenseEntity& );
   void operator=( const SenseEntity& );
};


SenseEntity::SenseEntity()
  : mySense(CUBIT_FORWARD),
    myParent(0), 
    nextInParent(0),
    prevInParent(0),
    myChild(0),
    nextInChild(0)
  {}

BasicTopologyEntity* SenseEntity::get_basic_topology_entity_ptr() 
  { return myChild; }

GroupingEntity* SenseEntity::get_grouping_entity_ptr()
  { return myParent; }

void SenseEntity::set_sense(CubitSense sense)
  { mySense = sense; }
  
CubitSense SenseEntity::get_sense() const 
  { return mySense; }

SenseEntity* SenseEntity::next()
  { return nextInParent; }

SenseEntity* SenseEntity::previous()
  { return prevInParent; }

SenseEntity* SenseEntity::next_on_bte()
  { return nextInChild; }
   
CubitStatus SenseEntity::gpe_insert_after(SenseEntity* prev_ptr)
{
  nextInParent = prev_ptr->nextInParent;
  if (nextInParent)
    nextInParent->prevInParent = this;
  prevInParent = prev_ptr;
  prev_ptr->nextInParent = this;
  return CUBIT_SUCCESS;
}  
   
CubitStatus SenseEntity::gpe_insert_before(SenseEntity* next_ptr)
{
  prevInParent = prevInParent;
  if (prevInParent)
    prevInParent->nextInParent = this;
  nextInParent = next_ptr;
  next_ptr->prevInParent = this;
  return CUBIT_SUCCESS;
}  

CubitStatus SenseEntity::gpe_remove()
{
  if (nextInParent)
    nextInParent->prevInParent = prevInParent;
  if (prevInParent)
    prevInParent->nextInParent = nextInParent;
  prevInParent = nextInParent = 0;
  return CUBIT_SUCCESS;
}

void SenseEntity::swap_gpe_list_ptrs()
{
  SenseEntity* tmp = nextInParent;
  nextInParent = prevInParent;
  prevInParent = tmp;
}

CubitStatus SenseEntity::set_grouping_entity_ptr(GroupingEntity* gpe_ptr)
{
  myParent = gpe_ptr;
  return CUBIT_SUCCESS;
}

SenseEntity* SenseEntity::set_bte_next(SenseEntity* next_ptr)
{
  SenseEntity* old_val = nextInChild;
  nextInChild = next_ptr;
  return old_val;
}

CubitStatus SenseEntity::set_basic_topology_entity_ptr(
                                  BasicTopologyEntity* bte_ptr)
{
  if (myChild && bte_ptr) 
    { assert(0); return CUBIT_FAILURE; }
    
  myChild = bte_ptr;
  return CUBIT_SUCCESS;
}


// ********** BEGIN HELPER CLASSES         **********
// ********** END   HELPER CLASSES         **********

// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********


#endif


