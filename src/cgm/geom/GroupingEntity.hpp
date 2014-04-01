//-------------------------------------------------------------------------
// Filename      : GroupingEntity.hpp
//
// Purpose       : This file contains the declarations of the class 
//                 GroupingEntity.
//                 This class is the base class of all the grouping entities
//                 Body, Shell, Loop, Chain.
//
// Special Notes : Each GroupingEntity is associated with a set of 
//                 SenseEntity's (SE's). These SE's are ordered in a 
//                 list.  Hence the GroupingEntity interface
//                 not only provides the ability to get the entire list
//                 of SE's, but also allows you to ask for the "first"
//                 associated SE.  The SE interface thus provides a
//                 function to ask for the "next" SE.  The linked
//                 list of SE's ends when the next function returns a 
//                 NULL pointer. 
//
//                 The same is true of GroupingEntities themselves. They
//                 are an ordered list in the DAG, implying that a next()
//                 function needs to be provided.
//
//                 This is a pure virtual class.
//
// Creator       : Xuechen Liu 
//
// Creation Date : 07/11/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef GROUPING_ENTITY_HPP
#define GROUPING_ENTITY_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "TopologyEntity.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN MACROS DEFINITIONS     **********
// ********** END MACROS DEFINITIONS       **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class SenseEntity;
class BasicTopologyEntity;

// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT GroupingEntity : public TopologyEntity
{
   public :
   
   inline GroupingEntity() ;
     //- The default constructor
   
   virtual ~GroupingEntity();

   inline SenseEntity* get_first_sense_entity_ptr();
   inline SenseEntity* get_last_sense_entity_ptr();
     //R SenseEntity*
     //R- A pointer to the first of the list of SenseEntitys pointed
     //R- to by this object.
     //- This function returns a pointer to the first of the list of
     //- SenseEntity's pointed to by this GroupingEntity.
   
   CubitStatus get_sense_entity_list(DLIList<SenseEntity*>& list);
     //R CubitStatus
     //R- CUBIT_SUCCESS/CUBIT_FAILURE
     //O list
     //O- A list of SenseEntity pointers.
     //- This function returns a list of SenseEntity pointers
     //- associated with this GroupingEntity. If there are none, then
     //- CUBIT_FAILURE is returned.
   
   CubitStatus set_sense_entity_list(DLIList<SenseEntity*>& list,
                                     DLIList<SenseEntity*>& removed);
     //R CubitStatus
     //R- CUBIT_SUCCESS/CUBIT_FAILURE
     //I list
     //I- New set of child SenseEntitys for this GroupingEntity.
     //O removed
     //O- List of child SenseEntitys removed from this GroupingEntity
     //- Change/reorder child list.
     //- It is an error for any input SenseEntity to have a parent 
     //- GroupingEntity unless that parent GroupingEntity is this.
   
   
   inline GroupingEntity* next();
     //R GroupingEntity*
     //R- A GroupingEntity pointer
     //- This function returns a pointer to the GroupingEntity that is
     //- a child of this one in the DAG. If it is the end of the line,
     //- then a NULL pointer is returned.
     
    inline GroupingEntity* previous();
     //R GroupingEntity*
     //R- A GroupingEntity pointer
     //- This function returns a pointer to the GroupingEntity that is
     //- a child of this one in the DAG. If it is the end of the line,
     //- then a NULL pointer is returned.
   
    inline BasicTopologyEntity* get_basic_topology_entity_ptr();
     //R BasicTopologyEntity*
     //R- A pointer to the BasicTopologyEntity which the current grouping
     //R- entity is associated with.
     //- This function returns a pointer to the BasicTopologyEntity which
     //- the current sense entity is associated with.
		 
   CubitStatus add_sense_entity(SenseEntity *sense_entity_ptr,
                                SenseEntity *after_this = 0) ;
     //R CubitStatus
     //I senseEntityPtr
     //I- The pointer to a SenseEntity which will be added to the 
     //I- list of sense entities of the grouping entity.
     //- This function is used to add a SenseEntity to the list of
     //- sense entities of a grouping entity. If the input sense entity
     //- is not the appropriate type, nothing is done and the function
     //- returns CUBIT_FAILURE. If the sense entity is added successfully,
     //- the function returns CUBIT_SUCCESS.
     
   CubitStatus remove_sense_entity(SenseEntity* sense_entity_ptr);
   
   void reverse_direction();  
     // reverse order and sense of child sense entities

   virtual int get_parents( DLIList<ModelEntity*>* list = 0 ) const;
   virtual int get_children(DLIList<ModelEntity*>* list = 0 ) const;

protected :
   
   virtual CubitBoolean query_append_parents( DLIList<ModelEntity*>& list );
   virtual CubitBoolean query_append_children(DLIList<ModelEntity*>& list );

   virtual CubitStatus remove_child_link( ModelEntity* child_ptr );

   CubitStatus disconnect_all_children( DLIList<ModelEntity*>* children = 0 );
   CubitStatus disconnect_all_parents( DLIList<ModelEntity*>* parents = 0 );
   
private :
   
    // for use by BasicTopologyEntity only
  friend class BasicTopologyEntity;
  inline CubitStatus remove_from_list();
  inline void set_basic_topology_entity_ptr( BasicTopologyEntity* );
  inline CubitStatus insert_after(GroupingEntity* next);
  inline CubitStatus insert_before(GroupingEntity* prev);

  BasicTopologyEntity* myParent;
  GroupingEntity* nextInParent;
  GroupingEntity* prevInParent;
  SenseEntity* firstSenseEntity;
  SenseEntity* lastSenseEntity;
  
  GroupingEntity( const GroupingEntity& );
  void operator=( const GroupingEntity& );
};

// ********** BEGIN HELPER CLASSES         **********
// ********** END   HELPER CLASSES         **********

// ********** BEGIN INLINE FUNCTIONS       **********

inline GroupingEntity::GroupingEntity()
  : myParent(0), 
    nextInParent(0),
    prevInParent(0),
    firstSenseEntity(0),
    lastSenseEntity(0)
  {}

inline SenseEntity* GroupingEntity::get_first_sense_entity_ptr()
  { return firstSenseEntity; }

inline SenseEntity* GroupingEntity::get_last_sense_entity_ptr()
  { return lastSenseEntity; }
   
inline GroupingEntity* GroupingEntity::next()
  { return nextInParent; }
   
inline GroupingEntity* GroupingEntity::previous()
  { return prevInParent; }
   
inline BasicTopologyEntity* GroupingEntity::get_basic_topology_entity_ptr()
  { return myParent; }
    
inline CubitStatus GroupingEntity::remove_from_list()
{
  if (nextInParent)
    nextInParent->prevInParent = prevInParent;
  if (prevInParent)
    prevInParent->nextInParent = nextInParent;
  prevInParent = nextInParent = 0;
  return CUBIT_SUCCESS;
}

inline CubitStatus GroupingEntity::insert_after( GroupingEntity* next_ptr )
{
  prevInParent = next_ptr;
  nextInParent = next_ptr->nextInParent;
  if (nextInParent)
    nextInParent->prevInParent = this;
  next_ptr->nextInParent = this;
  return CUBIT_SUCCESS;
}

inline CubitStatus GroupingEntity::insert_before( GroupingEntity* prev_ptr )
{
  nextInParent = prev_ptr;
  prevInParent = prev_ptr->prevInParent;
  if (prevInParent)
    prevInParent->nextInParent = this;
  prev_ptr->prevInParent = this;
  return CUBIT_SUCCESS;
}

inline void GroupingEntity::set_basic_topology_entity_ptr( BasicTopologyEntity* bte_ptr )
{
  assert(!myParent || !bte_ptr);
  myParent = bte_ptr;
}

 
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

