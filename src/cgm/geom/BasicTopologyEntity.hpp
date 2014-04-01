//-------------------------------------------------------------------------
// Filename      : BasicTopologyEntity.hpp
//
// Purpose       : This is the interface to the BasicTopologyEntity base class.
//                 The main characteristic of specific BasicTopoEntities is
//                 that they form the basis of a topological description
//                 of a solid model.  Having ordered lists of these would
//                 be a complete and unambiguous description of a solid
//                 model.
//
//                 The other characteristic of BasicTopoEntities is that 
//                 they contain ("own") various mesh-related entities. 
//
// Special Notes : Although each BasicTopologyEntity is associated with a set of 
//                 GroupingEntity's, there is no member data in a BasicTopologyEntity
//                 that stores this list.  Instead, this connection is
//                 established within the Directed Acyclic Graph (DAG)
//                 datastructure that captures the relationships between
//                 the various TopologyEntity's in the Model.  Each
//                 TopologyEntity has a member datum which is a pointer
//                 to a ModelEntity.  This
//                 pointer is a link from the TopologyEntity to the node
//                 in the Model DAG that represents its "position" and 
//                 "links" within the Model. 
//
//                 Each BasicTopologyEntity (with the exception of RefVertex) is
//                 associated with a set of GroupingEntity's (GrE's).  These 
//                 GrE's are ordered in a list.  Hence, the BasicTopologyEntity 
//                 interface not only provides the ability to get the entire 
//                 list of GrE's, but also allows you to ask for the "first"
//                 associated GrE. By extension, the GrE interface, provides 
//                 a function to ask for the "next" GrE.  The linked 
//                 list of GrE's ends when the "next" function returns a 
//                 NULL pointer. 
//
//                 Each BasicTopologyEntity HasA GeometryEntity pointer.  This 
//                 pointer is stored explicitly as a member datum because 
//                 Geometry Entity's are *not* represented by nodes in the 
//                 Model DAG.
//
//                 This is a pure virtual class, to prevent instantiation.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/14/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef BASIC_TOPOLOGY_ENTITY_HPP
#define BASIC_TOPOLOGY_ENTITY_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "TopologyEntity.hpp"
#include "RefEntity.hpp"
#include "CubitBox.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********

template <class X> class DLIList;
class GroupingEntity;
class SenseEntity;
class GeometryEntity;

// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN MACRO DEFINITIONS     **********
// ********** END MACRO DEFINITIONS       **********

// ********** BEGIN ENUM DEFINITIONS     **********
// ********** END ENUM DEFINITIONS       **********

class CUBIT_GEOM_EXPORT BasicTopologyEntity : public TopologyEntity,
                            public RefEntity
{
public:
  static const char* get_class_name()
    { return "BasicTopologyEntity"; }

  virtual const char* class_name() const
    { return get_class_name(); }
  
  virtual DagType dag_type() const = 0;
  
   
  inline BasicTopologyEntity() ;  
    //- Default constructor.
   
  virtual ~BasicTopologyEntity();  

 CubitStatus get_grouping_entity_list(
      DLIList<GroupingEntity*>& groupingEntityList) const;
    //R  CubitStatus
    //R- CUBIT_SUCCESS/FAILURE.
    //O  groupingEntityList
    //O- The list of GroupingEntity pointers associated with this
    //O- BasicTopologyEntity.
    //-  This function returns a list of pointers to GroupingEntity's
    //-  associated with this BasicTopologyEntity.
    
  inline GroupingEntity* get_first_grouping_entity_ptr() const;
    //R GroupingEntity* 
    //R- The child GroupingEntity pointer or NULL if none.
		
	CubitStatus get_sense_entity_list( DLIList<SenseEntity*>& senseEntityList ) const;
	  //R- CubitStatus
		//O senseEntityList
		//O- The parent SenseEntity pointers for this BTE.
  
  inline SenseEntity* get_first_sense_entity_ptr() const;
  
  CubitStatus add_grouping_entity(GroupingEntity*) ;
    //R  CubitStatus
    //R- CUBIT_SUCCESS/FAILURE.
    //I  GroupingEntity*
    //I- A pointer to a GroupingEntity which will be added to
    //I- the list of grouping entities associated with this BasicTopologyEntity.
    //-  This function is used to add a GroupingEntity to the
    //-  list of grouping entities associated with this BasicTopologyEntity.

  CubitStatus remove_grouping_entity(GroupingEntity*);
    //R  CubitStatus
    //R- CUBIT_SUCCESS/FAILURE.
    //I  GroupingEntity*
    //I- A pointer to a GroupingEntity which will be removed from
    //I- the list of grouping entities associated with this BasicTopologyEntity.
    //-  This function is used to remove a GroupingEntity from the
    //-  list of grouping entities associated with this BasicTopologyEntity.
  
  CubitStatus set_grouping_entity_list( DLIList<GroupingEntity*>& new_list,
                                        DLIList<GroupingEntity*>& removed_list );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I new_list
    //I- The list of child grouping entities for this BTE
    //O removed_list
    //O- Any grouping entities disconnected from this BTE
    //- Make the child list of grouping entities for this BTE be
    //- the passed list, in the same order as the passed list.
    //- Pass back any grouping entities that were children of this
    //- BTE but were not in new_list.
  
  CubitStatus add_sense_entity(SenseEntity*);
    //R  CubitStatus
    //R- CUBIT_SUCCESS/FAILURE.
    //I  SenseEntity*
    //I- A pointer to a SenseEntity which will be added to
    //I- the list of sense entities associated with this BasicTopologyEntity.
    //-  This function is used to add a SenseEntity to the
    //-  list of sense entities associated with this BasicTopologyEntity.

  CubitStatus remove_sense_entity(SenseEntity*);
    //R  CubitStatus
    //R- CUBIT_SUCCESS/FAILURE.
    //I  SenseEntity*
    //I- A pointer to a SenseEntity which will be removed from
    //I- the list of sense entities associated with this BasicTopologyEntity.
    //-  This function is used to remove a SenseEntity from the
    //-  list of sense entities associated with this BasicTopologyEntity.
  
  SenseEntity* find_sense_entity(GroupingEntity* gpe) const;
    //R SenseEntity*
    //R-  A parent SenseEntity of this BasicTopologEntity, or NULL.
    //I GroupingEntity*
    //I- A immediate parent grouping entity of this basic topology entity.
    //- Find the sense entity connecting this BTE to the passed 
    //- grouping entity.  Returns NULL if more than one sense entity.
  
  SenseEntity* find_sense_entity(BasicTopologyEntity* bte) const;
    //R SenseEntity*
    //R-  A parent SenseEntity of this BasicTopologEntity, or NULL.
    //I BasicTopologyEntity*
    //I- A immediate parent basic topology entity of this basic 
    //I- topology entity.
    //- Find the sense entity connecting this BTE to the passed 
    //- BTE.  Returns NULL if more than one sense entity.
    
  CubitStatus get_sense_entities( DLIList<SenseEntity*>& result,
                                  GroupingEntity* in_this);
  CubitStatus get_sense_entities( DLIList<SenseEntity*>& result,
                                  BasicTopologyEntity* in_this);
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I in_this
    //I- Subset of parent sense entities joining this to 'in_this'
    //I- Must be an immediate parent.
    //O result
    //O- List of parent sense entities.
    //- Get parent sense entities, optionally limiting the result
    //- to the subset connecting this to a passed entity.
  
  CubitBoolean is_nonmanifold( GroupingEntity* in_this_parent );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I in_this_parent
    //I- Immediate parent grouping entity to test with respect to.
    //- Result is false if there is exactly one sense entity
    //- connecting this entity to the parent grouping entity.  
    //- Result is true if there are more than one sense entities
    //- connecting this entity to the parent grouping entity.
    //- Result is undefined if there are no sense entities connecting
    //- this entity to the passed parent.
   
  GeometryEntity* get_geometry_entity_ptr() const;
    //R GeometryEntity*
    //R- A pointer to the GeometryEntity to which the current 
    //R- BasicTopologyEntity points. 
    //-  This function returns a pointer to the GeometryEntity
    //-  which the current BasicTopologyEntity points to.
   
  GeometryType geometry_type() const;
    //R GeometryType
    //R- An enumerated type describing the underlying geometry
    //R- representation
    //-  This function returns the type of geometry representation
    //- underlying this entity

  virtual CubitBox bounding_box();
    //- Returns the bounding box of this entity
   
  CubitStatus set_geometry_entity_ptr(GeometryEntity* geometryEntityPtr) ;
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I  geometryEntityPtr
    //I- A pointer to the GeometryEntity that will be associated with
    //I- this RefEntity.
    //- This function sets the GeometryEntity associated with this 
    //- RefEntity.
    //- CUBIT_FAILURE is returned if a problem was detected.
   
  double measure();
    //R double
    //R- The numeric value of the measure (its units depend on the dimension
    //R- of the RefEntity being "measured")
    //- A generic geometric extent function.
    //- Returns volume for RefVolumes, area for RefFaces, length for RefEdge,
    //- and 1.0 for RefVertices
   
   virtual int get_parents( DLIList<ModelEntity*>* list = 0 ) const;
   virtual int get_children(DLIList<ModelEntity*>* list = 0 ) const;
   
protected: 
   virtual CubitBoolean query_append_parents( DLIList<ModelEntity*>& list );
   virtual CubitBoolean query_append_children(DLIList<ModelEntity*>& list );
   
   virtual CubitStatus remove_child_link(ModelEntity* entity_ptr);
   
   CubitStatus disconnect_all_children( DLIList<ModelEntity*>* children = 0 );
   CubitStatus disconnect_all_parents( DLIList<ModelEntity*>* parents = 0 );
   
private:

  SenseEntity* firstSenseEntity;
  SenseEntity* lastSenseEntity;
  GroupingEntity* firstGroupingEntity;
  GroupingEntity* lastGroupingEntity;
  
  BasicTopologyEntity( const BasicTopologyEntity& );
  void operator=( const BasicTopologyEntity& );
} ;


// ********** BEGIN INLINE FUNCTIONS       **********
BasicTopologyEntity::BasicTopologyEntity()
  : firstSenseEntity(0),
    lastSenseEntity(0),
    firstGroupingEntity(0),
    lastGroupingEntity(0)
  {}

SenseEntity* BasicTopologyEntity::get_first_sense_entity_ptr() const
  { return firstSenseEntity; }
  
GroupingEntity* BasicTopologyEntity::get_first_grouping_entity_ptr() const
  { return firstGroupingEntity; }

// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif


