//-------------------------------------------------------------------------
// Filename      : RefEntity.hpp
//
// Purpose       : This file consolidates functionality from several
//                 base classes in CUBIT, and serves as the entry point
//                 to geometry from the meshing and other CUBIT classes
//
// Creator       : Tim Tautges (in it's new state)
//
// Creation Date : 9/99
//-------------------------------------------------------------------------

#ifndef REFENTITY_HPP
#define REFENTITY_HPP

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitObservable.hpp"
#include "CubitEntity.hpp"
#include "CubitAttribUser.hpp"
#include "CubitColorConstants.hpp"
#include "ToolDataUser.hpp"
#include "CubitEventDefines.h"
#include "DagType.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********

class CubitString;
template <class X> class DLIList;
class RefVolume;
class RefFace;
class RefEdge;
class RefVertex;
class RefEntity;
class MeshEntity;
class CubitEntity;

class Body;
class RefGroup;
class RefEdge;
class RefFace;
class RefVolume;
class RefVertex;

// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN MACRO DEFINITIONS     **********
// ********** END MACRO DEFINITIONS       **********

// ********** BEGIN ENUM DEFINITIONS     **********

enum AutoMergeStatus {
  AUTO_MERGE_AUTO = 2,
  AUTO_MERGE_ON   = 1,
  AUTO_MERGE_OFF  = 0
};


// ********** END ENUM DEFINITIONS       **********

class CUBIT_GEOM_EXPORT RefEntity: public CubitEntity,
                 public CubitObservable,
                 public ToolDataUser,
                 public CubitAttribUser
{
public:
  
  RefEntity() ;
    //- The default constructor
  
  virtual ~RefEntity ();
    //- A pure virtual destructor - this ensures that the class is
    //- non-instantiable.

    // Class name related functions
  static const char* get_ref_class_name(const type_info& ref_type);
  
    // Name-related functions
  static RefEntity* get_by_name(const CubitString& name);
  CubitString entity_name() const;
  CubitStatus entity_name (CubitString name);
  void entity_names(DLIList<CubitString*> &names) const;
  int num_names() const;
  CubitStatus generate_default_name ( CubitString &name );
  CubitStatus assign_default_name( CubitBoolean user_setting = CUBIT_FALSE );
  CubitStatus remove_entity_name(CubitString const & name);
  CubitStatus remove_entity_names();
  void merge_entity_names(RefEntity* dead_entity);
  void switch_entity_names(RefEntity *other_entity);
    //- Functions related to naming of RefEntity objects
  
  virtual void marked(int value);
  virtual int  marked ();
    //- Generic flag for temporary use by algorithms.
    //- Value is volatile and may change unexpectedly.

  virtual void is_mergeable(AutoMergeStatus val);
  AutoMergeStatus merge_status() const;
  bool is_mergeable();
  CubitBoolean is_merged();
  // return true if this RefEntity has multiple TopologyBridges.

  void update_auto_merge_state();
    //- Set/Get the mergeable status for this entity.
    
  virtual bool children_mergeable();
    //- Get whether all child entities are mergeable.
  
  virtual int can_modify();
    //- Allow unmerging and other operations.  Default in 
    //- RefEntity always returns true.  Provided for
    //- derived classes to override.
  
  virtual int dimension() const; 
    // Returns the geometric dimension of the entity. 
    // vertex == 0, edge == 1, etc.
  
  virtual void get_child_ref_entities(DLIList<RefEntity*>& entity_list); 
    //- Appends all immediate (child) RefEntities owned by this RefEntity to 
    //- entity_list. (The query goes down just one dimension.)
  
  void get_all_child_ref_entities(DLIList<RefEntity*>& entity_list);
    //- Appends all child RefEntities owned by this entity to entity_list.
    //- (The query recurses all the way down to RefEntities of dimension 0).
  
  static void get_all_child_ref_entities(DLIList<RefEntity*>& input_list,
				  DLIList<RefEntity*>& output_list );
  //- Appends all child RefEntities owned by entities in input_list to output_list


  static void gather_bdy_entities( DLIList<RefEntity*> &entity_list, 
                                   DLIList<RefEntity*> &bdy_list );
    // Gather the boundary entities of the entity_list into the bdy_list. 
    // Entities appear once in bdy_list, and will not appear in bdy_list
    // if they are already in the entity_list.
    // Uses listMark.

  virtual void get_parent_ref_entities(DLIList<RefEntity*>& entity_list);
    //- Appends all RefEntities that own this (parent RefEntities) to 
    //- entity_list.
    //- (The query goes up just one dimension. For example, if this is a
    //- vertex, the resulting list contains only RefEdges).
  
  void get_all_parent_ref_entities(DLIList<RefEntity*>& entity_list,
                                   const int get_bodies = CUBIT_FALSE );
    //- Appends all parent RefEntities owned by this entity to entity_list.
    //- Recurses up to RefVolumes, or RefBodies if get_bodies is true.
  
  static void change_to_parent_ref_entities( DLIList<RefEntity*>& ancestors );
    //- Modify the input list to contain the list of RefEntities that are
    //- the parents of each of the RefEntities in the original list.
  
  RefEntity *join( RefEntity* ref_entity_2, DLIList<RefEntity*> &join_set );
  static RefEntity *join( DLIList<RefEntity*> &ref_entities, 
                          DLIList<RefEntity*> &join_set );
    // RefEntity* join( RefEntity* ref_entity_2 );
    //- Computes the geometric "join" of elements (elements on the list or 
    //- this and ref_entity_2).  
    //- Definition "Join" = The lowest dimensional entitity that
    //- is a higher dimensional ancestor of them all. 
    //- Note join could be one of the entities itself, NULL, 
    //- or multiple elements.
    //- E.g. The join of a vertex and a curve containing the vertex is 
    //- the curve. The join of two curves of a split cylinder is 
    //- both containing surfaces.
    //- The join of two entities in separate, unmerged volumes is null.
    //- Returns the first element of the join_set, or NULL if set is empty.
  
  RefEntity *meet( RefEntity* ref_entity_2, DLIList<RefEntity*> &join_set );
  static RefEntity *meet( DLIList<RefEntity*> &ref_entities, 
                          DLIList<RefEntity*> &join_set );
    //- like join, except returns the lower order entities common to the input 
    //- entities

  int valence(RefEntity *parent = NULL);
    //- the valence of this entity with respect to parent (absolute
    //- valence if parent is null)

  virtual CubitBoolean is_child(RefEntity *entity);
  virtual CubitBoolean is_parent(RefEntity *entity);
    //- Return TRUE if this is entity, or a direct child (parent) of entity.

  int num_parent_ref_entities();
    //- returns the number of parent entities of this; also useful for determining
    //- whether an entity is free or not; returns -1 on error
  
//   void common_draw_label(int color, const int label_style);
//     //- Common code used for drawing labels. Actually, the entire
//     //- draw_label function could be implemented here if the
//     //- is_*_labeling_on() functions were virtual or better implemented.

  //CubitBoolean is_free_ref_entity();
    //- return CUBIT_TRUE if this ref entity has no non-virtual parents

  virtual CubitVector center_point();
    //- Return the approximate (spatial) center of this RefEntity
  
  virtual double measure();
    //- A generic geometric extent function.
    //- Returns volume for RefVolumes, area for RefFaces, length for RefEdge,
    //- and 1.0 for RefVertices
    //- A RefGroup calculates the maximum dimension of its contained
    //- entities and returns the sum of the measures() of all entities
    //- of that dimension.
    //- Default return value is 0.0 for all other entities.
  
  virtual CubitString measure_label();
    //- Returns the type of measure: (volume, area, length, or N/A)
  
  virtual int validate();
    //- Perform checks to see if entity valid.
  
  virtual DagType dag_type() const = 0;
  virtual const type_info& entity_type_info() const = 0;
  static DagType dag_type( const type_info& );
    //- Return the type of the entity
  
  DagType get_parent_ref_entity_type() const;
  DagType get_child_ref_entity_type() const;
  static DagType get_parent_ref_entity_type( DagType child_type );
  static DagType get_child_ref_entity_type( DagType parent_type );
  
  // send event to all observers (static and non-static) for 
  // this entity and all children
  void notify_sub_all_observers(const CubitEvent& event);

  
  void notify(RefEntity* partner, CubitEventType event);
    //R void
    //I partner
    //I- The merge partner for this object
    //I event
    //I- The type of event
    //- This function takes actions depending on the type of event it
    //- is notified of.
    //-   COMPARISON_FOUND:
    //-     Make temporary TDCompare objects and attach to "this" and
    //-     the "partner" object.
  
  void add_compare_data(RefEntity* partner) ;
    //R void
    //I partner
    //I- The compare partner for this object
    //- This function makes the connection between the two RefEntities,
    //- this and partner. At the end of this function the two entities
    //- would know who they compare with.
  
  void remove_compare_data() ;
    //R void
    //- This function clears the compare related temporary data.
    
  RefEntity* get_compare_partner() ;
    //R RefEntity*
    //R- The partner set in add_compare_data(), or NULL if none
    //R- has been set.

//========  Change Code by RY of Cat,  5/7/99 11:08:12 AM  ========
  void get_related_entity_list(const type_info& related_entity_type,
			       DLIList<RefEntity*>& entity_list);
  //- to parse group in <ref_entity> commands
//========  Change End by RY of Cat,  5/7/99 11:08:12 AM  ========

  virtual void set_id(int i );
    //- set the id of this entity to i
  
  void set_id(int i, CubitBoolean emit_event );

  static const type_info& get_entity_type_info(const char* entity_type);
    //- Returns the type of a class given the class name.
    
  static DagType dag_type( const char* cli_type_name );

  virtual void color(int value);
  virtual int color() const;

protected :

  int          autoMergeStatus;//- Whether entity will participate
                               //- in 'merge all XXX' operation.
  int markedFlag;      //- Scratch flag for algorithm use.
                       //- NOTE: should be Bit markedFlag: 8 ntfolwe
  Bit listFlag : 1;    //- Scratch flag for low-level use

  int mColor; // color of this entity
   
private:

  RefEntity( const RefEntity& );
  void operator=( const RefEntity& );

  void list_mark(int value);
  int  list_mark();
    //- Generic flag for marking whether an entity is in a list or not.
    //- For internal use by RefEntity.

};

// ********** BEGIN INLINE FUNCTIONS       **********

inline void
RefEntity::marked(int value)
{markedFlag = value;}

inline int
RefEntity::marked()
{return static_cast<int>(markedFlag);}

inline void
RefEntity::list_mark(int value)
{listFlag = value;}

inline int
RefEntity::list_mark()
{return listFlag;}

// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN HELPER CLASS DECLARATIONS **********
// ********** END   HELPER CLASS DECLARATIONS **********

#endif

