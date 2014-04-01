//- Class:       RefGroup
//- Description: The RefGroup class contains the reference entity structure
//-              for bodies. It contains the standard information such as
//-              meshed status, etc.
//- Owner:       Ted Blacker
//- Checked by:  Greg Sjaardema, 2/28/94
//- Version: $Id: 

#ifndef REFGROUP_HPP
#define REFGROUP_HPP

#include "RefCollection.hpp"
#include "CubitObserver.hpp"

class CubitBox;
template <class X> class DLIList;
class RefGroup;
class RefVolume;
class RefFace;
class RefEdge;
class RefVertex;
class RefEntity; 

class CUBIT_GEOM_EXPORT RefGroup : public RefCollection,
                 public CubitObserver
{ 	

public:

  friend class RefEntityFactory;
    //- the factory is allowed to call the (private) constructors

  friend class TDCAGE;
  
protected:  /* Delete using RefGroup::delete_group(..) */
  virtual ~RefGroup();
    //- Class destructor

public:

  static const char* get_class_name()
     {
       return "Group";
     }

  virtual const char* class_name() const
     {
       return get_class_name();
     }
  
  int maximum_dimension();
    //- This routine returns the maximum dimension of its owned subentities.
    //- The 'only' kludge is that if one of its subentities is a RefGroup, it
    //- must call 'maximum_dimension' on that entity instead of dimension().

  virtual const type_info& entity_type_info() const
      { return typeid(RefGroup); }
    //- return the type for this geometryEntity
    
  DagType dag_type() const
    { return DagType::invalid_type(); }

  virtual CubitStatus add_ref_entity(RefEntity *ref_entity);
  virtual CubitStatus add_ref_entity(DLIList<RefEntity*>& entity_list);
  virtual CubitStatus remove_ref_entity(RefEntity *ref_entity,
                                        const CubitBoolean from_observable = CUBIT_FALSE);
  virtual CubitStatus remove_ref_entity(DLIList<RefEntity*>& entity_list,
                                        const CubitBoolean from_observable = CUBIT_FALSE);
  
    //- add or remove one or more ref entities to/from this group

  virtual void get_sub_entities(DLIList<RefEntity*> &entity_list);
    //- appends all ref entities owned by this entity and recurses
    //- down to dimension 0

  virtual void is_mergeable(AutoMergeStatus val);
  bool is_mergeable() { return RefEntity::is_mergeable(); }
  
  int remove_all_ref_entities();
  //- remove all geometry entities in the group.  
  //- Returns number of items removed.

  virtual void get_parent_ref_entities(DLIList<RefEntity*>&);
    //- appends all ref entities that own this to entity_list.
    //- Goes up just one dimension.

  virtual void get_child_ref_entities(DLIList<RefEntity*>& entity_list);
    //- appends all immediate ref entities owned by this entity on entity_list
    //- Goes down just one dimension.
  
  virtual void get_child_entities(DLIList<CubitEntity*>& cub_entity_list);
    //- appends all immediate ref entities owned by this entity on entity_list
    //- Goes down just one dimension.

  void expand_group( DLIList<RefEntity*>& entity_list ); 
    //- appends all the ref entities owned by this group.  It will go down
    //- until there are no ref-groups in the entity_list.

  virtual CubitBox bounding_box();
  virtual CubitVector center_point();
  
  int subtract(RefGroup *group_to_subtract, RefGroup *target_group);
    //- subtract group_to_subtract from this group
  
  int intersect(RefGroup *other_group, RefGroup *target_group);
    //- intersect other_group with this one
  
  int unite(RefGroup *other_group, RefGroup *target_group);
    //- unite other_group with this one

  virtual int validate();
    //- Return number of problems detected, 0 if none.
  
  //virtual void draw (int color);
    //- draw the group's contained entities

  static CubitStatus delete_group (RefGroup *group_ptr,
                                   CubitBoolean propagate = CUBIT_FALSE);
    //- deletes a specified group from the global group list
    //- if the boolean "propagate" is true, the groups sub_groups are 
    //- deleted also; as well as the sub_group's groups, etc..

  static void delete_all_groups ();
    //- deletes all the groups from the model (except for 'picked_group' and 'drawn_group')

  static void get_contained_groups (RefGroup *group_ptr,
                                    DLIList<RefGroup*> &contained_groups);
    //- gets the groups owned by group_ptr, as well as any other groups
    //- owned by these groups, etc..  Current group (group_ptr) is also 
    //- added to the list.

  static void get_groups_within( CubitEntity* cubit_entity_ptr, 
                                 DLIList<RefGroup*> &groups_within,
                                 const CubitBoolean recursive = CUBIT_TRUE);
  static void get_groups_within( RefEntity* ref_entity_ptr, 
                                 DLIList<RefGroup*> &groups_within,
                                 const CubitBoolean recursive = CUBIT_TRUE);
  static void get_groups_within( RefGroup* ref_group_ptr, 
                                 DLIList<RefGroup*> &groups_within,
                                 const CubitBoolean recursive = CUBIT_TRUE);
    //- Finds those groups which contain the input entity, at any level
    //-  (i.e., if group 2 contains group 3 which contains the entity, then
    //-  group 2 and group 3 will be returned in the output list).

  CubitStatus notify_observer(CubitObservable *observable,
                              const CubitEvent &observer_event,
                              CubitBoolean from_observable = CUBIT_FALSE);
    //- handle notify observer function

protected:

  DLIList<RefEntity*>  entityList;
  
  int recursionMark;
  
  RefGroup(int proe_type); 
    // For Pro/E parts and assemblies (avoids notify). 1=assembly, 2=part
    //- Class contructors

  RefGroup(const char* name = NULL, int id = 0);
  RefGroup (DLIList<RefEntity*>& entity_list);
    //- need constructors accessable from derived classes
  
private:

};

#endif // REFGROUP_HPP

