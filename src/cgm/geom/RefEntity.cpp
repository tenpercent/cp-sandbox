// ********** BEGIN STANDARD INCLUDES      **********

#include <assert.h>

// ********** END STANDARD INCLUDES        **********

// ********** BEGIN MOTIF INCLUDES         **********
// ********** END MOTIF INCLUDES           **********

// ********** BEGIN OPEN INVENTOR INCLUDES **********
// ********** END OPEN INVENTOR INCLUDES   **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "RefEntity.hpp"
#include "CubitDefines.h"
#include "CubitMessage.hpp"
#include "CubitString.hpp"

#include "Body.hpp"
#include "BodySM.hpp"
#include "RefGroup.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefEntityName.hpp"
#include "RefCollection.hpp"
#include "RefVertex.hpp"

#include "CollectionEntity.hpp"
#include "IdSetEvent.hpp"

#include "ModelQueryEngine.hpp"

#include "TDCompare.hpp"


#include "GeometryQueryTool.hpp"
#include "RefEntityFactory.hpp"

#include "DLIList.hpp"
#include "CubitUtil.hpp"

// need to include VGE here for free entity test; since virtual geom
// is so tied into ref entities, I think it's ok (tjt)
//#include "VirtualGeometryEngine.hpp"
// Dependency on virtual classes removed, so the VGE does not need
// to be included. 10-18-01 (wjg)

// The following include generates the inline functions for casting
// a list of A objects to a list of B objects.
#include "CastTo.hpp"

#include "MergeTool.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

RefEntity::RefEntity()
{
  autoMergeStatus      = AUTO_MERGE_AUTO | AUTO_MERGE_ON;
  markedFlag           = CUBIT_FALSE;
  listFlag             = CUBIT_FALSE;
  mColor = CUBIT_DEFAULT_COLOR;
}

RefEntity::~RefEntity()
{
  // Remove the name of this entity from the entity name map
  RefEntityName::instance()->remove_refentity_name(this, CUBIT_FALSE);
  
  remove_from_observers();
}

RefEntity* RefEntity::get_by_name(const CubitString& name)
{
  return RefEntityName::instance()->get_refentity(name);
}

CubitStatus RefEntity::entity_name (CubitString name)
{
    // Add the new {RefEntity, Name} tuple to the RefEntityNameMap class.
    // Note that if this name already exists for another RefEntity, then
    // it will return with an error. 
  CubitStatus success = RefEntityName::instance()->
    add_refentity_name(this, name);
  
  if (success == CUBIT_SUCCESS)
    CubitObserver::notify_static_observers(this, ENTITY_NAME_CHANGED);

  return success;
}

CubitString RefEntity::entity_name() const
{
  DLIList<CubitString*> names;
  
    // Get the name(s) associated with this RefEntity.
  RefEntityName::instance()->get_refentity_name(this, names, CUBIT_TRUE);
  
  CubitString name;
  
    // If there are no names, then create the default name and return it.
  if (names.size() == 0) 
  {
    name = class_name();
    name += " ";
    name += CubitString(entityId);
  } 
  else 
  {
    name = *(names.get());
  }
  
  return name;
}

int RefEntity::num_names() const
{
  DLIList<CubitString*> names;
  
  // Get the name(s) associated with this RefEntity.
  RefEntityName::instance()->get_refentity_name(this, names, CUBIT_TRUE);
  return names.size();
}
  
void RefEntity::entity_names(DLIList<CubitString*>& names) const
{
  RefEntityName::instance()->get_refentity_name(this, names);
}

CubitStatus RefEntity::remove_entity_name(CubitString const & name)
{
  RefEntity* entity = RefEntityName::instance()->get_refentity(name);
  if (entity != this)
    return CUBIT_FAILURE;
  
  RefEntityName::instance()->remove_refentity_name(name);
  CubitObserver::notify_static_observers(this, ENTITY_NAME_CHANGED);
  return CUBIT_SUCCESS;
}

CubitStatus RefEntity::remove_entity_names()
{
  RefEntityName::instance()->remove_refentity_name(this, CUBIT_TRUE);
  CubitObserver::notify_static_observers(this, ENTITY_NAME_CHANGED);
  return CUBIT_SUCCESS;
}

// All this function basically does is assign the default name.
// It will do this as an attribute so they get actuated to the
// solid model.
CubitStatus RefEntity::generate_default_name( CubitString &name )
{
  if ( CAST_TO( this, Body ) )
     name += CubitString("bod");
  else if ( CAST_TO( this, RefVolume ) )
     name += CubitString("vol");
  else if ( CAST_TO( this, RefFace ) )
     name += CubitString("sur");
  else if ( CAST_TO( this, RefEdge ) )
     name += CubitString("cur");
  else if ( CAST_TO( this, RefVertex ) )
     name += CubitString("ver");
  else if ( CAST_TO( this, RefGroup ) ) {
     //PRINT_INFO( "Debug: RefGroup entity name set in RefEntity::generate_default_name\n" );
    name += CubitString("gro");
  }
  else
  {
    PRINT_ERROR("Invalid entity: %s.\n",
                entity_name().c_str() );
    return CUBIT_FAILURE;
  }
    
  name += CubitString(entityId);
  return CUBIT_SUCCESS;
}

CubitStatus RefEntity::assign_default_name( CubitBoolean user_setting )
{

  if ( RefEntityName::instance()->get_generate_default_names() 
    || user_setting == CUBIT_TRUE )
  {
      // first generate the default name
    CubitString name;
    CubitStatus result = generate_default_name( name );
    if (result == CUBIT_FAILURE) return result;

      // now assign it to this entity and return
    RefEntityName::instance()->add_refentity_name(this, name);
    return CUBIT_SUCCESS;
  }
  return CUBIT_SUCCESS;
}


void RefEntity::merge_entity_names(RefEntity* dead_entity)
{
  RefEntityName::instance()->merge_refentity_names(this, dead_entity);
}

void RefEntity::switch_entity_names(RefEntity* other_entity)
{
  RefEntityName::instance()->switch_refentity_names(this, other_entity);
}


//-------------------------------------------------------------------------
// Purpose       : Target types for parent/child queries
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/24/03
//-------------------------------------------------------------------------
DagType RefEntity::get_parent_ref_entity_type() const
  { return get_parent_ref_entity_type(dag_type()); }
DagType RefEntity::get_parent_ref_entity_type(DagType my_type)
{
  int dimension = my_type.dimension();

    // special case for bodies
  if (dimension == 3)
  {
    if (my_type.functional_type() == DagType::BasicTopologyEntity_TYPE)
      return DagType::body_type();
    else
      return DagType::invalid_type();
  }
  
  return DagType( dimension + 1, DagType::BasicTopologyEntity_TYPE);
}
DagType RefEntity::get_child_ref_entity_type() const
  { return get_child_ref_entity_type(dag_type()); }
DagType RefEntity::get_child_ref_entity_type( DagType my_type )
{
  int dimension = my_type.dimension();
  if (my_type.functional_type() == DagType::BasicTopologyEntity_TYPE)
    return DagType(dimension - 1, DagType::BasicTopologyEntity_TYPE);
  else
    return DagType(dimension, DagType::BasicTopologyEntity_TYPE);
}

DagType RefEntity::dag_type( const type_info& type )
{
  if (type == typeid(Body))
    return DagType::body_type();
  else if (type == typeid(RefVolume))
    return DagType::ref_volume_type();
  else if (type == typeid(RefFace))
    return DagType::ref_face_type();
  else if (type == typeid(RefEdge))
    return DagType::ref_edge_type();
  else if (type == typeid(RefVertex))
    return DagType::ref_vertex_type();
  
  assert(0);
  return DagType::invalid_type();
}

void RefEntity::get_child_ref_entities(DLIList<RefEntity*>& entity_list)
{
  // First get the type of RefEntity that is a child of "this" one
  DagType child_type = get_child_ref_entity_type();
  
  DLIList<ModelEntity*> tempList ;
  
    // Now retrieve the appropriate type of child entities of this one
    // if the child_type is a valid type
  if (child_type.is_valid())
  {
    ModelEntity* modelEntityPtr = CAST_TO(this, ModelEntity) ;
    
    assert(modelEntityPtr != 0);
    
    CubitStatus result = ModelQueryEngine::instance()->
        query_model( *modelEntityPtr, child_type, tempList );
    if (result == CUBIT_FAILURE)
    {
      PRINT_ERROR("In RefEntity::get_child_ref_entities\n");
      PRINT_ERROR("       Query failed for unknown reason.\n");
      return;
    }
   
    CAST_LIST(tempList, entity_list, RefEntity) ;
  }
}


//-------------------------------------------------------------------------
// Purpose       : Get all types of child RefEntity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
void RefEntity::get_all_child_ref_entities(DLIList<RefEntity*>& entity_list)
{
  
  ModelQueryEngine *const mqe = ModelQueryEngine::instance();
  DagType child_type = get_child_ref_entity_type();
  DLIList<ModelEntity*> query_output;
  ModelEntity* this_me = dynamic_cast<ModelEntity*>(this);
  assert(!!this_me);
 
    //While there are more child types
  while (child_type.is_valid())
  {
    mqe->query_model(*this_me, child_type, query_output);
    query_output.reset();
    for (int i = query_output.size(); i--; )
    {
      RefEntity* ref_ptr = dynamic_cast<RefEntity*>(query_output.get_and_step());
      assert(!!ref_ptr);
      entity_list.append(ref_ptr);
    }
    
    child_type = get_child_ref_entity_type( child_type );
  }
}



//-------------------------------------------------------------------------
// Purpose       : Get all parent RefEntities for each passed entity.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
void RefEntity::get_all_child_ref_entities( DLIList<RefEntity*>& input_list,
					                                  DLIList<RefEntity*>& output_list )
{
  DLIList<ModelEntity*> query_output;
  ModelQueryEngine *const mqe = ModelQueryEngine::instance();
  int i;
  
    // Find parent-most type from all passed entities
  DagType target_type = DagType::ref_vertex_type();
  for (i = input_list.size(); i--; )
  {
    DagType type = input_list.step_and_get()->dag_type();
    if (type.is_valid() && type > target_type)
      target_type = type;
  }
  
  target_type = get_child_ref_entity_type( target_type );
  while (target_type.is_valid())
  {
#ifdef BOYD17 
    ModelQueryEngine::BeginQuery lock;
#endif
    
    input_list.reset();
    for (i = input_list.size(); i--; )
    {
      RefEntity* input_ptr = input_list.get_and_step();
      DagType source_type = input_ptr->dag_type();

      if (source_type.is_valid() && source_type > target_type)
      {
        ModelEntity* me_ptr = dynamic_cast<ModelEntity*>(input_ptr);
        assert(!!me_ptr);
        mqe->query_model (*me_ptr, target_type, query_output);
        
        query_output.size();
        for ( int j = query_output.size(); j--; )
        {
          RefEntity* ref_ent = dynamic_cast<RefEntity*>(query_output.get_and_step());
          output_list.append( ref_ent );
        }
      }
    }
    
    target_type = get_child_ref_entity_type( target_type );

  }
}


void RefEntity::get_parent_ref_entities(DLIList<RefEntity*>& entity_list)
{

  // First get the type of RefEntity that is a child of "this" one
  DagType parent_type = get_parent_ref_entity_type();;
  
  DLIList<ModelEntity*> tempList ;
  
  // Now retrieve the appropriate type of child entities of this one
  // if the child_type is a valid type
  if (parent_type.is_valid())
  {
    ModelEntity* modelEntityPtr = CAST_TO(this, ModelEntity) ;
    
    assert(modelEntityPtr != 0);
    
    CubitStatus result = ModelQueryEngine::instance()->
        query_model( *modelEntityPtr, parent_type, tempList );
    if (result == CUBIT_FAILURE)
    {
      PRINT_ERROR("In RefEntity::get_parent_ref_entities\n");
      PRINT_ERROR("       Query failed for unknown reason.\n");
      return;
    }
   
    CAST_LIST(tempList, entity_list, RefEntity) ;
  }
}

//-------------------------------------------------------------------------
// Purpose       : Get all types of parent RefEntity
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
void RefEntity::get_all_parent_ref_entities(DLIList<RefEntity*>& entity_list,
                                            const int get_bodies)
{
  
  ModelQueryEngine *const mqe = ModelQueryEngine::instance();
  DagType parent_type = get_parent_ref_entity_type();
  DLIList<ModelEntity*> query_output;
  ModelEntity* this_me = dynamic_cast<ModelEntity*>(this);
  assert(!!this_me);
 
    //While there are more parent types
  while ( parent_type.is_valid() && 
          (get_bodies || (parent_type != DagType::body_type())) )
  {
    mqe->query_model(*this_me, parent_type, query_output);
    query_output.reset();
    for (int i = query_output.size(); i--; )
    {
      RefEntity* ref_ptr = dynamic_cast<RefEntity*>(query_output.get_and_step());
      assert(!!ref_ptr);
      entity_list.append(ref_ptr);
    }
    
    parent_type = get_parent_ref_entity_type( parent_type );
  }
}


void RefEntity::change_to_parent_ref_entities( DLIList<RefEntity*>& ancestors )
{
  DLIList<RefEntity*> new_ancestors, temp_list;
  if ( ancestors.size() == 0 ) return;
  
  ancestors.reset();
  do 
  {
    temp_list.clean_out();
    ancestors.get_and_step()->get_parent_ref_entities( temp_list );
    new_ancestors.merge_unique( temp_list );
  } while ( !ancestors.is_at_beginning() );
  
    // Overwrite the input list
  ancestors = new_ancestors;
}

int RefEntity::dimension() const
{
    // Virtual function. This is the default return value. 
    // Override in subclass if different return value is needed.
  return -1;
}

 
RefEntity *RefEntity::join( DLIList<RefEntity*> &ref_entities, 
                            DLIList<RefEntity*> &join_set )
{
  join_set.clean_out();

    // Trivial cases 
  
    // empty set
  if ( !ref_entities.size() )
    return NULL;
  
    // All entities are the same, e.g. list is size 1
  int i, j;
  int all_same = CUBIT_TRUE;
  RefEntity *e1, *e2;
  for (i = ref_entities.size()-1; all_same && i--; )
  {
    e1 = ref_entities.get_and_step();
    e2 = ref_entities.get();
    all_same = ( e1 == e2 );
  }
  if ( all_same )
  {
    join_set.append( ref_entities.get() );
    return join_set.get();
  }  
  
    // Generic test, build ancestor lists until dimensions are equal
  const int size = ref_entities.size();

    // Lists of ancestors
  DLIList<RefEntity*> *ancestors = new DLIList<RefEntity*> [ size ];
    // current dimension of ancestors of each entity
  int *dimension = new int [size];
    // Set ancestors and dimension: Start with the entities themselves
  int d1, max_dimension = -1;
  ref_entities.reset();
  for ( i = 0; i < size; i++ )
  {
    e1 = ref_entities.get_and_step();
    ancestors[i].append( e1 );
    dimension[i] = d1 = e1->dimension();
    if ( d1 > max_dimension )
      max_dimension = d1;
  }
  
  do
  {

    // Bump dimensions up to max dimension
    for ( i = 0; i < size; i++ )
    {
      while ( dimension[i] < max_dimension ) 
      {
        change_to_parent_ref_entities( ancestors[i] );
          // no ancestor of max dimension -> join is nothing
        if ( !ancestors[i].size() )
          goto clean_return;
        dimension[i] = d1 = ancestors[i].get()->dimension();
          // if somehow a dimension was skipped, go through the list again.
        if ( d1 > max_dimension )
        {
          max_dimension = d1;
          i = 0; // restart loop
        }
      }
    }

      // dimensions are the same, see if any entities common to *all* lists
    for ( j = ancestors[0].size(); j--; )
    {
      e1 = ancestors[0].get_and_step();
      int in_all_lists = CUBIT_TRUE;
      for ( i = 1; in_all_lists && i < size; i++ )
      {
        in_all_lists = ancestors[i].is_in_list( e1 );
      }
      if ( in_all_lists )
        join_set.append( e1 );
    } 
    
      // iterate
    max_dimension++;
    
      // they had common entities in the current dimension - return
  } while ( join_set.size() == 0 );
  
    clean_return:
  delete [] dimension;
  delete [] ancestors;

  if ( join_set.size() )
  {
    join_set.reset();
    return join_set.get();
  }
  return NULL;
}

  
RefEntity *RefEntity::join( RefEntity* ref_entity_2, DLIList<RefEntity*> &join_set )
{
  join_set.clean_out();
  assert( this != 0 );
  
    // Trivial case :-)
    // Both entities are the same, or 
    // The join of the empty set and any thing is defined to be the thing.
  if (this == ref_entity_2 || ref_entity_2 == NULL ) 
  { 
    join_set.append( this );
    return this;
  }
  
    // call generic function that takes a list
  DLIList<RefEntity*> both_entities;
  both_entities.append( this );
  both_entities.append( ref_entity_2 );
  return join( both_entities, join_set );
}

RefEntity *RefEntity::meet( DLIList<RefEntity*> &ref_entities, 
                                  DLIList<RefEntity*> &join_set )
{
  join_set.clean_out();

    // Trivial cases 
  
    // empty set
  if ( !ref_entities.size() )
    return NULL;
  
    // All entities are the same, e.g. list is size 1
  int i, j;
  int all_same = CUBIT_TRUE;
  RefEntity *e1, *e2;
  for (i = ref_entities.size()-1; all_same && i--; )
  {
    e1 = ref_entities.get_and_step();
    e2 = ref_entities.get();
    all_same = ( e1 == e2 );
  }
  if ( all_same )
  {
    ref_entities.get()->get_child_ref_entities(join_set);
    return join_set.get();
  }  

    // they aren't all the same; get all the children, intersect the lists,
    // and remove any related entities
  DLIList<RefEntity*> temp_list;
  ref_entities.get()->get_all_child_ref_entities(join_set);
  join_set.append(ref_entities.get_and_step());
  
  for (i = ref_entities.size(); i > 1; i--) {
    temp_list.clean_out();
    ref_entities.get()->get_all_child_ref_entities(temp_list);
    temp_list.append(ref_entities.get_and_step());
    join_set.intersect(temp_list);

      // if the join set is null at any time, there's no common entity
    if (join_set.size() == 0) return NULL;
  }

  temp_list.clean_out();

    // clean out any duplicates before checking for relations
  temp_list.merge_unique(join_set);
  join_set = temp_list;
  temp_list.clean_out();
  
    // now remove any related entities
  join_set.reset();
  join_set.step();
  for (i = join_set.size()-1; i > 0; i--) {
    RefEntity *entity = join_set.get();
    if (temp_list.move_to(entity)) continue;
    
    for (j = join_set.size()-i; j > 0; j--) {
      RefEntity *other_entity = join_set.prev(j);
      if (temp_list.move_to(other_entity)) continue;
      if (entity->is_child(other_entity)) temp_list.append(entity);
      else if (entity->is_parent(other_entity)) temp_list.append(other_entity);
    }
    join_set.step();
  }
  
  if (temp_list.size() > 0) join_set -= temp_list;

  return (join_set.size() ? join_set.get() : NULL);
}

RefEntity *RefEntity::meet( RefEntity* ref_entity_2, DLIList<RefEntity*> &join_set )
{
  join_set.clean_out();
  assert( this != 0 );
  
    // Trivial case :-)
    // Both entities are the same, or 
    // The join of the empty set and any thing is defined to be the thing.
  if (this == ref_entity_2 || ref_entity_2 == NULL ) 
  { 
    join_set.append( this );
    return this;
  }
  
    // call generic function that takes a list
  DLIList<RefEntity*> both_entities;
  both_entities.append( this );
  both_entities.append( ref_entity_2 );
  return meet( both_entities, join_set );
}

int RefEntity::valence(RefEntity *parent) 
{
  DLIList<RefEntity*> parents;
  get_parent_ref_entities(parents);
  if (parent == NULL) return parents.size();
  
  int val = parents.size();
  int i;
  BasicTopologyEntity *topo_ent = CAST_TO(this, BasicTopologyEntity);
  for (i = val; i > 0; i--) {
    RefEntity *other_ent = parents.get_and_step();
    if (!topo_ent->is_directly_related(CAST_TO(other_ent, BasicTopologyEntity))) val--; 
  }
  
  return val;
}
  
CubitVector RefEntity::center_point()
{
  return bounding_box().center();
}

// autoMergeStatus :
// bit 0:  1 = mergeable, 0 = not mergeable
// bit 1:  1 = auto,  0 = explicit

void RefEntity::is_mergeable(AutoMergeStatus val)
{
  AutoMergeStatus old_status = (AutoMergeStatus)(autoMergeStatus & 1);
  autoMergeStatus = (int)val;
  
    // always want to update children, but update_auto_merge_state
    // recursively updates the children, so don't update the
    // children twice if we call update_auto_merge_state.
  if( val == AUTO_MERGE_AUTO )
  {
    this->update_auto_merge_state();
  }
  else if( old_status != val )
  {
    DLIList<RefEntity*> children;
    get_child_ref_entities( children );
    for( int i = children.size(); i--; )
      children.get_and_step()->update_auto_merge_state();
  }
}

bool RefEntity::is_mergeable()
{
    // return 0-bit of autoMergeState
  return (bool)(autoMergeStatus & 1);
}

CubitBoolean RefEntity::is_merged()
{
  TopologyEntity* topo = CAST_TO(this, TopologyEntity);
  if (!topo)
    return CUBIT_FALSE;

  return MergeTool::instance()->entity_merged(topo);
}

AutoMergeStatus RefEntity::merge_status() const
{
  return autoMergeStatus > 1 ? AUTO_MERGE_AUTO : (AutoMergeStatus)autoMergeStatus;
}

bool RefEntity::children_mergeable()
{
   DLIList<RefEntity*> children;
   get_all_child_ref_entities( children );
   
   for( int i = children.size(); i--; )
     if( children.get_and_step()->autoMergeStatus == 0 ) 
       return false;

   return true;
}

void RefEntity::update_auto_merge_state()
{
  int i;
  
  if( autoMergeStatus & 2 )
  {
    int mergeable = 1;
    
    DLIList<RefEntity*> list;

    get_parent_ref_entities( list );
    for( i = list.size(); i--; )
      if( ! list.get_and_step()->is_mergeable() )
        mergeable = 0;
  
    int old_value = autoMergeStatus & 1;
    if( old_value != mergeable )
    {
      autoMergeStatus = (autoMergeStatus & 2) | mergeable;
  
      list.clean_out();
      get_child_ref_entities( list );
    
      for( i = list.size(); i--; )
        list.get_and_step()->update_auto_merge_state();
    }    
  }
}

int RefEntity::can_modify() { return 1; }

double RefEntity::measure()
{
  return 0.0;
}

CubitString RefEntity::measure_label()
{
  return "N/A";
}


void RefEntity::notify_sub_all_observers(const CubitEvent& event)
{
  DLIList<RefEntity*> entity_list;
  get_all_child_ref_entities( entity_list );
  for ( int i = entity_list.size(); i>0; i-- )
    entity_list.get_and_step()->notify_all_observers(event);
  this->notify_all_observers(event);
}


//-------------------------------------------------------------------------
// Purpose       : This function takes actions depending on the type of 
//                 event it is notified of.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 11/25/96
//-------------------------------------------------------------------------

void RefEntity::notify(RefEntity* partner, CubitEventType event)
{
  if ( event == COMPARISON_FOUND )
  {
    add_compare_data(partner) ;
    MergeTool::instance()->compare_notify(this, COMPARISON_FOUND) ;
    MergeTool::instance()->compare_notify(partner, COMPARISON_FOUND) ;
  }
}

//-------------------------------------------------------------------------
// Purpose       : This function makes the connection between the two 
//                 RefEntities, this and partner. At the end of this 
//                 function the two entities would know who they 
//                 compare with.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 11/25/96
//-------------------------------------------------------------------------

void RefEntity::add_compare_data(RefEntity* partner) 
{
  TDCompare* compareDataPtr = (TDCompare*)(this->get_TD(&TDCompare::is_compare));
  if (compareDataPtr == NULL)
  {
    compareDataPtr = new TDCompare() ;
    this->add_TD(compareDataPtr) ;
  }
  compareDataPtr->set_compare_partner(partner) ;
  
  
  compareDataPtr = (TDCompare*)(partner->get_TD(&TDCompare::is_compare));
  if (compareDataPtr == NULL)
  {
    compareDataPtr = new TDCompare() ;
    partner->add_TD(compareDataPtr) ;
  }
  compareDataPtr->set_compare_partner(this) ;
}

//-------------------------------------------------------------------------
// Purpose       : This function clears the compare related temporary data. 
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 11/25/96
//-------------------------------------------------------------------------

void RefEntity::remove_compare_data()
{
  ToolData* tdPtr = this->get_TD(&TDCompare::is_compare) ;
  TDCompare* tdComparePtr = CAST_TO(tdPtr, TDCompare) ;
  
  if (tdComparePtr == NULL)
  {
    return ;
  }
  RefEntity* partner = tdComparePtr->get_compare_partner() ;

  if (partner != NULL)
    {
      PRINT_DEBUG_19(
		  "RefEntity::remove_compare_data - Removing TDCompare from"
		  " %s %d, via partner %s %d\n",
		  partner->class_name(),
		  partner->id(),
		  this->class_name(),
		  this->id());
      partner->delete_TD(&TDCompare::is_compare) ;
    }
  
  this->delete_TD(&TDCompare::is_compare) ;
//   partner->delete_TD(&TDCompare::is_compare) ;
}

RefEntity* RefEntity::get_compare_partner()
{
  ToolData* tdPtr = get_TD(&TDCompare::is_compare) ;
  TDCompare* tdComparePtr = CAST_TO(tdPtr, TDCompare) ;
  return tdComparePtr ? tdComparePtr->get_compare_partner() : 0;
}

int RefEntity::validate()
{
    //- This function determines whether the entity is valid.
    //- Several types of checks can be done, 
  
    // Check that measure is positive
  int error = 0;
  double this_measure = measure();
  if (this_measure <= 0.0) {
    PRINT_WARNING("\tWARNING: non-positive %s (%f) for %s, (%s %d)\n",
                  measure_label().c_str(), this_measure,
                  entity_name().c_str(), class_name(), id());
    error++;
  }    
  return error;
}


CubitBoolean RefEntity::is_child(RefEntity *entity)
{
    // same?
  if ( this == entity )
    return CUBIT_TRUE;
  
    // wrong dimensions?
  if (!entity || dimension() >= entity->dimension() )
    return CUBIT_FALSE;
  
    // Topological query.
    // Usually slightly faster to search up in dimension.
  TopologyEntity *topo_this = CAST_TO( this, TopologyEntity );
  TopologyEntity *topo_entity = CAST_TO( entity, TopologyEntity );  
  return topo_this->is_directly_related( topo_entity );
}

CubitBoolean RefEntity::is_parent(RefEntity *entity)
{
    // same?
  if ( this == entity )
    return CUBIT_TRUE;
  
    // wrong dimensions?
  if (!entity || dimension() <= entity->dimension() )
    return CUBIT_FALSE;

    // Topological query. 
    // Usually slightly faster to search up in dimension.
  BasicTopologyEntity *topo_this = CAST_TO( this, BasicTopologyEntity );
  BasicTopologyEntity *topo_entity = CAST_TO( entity, BasicTopologyEntity );  
  return topo_entity->is_directly_related( topo_this );
}

int RefEntity::num_parent_ref_entities()
{
    // First get the type of RefEntity that is a parent of "this" one
  DagType parent_type = get_parent_ref_entity_type();
  
  DLIList<ModelEntity*> tempList ;
#ifdef BOYD17 
  DLIList<RefEntity*> tmp_ref_list ;
#endif
  
    // Now retrieve the appropriate type of parent entities of this one,
    // if the parent_type is a valid type
  if (parent_type.is_valid())
  {
    ModelEntity* modelEntityPtr = CAST_TO(this, ModelEntity) ;
    
      //Make sure that we have a valid pointer
    assert(modelEntityPtr != NULL) ;
    
    CubitStatus result = ModelQueryEngine::instance()->
        query_model( *modelEntityPtr,
                     parent_type,
                     tempList );
    if (result == CUBIT_FAILURE)
    {
      PRINT_ERROR("In RefEntity::num_parent_ref_entities\n");
      PRINT_ERROR("       Query failed for unknown reason.\n");
      return -1;
    }
    return tempList.size();
  }

  return -1;
}
/*
CubitBoolean RefEntity::is_free_ref_entity() 
{
    // return true if this RefEntity has no non-virtual parents
    // first look for Body type or # of parents
  if (entity_type_info() == typeid(Body) ||
      num_parent_ref_entities() == 0) return CUBIT_TRUE;

  
    // now look at immediate parents to quickly rule out some entities
  DLIList<RefEntity*> parents;
  get_parent_ref_entities(parents);
  if (!VGE->check_for_virtual(CAST_TO(parents.get(), TopologyEntity)))
    return CUBIT_FALSE;
  
    // now look at all the parents; work back from the list, as this is most
    // likely to be non-virtual
  parents.clean_out();
  get_all_parent_ref_entities(parents);
  parents.last();
  RefEntity *ref_entity;
  int i;
  for (i = parents.size(); i > 0; i--) {
    ref_entity = parents.get_and_back();
    if (!VGE->check_for_virtual(CAST_TO(ref_entity, TopologyEntity)))
      return CUBIT_FALSE;
  }

    // if we've gotten here, we have no non-virtual parents
  return CUBIT_TRUE;
}
*/
void RefEntity::gather_bdy_entities( DLIList<RefEntity*> &entity_list, 
                                     DLIList<RefEntity*> &bdy_list ) 
{
  RefEntity *entity;
  DLIList<RefEntity*> tmp_bdy_list;
  int i;

  for ( i = entity_list.size(); i--; )
  {
    entity = entity_list.get_and_step();
    entity->list_mark( CUBIT_TRUE );
    //entity->get_all_child_ref_entities( tmp_bdy_list );
  }
    
  get_all_child_ref_entities(entity_list, tmp_bdy_list);


    // copy non-duplicate and non-entity_list entities
  for ( i = tmp_bdy_list.size(); i--; ) {
    entity = tmp_bdy_list.get_and_step();
    if ( !entity->list_mark() ) {
      bdy_list.append( entity );
      entity->list_mark( CUBIT_TRUE );
    }
  }
    // clean-up
  for ( i = bdy_list.size(); i--; ) {
    entity = bdy_list.get_and_step();
    entity->list_mark( CUBIT_FALSE );
  }
  for ( i = entity_list.size(); i--; ) {
    entity = entity_list.get_and_step();
    entity->list_mark( CUBIT_FALSE );
  }
}

//-------------------------------------------------------------------------
// Purpose       : Most RefEntites related_to another RefEntity can be done
//                 through query model as is in CommandHandler, but
//                 group in <ref_entity> crashes when query model is called, instead
//                 this function will be called
//
// Special Notes : 
//
// Creator       : RY (CAT)
//
// Creation Date : 5-99
//-------------------------------------------------------------------------
void RefEntity::get_related_entity_list(const type_info& related_entity_type,
			                DLIList<RefEntity*>& entity_list)
{
  if (related_entity_type == typeid(RefGroup)){
    DLIList <CubitObserver*> observer_list;
    this->get_observer_list (observer_list);
    for (int i = observer_list.size(); i > 0; i--){
      entity_list.append_unique (CAST_TO (observer_list.get(), RefEntity));
      observer_list.step();
    }
  }
}

void RefEntity::set_id(int i )
{
  set_id( i, CUBIT_TRUE );
}

void RefEntity::set_id(int i, CubitBoolean emit_event )
{
  if (entityId == i)
    return;
 
  int old_id = entityId;
  entityId = i;

  if( emit_event )
    CubitObserver::notify_static_observers( this, IdSetEvent( old_id, entityId ) );

  int old_max = RefEntityFactory::instance()->maximum_id(this);
  
  if (old_max < entityId)
      // Need to reset the maxId for this entitytype
    RefEntityFactory::instance()->incorporate_id(this);
  else if (old_max == old_id) {
      // We just reset the entity with the max id to something less
      // than that - should decrement the max id to restore the max to what
      // it was before this entity was created
    old_max--;
    RefEntityFactory::instance()->maximum_id(this->entity_type_info(), old_max);
  }
    
}

const type_info& RefEntity::get_entity_type_info(const char* entity_type)
{
  CubitString string(entity_type);
  string.to_lower();
  
  if( string == "group" )
     return typeid(RefGroup);
  else if (string == "body" )
     return typeid(Body);
  else if( string == "volume" )
     return typeid(RefVolume);
  else if( string == "surface" )
     return typeid(RefFace);
  else if( string == "curve" )
     return typeid(RefEdge);
  else if( string == "vertex" )
     return typeid(RefVertex);
  else
     return typeid(InvalidEntity);
}

DagType RefEntity::dag_type(const char* name)
{
  if (CubitUtil::compare(name,"body"))
    return DagType::body_type();
  else if (CubitUtil::compare(name,"volume"))
    return DagType::ref_volume_type();
  else if (CubitUtil::compare(name,"surface"))
    return DagType::ref_face_type();
  else if (CubitUtil::compare(name,"curve"))
    return DagType::ref_edge_type();
  else if (CubitUtil::compare(name,"vertex"))
    return DagType::ref_vertex_type();
  else
    return DagType();
}

#ifdef CAT
const char* RefEntity::get_ref_class_name(const type_info& ref_type)
{
  if( ref_type == typeid(RefGroup) )
     return RefGroup::get_class_name();
  else if( ref_type == typeid(Body) )
     return Body::get_class_name();
  else if( ref_type == typeid(RefVolume) )
     return RefVolume::get_class_name();
  else if( ref_type == typeid(RefFace) )
     return RefFace::get_class_name();
  else if( ref_type == typeid(RefEdge) )
     return RefEdge::get_class_name();
  else if( ref_type == typeid(RefVertex) )
     return RefVertex::get_class_name();
  else
     return NULL;
}
#endif


void RefEntity::color(int value)
{
  mColor = value;
  CubitObserver::notify_static_observers(this, ENTITY_GEOMETRY_COLOR_CHANGED);
}

int RefEntity::color() const
{
  return mColor;
}


