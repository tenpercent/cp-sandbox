//-------------------------------------------------------------------------
// Filename      : Body.cpp
//
// Purpose       : This file contains the implementation of the class
//                  Body. 
//
// Special Notes :
//
// Creator       : 
//
// Creation Date : 
//
// Owner         : 
//-------------------------------------------------------------------------

#include <assert.h>
#include "CubitDefines.h"
#include "CubitMessage.hpp"
#include "CubitBox.hpp"

#include "BodySM.hpp"
#include "Body.hpp"

#include "CoVolume.hpp"
#include "CoEdge.hpp"

#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "DLIList.hpp"
#include "RefEntityFactory.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryQueryEngine.hpp"
#include "CastTo.hpp"

//-------------------------------------------------------------------------
// Purpose       : Default constructor.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//-------------------------------------------------------------------------
Body::Body()
{
}

//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to an TopologyBridge
//                 as input.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//-------------------------------------------------------------------------
Body::Body(BodySM* OSMEPtr)
{
   // Set the new Body's OSME pointer
   if (OSMEPtr != NULL)
   { 
      set_topology_bridge(OSMEPtr) ;
   }
   else
   {
      PRINT_ERROR("In the Body(OSME*) constructor\n");
      PRINT_ERROR("       Input OSME pointer is NULL\n");
      assert(OSMEPtr != NULL);
   }

   // Set the Entity ID for this new Body
   entityId = RefEntityFactory::instance()->next_body_id();

   copied_from_body_id = 0;

     // read and initialize attributes
   auto_read_cubit_attrib();
   auto_actuate_cubit_attrib();

     // Assign a default entity name
   assign_default_name();   
}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------
Body::~Body()
{
}


//-------------------------------------------------------------------------
// Purpose       : Get BodySM ptr
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/23/03
//-------------------------------------------------------------------------
BodySM* Body::get_body_sm_ptr() const
{
  TopologyBridge* bridge = bridge_manager()->topology_bridge();
  return dynamic_cast<BodySM*>(bridge);
}

//* Method: copiedFromId
//** Sets the Id of the body that this body was copied from
//*==

void Body::copiedFromId( int Id )
{
    copied_from_body_id = Id;
}

//* Method: copiedFromId
//** Returns the Id of the body that this body was copied from
//*==

int Body::copiedFromId ()
{
    return copied_from_body_id;
}

/* This is wrong for misc. arrangements of non-manifold topology
 * -- j.kraftcheck
   
CubitBoolean Body::is_sheet_body()
{
    //This function just checks for bodies which have a single coedge
    // on a RefEdge, which would indicate an open sheet body.  This
    // function will not find sheet bodies which are closed...
  DLIList<CoEdge*> co_edges;
  this->co_edges( co_edges );
  DLIList<RefEdge*> edge_list_1, edge_list_2;
  int i;
  for( i = co_edges.size(); i > 0; i-- )
  {
    CoEdge* this_co_edge = co_edges.get_and_step();
    RefEdge* this_ref_edge = this_co_edge->get_ref_edge_ptr();
    if( edge_list_1.is_in_list( this_ref_edge ) )
    {
      edge_list_2.append( this_ref_edge );
    }
    else
    {
      edge_list_1.append_unique( this_ref_edge );
    }
  }
  
    //subtract edge_list_2 from edge_list_1
  DLIList<RefEdge*> temp_list = edge_list_1;
  for( i = temp_list.size(); i > 0; i-- )
  {
    RefEdge* this_edge = temp_list.get_and_step();
    if( edge_list_2.is_in_list( this_edge ) )
    {
      edge_list_1.remove( this_edge );
    }
  }
  
  if( edge_list_1.size() > 0 )
  {
    for( i = edge_list_1.size(); i > 0; i-- )
    {
      if( edge_list_1.get_and_step()->geometry_type() != POINT_CURVE_TYPE )
         return CUBIT_TRUE;
    }
  }
  return CUBIT_FALSE;
}
 */

//-------------------------------------------------------------------------
// Purpose       : Check if this body contains only sheet volumes.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/15/04
//-------------------------------------------------------------------------
CubitBoolean Body::is_sheet_body()
{
  DLIList<RefVolume*> volumes;
  ref_volumes(volumes);
  while (volumes.size())
    if (!volumes.pop()->is_sheet())
      return CUBIT_FALSE;
  return CUBIT_TRUE;
}
 
CubitVector Body::center_point()
{
  return bounding_box().center();
}

CubitPointContainment Body::point_containment( CubitVector &point )
{
  return get_body_sm_ptr()->point_containment( point );
}

CubitBoolean Body::get_mass_props(CubitVector &cofg)
{
  double vol;
  CubitStatus result = get_body_sm_ptr()->mass_properties(cofg, vol);
  return result ? CUBIT_TRUE : CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes : This function assumes that a Body's spatial extent is
//                 defined by the set of RefVolumes it "owns".
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/23/96
//-------------------------------------------------------------------------
CubitBox Body::bounding_box()
{
   // Get the list of RefVolumes. Perform the union of the bounding
   // boxes of all the RefVolumes to obtain the bounding box of 
   // the Body.
   DLIList<RefVolume*> ref_volume_list ;
   /*CubitStatus status = */
   ref_volumes(ref_volume_list);

   if ( ref_volume_list.size() == 0 )
   {
     CubitVector vec(0.0, 0.0, 0.0);
     return CubitBox( vec, vec );
   }
   RefVolume* ref_volume_ptr = NULL ;

   ref_volume_list.reset() ;
   CubitBox result_box = ref_volume_list.get_and_step()->bounding_box();
   for(int i = 1 ; i < ref_volume_list.size() ; i++)
   {
      ref_volume_ptr = ref_volume_list.get_and_step() ;
      result_box |= ref_volume_ptr->bounding_box() ;
   }

   return result_box ;
}


double Body::measure()
{
  DLIList<RefVolume*> volumes;
  ref_volumes(volumes);
  double volume = 0.0;
  for (int i = volumes.size(); i > 0; i--)
  {
     volume += volumes.get_and_step()->measure();
  }
  
  return volume;
}

CubitString Body::measure_label()
{
  return "volume";
}

int Body::validate()
{
  int error = 0;

      // Perform general RefEntity checks (measure > 0)
  error += RefEntity::validate();
  
    // check the body from acis
  BodySM* osme_ptr = get_body_sm_ptr();
  DLIList <TopologyEntity*> bad_entities;
  if ( osme_ptr != NULL )
  {
    
    error += osme_ptr->validate( entity_name(), bad_entities);
  }
  else 
  {
    PRINT_WARNING("\tWARNING: Null underlying solid modeling body for %s, (%s %d)\n",
                  entity_name().c_str(), class_name(), id());
    error++;
  }
  return error;
}

void Body::color(int value)
{
  int i;
  DLIList<RefVolume*> refVolumeList;
  ref_volumes ( refVolumeList );

  refVolumeList.reset();
  for (i = refVolumeList.size(); i--; )
    refVolumeList.get_and_step()->color(value);
}

int Body::color() const
{
  return RefEntity::color();
}


