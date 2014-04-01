//-------------------------------------------------------------------------
// Filename      : RefVertex.cpp
//
// Purpose       : This file contains the implementation of the class 
//                 RefVertex. A RefVertex is a TopologyEntity and represents
//                 a point or vertex in the topology of the Model. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96 
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------


#include <assert.h>

#include "CubitDefines.h"
#include "GeometryDefines.h"

#include "RefEntityFactory.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryQueryEngine.hpp"


#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "Point.hpp"
#include "DLIList.hpp"
#include "CastTo.hpp"
#include "ToolData.hpp"
#include "CoVertex.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor with a pointer to a Point
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------
RefVertex::RefVertex(Point* pointPtr)
{
   // Set the GeometryEntity pointer
   if (pointPtr != NULL)
   {
      set_geometry_entity_ptr(pointPtr) ;
   }
   else
   {
      PRINT_ERROR("In the RefVertex(Point*) constructor\n");
      PRINT_ERROR("       Input Point pointer is NULL\n");
      assert(CUBIT_FALSE);
   }
   
   // Initialize the member data
   initialize();
}

//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 11/06/96
//-------------------------------------------------------------------------
RefVertex::~RefVertex()
{
}

//-------------------------------------------------------------------------
// Purpose       : Return a pointer to the point associated with a vertex.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
Point* RefVertex::get_point_ptr() 
{
  return CAST_TO(get_geometry_entity_ptr(), Point) ;
}

const Point* RefVertex::get_point_ptr() const 
{
  return CAST_TO(get_geometry_entity_ptr(), Point) ;
}


//-------------------------------------------------------------------------
// Purpose       : Returns the global coordinates of this RefVertex.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 09/25/96
//-------------------------------------------------------------------------
CubitVector RefVertex::coordinates() const
{
  const Point* point_ptr = get_point_ptr();
  
  if (point_ptr != NULL)
  {
    return point_ptr->coordinates();
  }
  else
  {
    PRINT_ERROR("In RefVertex::coordinates\n"
                "       RefVertex %d has no GeometryEntity attached to it.\n"
                "       THIS IS A BUG - PLEASE REPORT IT.\n", id());
    assert( point_ptr != NULL );
    return CubitVector();
  }
}

//-------------------------------------------------------------------------
// Purpose       : Return the "center point" of the RefVertex -- its 
//                 coordinates.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/30/96
//-------------------------------------------------------------------------
CubitVector RefVertex::center_point()
{
   return this->coordinates();
}
void RefVertex::common_ref_edges(RefVertex *other_vertex,
                                DLIList <RefEdge*> &common_ref_edges,
                                RefFace *owning_face) 
{
    //- returns an edge sharing the other vertex and owned by owning_face 
    //- (if non-NULL)
    //- The edge must have this vertex and other_vertex as its
    //- start and end vertices.
  DLIList<RefEntity*> temp_list;
  join(other_vertex, temp_list);

  int i;
  for (i = temp_list.size(); i > 0; i--)
  {
    RefEdge *common_edge = CAST_TO(temp_list.get(), RefEdge);
      //This extra 'if' block is needed in case other_vertex == this
    if (common_edge &&  
        ((common_edge->start_vertex() == other_vertex &&
           common_edge->end_vertex() == this) ||
         ( common_edge->start_vertex() == this &&
           common_edge->end_vertex() == other_vertex )))
    {
      if (common_edge && (!owning_face || common_edge->is_child(owning_face)))
        common_ref_edges.append(common_edge);
    }
    temp_list.step();
  }
  return;
}

RefEdge *RefVertex::common_ref_edge(RefVertex *other_vertex, 
                                    RefFace *owning_face) 
{
    //- returns an edge sharing the other vertex and owned by owning_face
    //- (if non-NULL)
    //- The edge must have this vertex and other_vertex as its
    //- start and end vertices.
  DLIList<RefEntity*> temp_list;
  join(other_vertex, temp_list);

  int i;
  for (i = temp_list.size(); i > 0; i--) {
    RefEdge *common_edge = CAST_TO(temp_list.get(), RefEdge);
      //This extra 'if' block is needed in case other_vertex == this
    if ( ( common_edge->start_vertex() == other_vertex &&
           common_edge->end_vertex() == this) ||
         ( common_edge->start_vertex() == this &&
           common_edge->end_vertex() == other_vertex ) )
    {
      if (common_edge && (!owning_face || common_edge->is_child(owning_face)))
        return common_edge;
    }
    else
      temp_list.step();
  }
  
  return NULL;
}

//-------------------------------------------------------------------------
// Purpose       : Spatially compare two RefVertexes, notification is made if
//                 the flag notify_ref_entity is set.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 04/08/97
//-------------------------------------------------------------------------
CubitBoolean RefVertex::about_spatially_equal(
                                RefVertex* ref_vertex_ptr_2,
                                double tolerance_factor,
                                CubitBoolean notify_ref_entity )
{
   if( this == ref_vertex_ptr_2 )
   {
      if (notify_ref_entity)
        remove_compare_data();
      return CUBIT_TRUE;
   }
 
   const CubitVector vertex_1_position = this->coordinates();
   const CubitVector vertex_2_position = ref_vertex_ptr_2->coordinates();
   if (!GeometryQueryTool::instance()-> about_spatially_equal(
     vertex_1_position, vertex_2_position, tolerance_factor))
     return CUBIT_FALSE;
     
   if ( notify_ref_entity == CUBIT_TRUE )
     this->notify(ref_vertex_ptr_2, COMPARISON_FOUND);

   return CUBIT_TRUE;
}  

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********

//-------------------------------------------------------------------------
// Purpose       : Initializes member data
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/07/96
//-------------------------------------------------------------------------
void RefVertex::initialize()
{
   GeometryEntity* geom_ptr = get_geometry_entity_ptr();
   int saved_id = geom_ptr->get_saved_id();
   if ( !saved_id || RefEntityFactory::instance()->get_ref_vertex(saved_id) )
   {
     saved_id =  RefEntityFactory::instance()->next_ref_vertex_id();
     geom_ptr->set_saved_id(saved_id);
   }
   entityId = saved_id;

   // Initialize some attributes

     // read and initialize attributes
   auto_read_cubit_attrib();
   auto_actuate_cubit_attrib();
   
   // Assign a default entity name
   assign_default_name();

}
