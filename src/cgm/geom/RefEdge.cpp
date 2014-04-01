//-------------------------------------------------------------------------
// Filename      : RefEdge.cpp
//
// Purpose       : This file contains the implementation of the class 
//                 RefEdge. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96 
//
// Owner         :  Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"

#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "Chain.hpp"

#include "Curve.hpp"
#include "CoEdge.hpp"
#include "Loop.hpp"

#include "RefEntityFactory.hpp"
#include "GeometryQueryTool.hpp"

#include "CastTo.hpp"
#include "DLIList.hpp"

#include "CubitAttrib.hpp"
#include "CubitUtil.hpp"

#include "ToolData.hpp"
#include "CpuTimer.hpp"

bool RefEdge::mSuppressEdgeLengthWarning = false;

//-------------------------------------------------------------------------
// Purpose       : Constructor with a pointer to a Curve
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------
RefEdge::RefEdge(Curve* curvePtr)
{
    // Set the GeometryEntity pointer
  if (curvePtr != NULL)
  {
    set_geometry_entity_ptr(curvePtr) ;
  }
  else
  {
    PRINT_ERROR("In the RefEdge(Curve*) constructor\n");
    PRINT_ERROR("       Input Curve pointer is NULL\n");
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
// Creator       : Raikanta Sahu
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------
RefEdge::~RefEdge()
{
}

//-------------------------------------------------------------------------
// Purpose       : Return a pointer to the curve associated with a edge.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
Curve* RefEdge::get_curve_ptr() 
{
  return CAST_TO(get_geometry_entity_ptr(), Curve) ;
}

const Curve* RefEdge::get_curve_ptr() const 
{
  return CAST_TO(get_geometry_entity_ptr(), Curve) ;
}

//-------------------------------------------------------------------------
// Purpose       : Get Chain associated with this RefEdge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/22/03
//-------------------------------------------------------------------------
Chain* RefEdge::get_chain_ptr() 
{
  GroupingEntity* gpe_ptr = get_first_grouping_entity_ptr();
  if ( !gpe_ptr ||        // no chain
       gpe_ptr->next() )  // multiple chains
  {
    return 0;
  }
  
  return dynamic_cast<Chain*>(gpe_ptr);
}

CubitStatus RefEdge::get_point_direction( CubitVector& origin, CubitVector& direction )
{
   Curve* curve_ptr = get_curve_ptr();

   if( curve_ptr != NULL )
   {

      if( curve_ptr->geometry_type() != STRAIGHT_CURVE_TYPE ) 
         return CUBIT_FAILURE;

      if( curve_ptr->get_point_direction( origin, direction ) == CUBIT_FAILURE )
      {
         origin = start_coordinates();
         direction = end_coordinates() - origin;
         direction.normalize();
      }
   }
   else 
   {
      PRINT_WARNING("In RefEdge::get_point_direction\n"
                    "         %s (curve %d) is not associated with a valid\n"
                    "         underlying geometric Curve\n",
                    entity_name().c_str(), id()) ;
      return CUBIT_FAILURE ;
   }
   
   if (curve_ptr->bridge_sense() == CUBIT_REVERSED)
     direction = -direction;
   
   return CUBIT_SUCCESS;

}

CubitStatus RefEdge::get_center_radius( CubitVector& center, double& radius )
{
   Curve* curve_ptr = get_curve_ptr();

   if( curve_ptr != NULL )
   {

      if( curve_ptr->geometry_type() != ARC_CURVE_TYPE &&
          curve_ptr->geometry_type() != ELLIPSE_CURVE_TYPE )
         return CUBIT_FAILURE;

      if( curve_ptr->get_center_radius( center, radius ) == CUBIT_FAILURE )
         return CUBIT_FAILURE;
   }
   else 
   {
      PRINT_WARNING("In RefEdge::get_center_radius\n"
                    "         %s (curve %d) is not associated with a valid\n"
                    "         underlying geoemtric Curve\n",
                    entity_name().c_str(), id()) ;
      return CUBIT_FAILURE ;
   }
   
   return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Finds the closest point on the RefEdge to the input
//                 point and modifies the coordinate values of the input
//                 point to be those of the closest point.
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/25/97
//-------------------------------------------------------------------------
void RefEdge::move_to_curve( CubitVector& vector )
{
    // Get the Curve associated with this RefEdge
  Curve* curvePtr = this->get_curve_ptr();
  
    // Move the point to the Curve (the following call modifes the values
    // of the input "vector", if necessary)
  CubitVector closest_point;
  curvePtr->closest_point(vector, closest_point);
  vector.set( closest_point.x(),
              closest_point.y(),
              closest_point.z() );
}

CubitStatus RefEdge::get_interior_extrema(DLIList<CubitVector*>& interior_points,
                                          CubitSense& return_sense) const
{
  CubitStatus result;
  
    // Cast away the constness of the Curve
  Curve* curve = const_cast<Curve*>(get_curve_ptr());
  result = curve->get_interior_extrema(interior_points, return_sense);
  if (curve->bridge_sense() == CUBIT_REVERSED)
    return_sense = CubitUtil::opposite_sense(return_sense);
  return result;
}

CubitStatus RefEdge::closest_point( CubitVector const& location, 
                                    CubitVector& closest_location,
                                    CubitVector* tangent_ptr,
                                    CubitVector* curvature_ptr)
{
  Curve *curve_ptr = get_curve_ptr();
  CubitStatus result = curve_ptr->closest_point(location, closest_location, 
                                           tangent_ptr, curvature_ptr);
  if (tangent_ptr && curve_ptr->bridge_sense() == CUBIT_REVERSED)
    *tangent_ptr = -(*tangent_ptr);
    
  return result;
}

CubitStatus RefEdge::closest_point_trimmed( CubitVector const& location, 
                                    CubitVector& closest_location)
{
  Curve *curve_ptr = get_curve_ptr();
  return curve_ptr->closest_point_trimmed(location, closest_location);
}

CubitPointContainment RefEdge::point_containment( const CubitVector &point )
{
  Curve *curve_ptr = get_curve_ptr();
  return curve_ptr->point_containment(point);
}


//-------------------------------------------------------------------------
// Purpose       : These functions return the start and end global coordinate
//                 locations of the RefEdge.
//
// Special Notes :
// These coordinates remain consistent throughout the life of the
// RefEdge, regardless of the fact that the actual start and
// end RefVerex'es may change (e.g., after a merge operation).
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 07/03/97
//-------------------------------------------------------------------------
CubitVector RefEdge::start_coordinates()
{
  return this->start_vertex()->coordinates();
}

CubitVector RefEdge::end_coordinates()
{
  return this->end_vertex()->coordinates();
}

void RefEdge::tangent ( const CubitVector &point, CubitVector& tangent_vec )
{
    // Get the tangent at the point closest to "point" on the RefEdge.
    // The tangent is always pointing in the positive direction of the
    // RefEdge
  Curve* curve_ptr = this->get_curve_ptr();
  CubitVector closest_point;
  curve_ptr->closest_point( point, closest_point, &tangent_vec );
  if (curve_ptr->bridge_sense() == CUBIT_REVERSED)
    tangent_vec = -tangent_vec;
}
//-------------------------------------------------------------------------
// Purpose       : Get the correct tangent with respect to
//                 the ref_face_ptr.
// Special Note  : This tangent function is not the safest method for getting
//               : the tangent on a surface.  In face if there could be
//                 another direction, this function could fail...(assert).
//                 This function assumes that there is only 1 co-edge per
//                 this ref_face for this REfEdge.
//
// Creator       : David White
//
// Creation Date : 03/15/97
//-------------------------------------------------------------------------
CubitStatus RefEdge::tangent( const CubitVector &point, 
                              CubitVector& tangent_vec,
                              RefFace *ref_face_ptr )
{
    //Get the tangent for this edge.
  tangent( point, tangent_vec );
  
    //Now allign the tangent with the face.
    //For this function we assume that there must only be one
    //co edge for this ref_edge.
  DLIList<CoEdge*> co_edge_list;
  get_co_edges( co_edge_list, ref_face_ptr );
 
  assert( co_edge_list.size() == 1 ); 
  
  if ( co_edge_list.get()->get_sense() == CUBIT_REVERSED )
      tangent_vec = -tangent_vec;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Get the correct tangent with respect to
//                 the ref_face_ptr and the next_ref_edge
// Special Note  : For this function next_ref_edge must follow RefEdge.
//               : This tangent function is the safest method for getting
//               : the tangent on a surface.  The others are less restrictive
//               : about the tangents direction...  (forward, reverse)
//
// Creator       : David White
//
// Creation Date : 03/15/97
//-------------------------------------------------------------------------
CubitStatus RefEdge::tangent( const CubitVector &point, 
                              CubitVector& tangent_vec,
                              RefEdge *next_ref_edge, 
                              RefFace *ref_face_ptr )
{
  CoEdge *co_edge_this = NULL;
  CoEdge *co_edge_next = NULL;
    //First get the two coedges that corrispond to
    //this ref_edge an the next one, with reference to
    //the ref_face_ptr.
  CubitStatus status = get_two_co_edges( next_ref_edge,
                                         ref_face_ptr,
                                         co_edge_this,
                                         co_edge_next );
  if (status == CUBIT_FAILURE )
      return status;
  
  assert(co_edge_this != NULL );
  assert(co_edge_next != NULL );
  
    //Now get the tangent from this curve.
  tangent ( point, tangent_vec );
  
    //with the go_edge data we have, we can get the tangent
    //going in the right direction...
  if ( co_edge_this->get_sense() == CUBIT_REVERSED )
      tangent_vec = -tangent_vec;
  return CUBIT_SUCCESS;
}

double RefEdge::get_arc_length ()
{
    // Get the Curve associated with this RefEdge
  Curve* curve_ptr = get_curve_ptr();
  
  assert(curve_ptr != NULL);
  
    // Get the length
  return curve_ptr->get_arc_length();
}

double RefEdge::get_arc_length ( const CubitVector &point1,
                                 const CubitVector &point2 )
{
    // Get the Curve associated with this RefEdge
  Curve* curve_ptr = get_curve_ptr();
  
  assert(curve_ptr != NULL);
  
    // Get the length between the 2 points
  return curve_ptr->get_arc_length(point1, point2);
}

double RefEdge::get_arc_length ( const CubitVector &point1,
                                 int which_end )
{
  Curve* curve_ptr = get_curve_ptr();
  
  assert (curve_ptr != NULL);
  
  if (curve_ptr->bridge_sense() == CUBIT_REVERSED)
    which_end = 1 - which_end;
  
  return curve_ptr->get_arc_length(point1, which_end);
}

double RefEdge::get_chord_length()
{
  CubitVector start_pos(start_vertex()->coordinates());
  CubitVector end_pos(start_vertex()->coordinates());
  double distance = (start_pos - end_pos).length();
  return distance;
}


//-------------------------------------------------------------------------
// Purpose       : Return the "actual" center point (the midpoint along 
//                 the arc length) of the RefEdge.
//
// Special Notes :
//
// Creator       : Malcolm Panthaki
//
// Creation Date : 2/27/97
//-------------------------------------------------------------------------
CubitVector RefEdge::center_point()
{
    // Get the parameter range of the RefEdge
  Curve* curve_ptr = get_curve_ptr();
  
  assert (curve_ptr != NULL);
  
  return curve_ptr->center_point();
}

CubitStatus RefEdge::mid_point ( const CubitVector &point1,
                                 const CubitVector &point2,
                                 CubitVector& mid_point )
{
    // Get the Curve associated with this RefEdge
  Curve* curve_ptr = get_curve_ptr();
  
  assert(curve_ptr != NULL);
  
    // Get the global location of parameter3
  return curve_ptr->mid_point(point1, point2, mid_point);
  
}

CubitStatus RefEdge::mid_point (  CubitVector& mid_point )
{
    // Get the Curve associated with this RefEdge
  Curve* curve_ptr = get_curve_ptr();
  
  assert(curve_ptr != NULL);
  
    // Get the global location of parameter3
  return curve_ptr->mid_point(mid_point);
  
}
CubitStatus RefEdge::position_from_fraction( double fraction_along_curve,
                                             CubitVector& output_position )
{
  
    //Get the Curve of this RefEdge.
  Curve* curve_ptr = this->get_curve_ptr();
  
  assert( fraction_along_curve < 1.0000001 && fraction_along_curve > -0.0000001 );
    //Now get the postion from this fraction value.
  
  if (curve_ptr->bridge_sense() == CUBIT_REVERSED)
    fraction_along_curve = 1.0 - fraction_along_curve;  
    
  CubitStatus result = curve_ptr->position_from_fraction( fraction_along_curve,
                                                          output_position );
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Get the parameter value of the start of the RefEdge. 
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/27/97
//-------------------------------------------------------------------------
double RefEdge::start_param()
{
    // Get the Curve associated with this RefEdge
  Curve* curve_ptr = get_curve_ptr();
  
  assert(curve_ptr != NULL);
  
  if (curve_ptr->bridge_sense() == CUBIT_REVERSED)
    return -(curve_ptr->end_param());
  else
    return curve_ptr->start_param();
}

//-------------------------------------------------------------------------
// Purpose       : Get the parameter value of the end of the RefEdge. 
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 2/27/97
//-------------------------------------------------------------------------
double RefEdge::end_param()
{
    // Get the Curve associated with this RefEdge
  Curve* curve_ptr = get_curve_ptr();
  
  assert(curve_ptr != NULL);
  
  if (curve_ptr->bridge_sense() == CUBIT_REVERSED)
    return -(curve_ptr->start_param());
  else
    return curve_ptr->end_param();
}

CubitBoolean RefEdge::get_param_range( double& start_param, double& end_param )
{
    // Get the Curve of this RefEdge
  Curve* curvePtr = this->get_curve_ptr();
  
    // Now get the parameter values of the start and end locations
    // of this RefEdge
  CubitBoolean result = curvePtr->get_param_range( start_param, end_param );
  
  if (curvePtr->bridge_sense() == CUBIT_REVERSED)
  {
    double tmp_start_param = start_param;
    start_param = -end_param;
    end_param = -tmp_start_param;
  }
  
  return result;
}

double RefEdge::u_from_position (const CubitVector& input_position)
{
    // Get the Curve of this RefEdge
  Curve* curvePtr = this->get_curve_ptr();
  
    // Now get the parameter values of the start and end locations
    // of this RefEdge
  double param = curvePtr->u_from_position(input_position);
  if (curvePtr->bridge_sense() == CUBIT_REVERSED)
    param = -param;
  
  return param;
}

CubitStatus RefEdge::position_from_u (double u_value,
                                      CubitVector& output_position)
{
    // Get the Curve of this RefEdge
  Curve* curvePtr = this->get_curve_ptr();
  
  if (curvePtr->bridge_sense() == CUBIT_REVERSED)
    u_value = -u_value;
  
    // Now get the parameter values of the start and end locations
    // of this RefEdge
  return curvePtr->position_from_u(u_value, output_position);
}

double RefEdge::u_from_arc_length ( double root_param, double arc_length )
{
    // Get the Curve of this RefEdge
  Curve* curvePtr = this->get_curve_ptr();
  
  if (curvePtr->bridge_sense() == CUBIT_REVERSED)
    return -(curvePtr->u_from_arc_length(-root_param, -arc_length));
  else
    return curvePtr->u_from_arc_length(root_param, arc_length);
}

double RefEdge::fraction_from_arc_length(RefVertex *root_vertex,
                                         double     length)
{
  if (root_vertex != start_vertex() && root_vertex != end_vertex())
    return -1.0;
 
  if (geometry_type() == POINT_CURVE_TYPE || get_arc_length() < GEOMETRY_RESABS)
    return 0.0;
 
  if (length >= get_arc_length())
    return 1.0;

  if (root_vertex == start_vertex())
    return length/get_arc_length();

  else
    return 1-length/get_arc_length();
}

CubitStatus RefEdge::point_from_arc_length ( const CubitVector& root_point, 
                                             double arc_length,
                                             CubitVector& new_point )
{
    // Get the Curve of this RefEdge
  Curve* curvePtr = this->get_curve_ptr();
  
  if (curvePtr->bridge_sense() == CUBIT_REVERSED)
    arc_length = -arc_length;
    
    // Now get the parameter values of the start and end locations
    // of this RefEdge
  return curvePtr->point_from_arc_length (root_point, arc_length,
                                          new_point );
}

CubitStatus RefEdge::point_from_arc_length ( double root_param, 
                                             double arc_length,
                                             CubitVector& new_point )
{
    // Get the Curve of this RefEdge
  Curve* curvePtr = this->get_curve_ptr();
  
  if (curvePtr->bridge_sense() == CUBIT_REVERSED)
    arc_length = -arc_length;
    
    // Now get the parameter values of the start and end locations
    // of this RefEdge
  return curvePtr->point_from_arc_length (root_param, arc_length,
                                          new_point );
}

double RefEdge::length_from_u( double parameter1,
                               double parameter2 )
{
    // Get the Curve of this RefEdge
  Curve* curvePtr = this->get_curve_ptr();
  
  if (curvePtr->bridge_sense() == CUBIT_REVERSED)
    return -(curvePtr->length_from_u(-parameter1, -parameter2));
  else
    return curvePtr->length_from_u(parameter1, parameter2);
}

CubitBoolean RefEdge::is_periodic( )
{
  double temp_val;
  return this->is_periodic(temp_val);
}

CubitBoolean RefEdge::is_periodic( double& period)
{
    // Get the Curve of this RefEdge
  Curve* curvePtr = this->get_curve_ptr();
  
    // Now get the parameter values of the start and end locations
    // of this RefEdge
  CubitBoolean periodic = curvePtr->is_periodic(period);
  
  if (curvePtr->bridge_sense() == CUBIT_REVERSED)
    period = -period;
  
  return periodic;
}

//-------------------------------------------------------------------------
// Purpose       : Get parent CoEdges
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/29/03
//-------------------------------------------------------------------------
CubitStatus RefEdge::get_co_edges( DLIList<CoEdge*>& co_edges_found_list,
                                   RefFace *input_ref_face_ptr )
{
  for (SenseEntity* coedge_ptr = get_first_sense_entity_ptr();
       coedge_ptr;
       coedge_ptr = coedge_ptr->next_on_bte())
  {
    if (!input_ref_face_ptr ||
      coedge_ptr->get_parent_basic_topology_entity_ptr() == input_ref_face_ptr)
    {
      co_edges_found_list.append (dynamic_cast<CoEdge*>(coedge_ptr));
    }
  }
  
  return CUBIT_SUCCESS;
}

double RefEdge::angle_between( RefEdge *other_edge_ptr,
                               RefFace *face_ptr )
{
  DLIList<DLIList<RefEdge*>*> ref_edge_loop;
  DLIList<RefEdge*> *ref_edge_list;
  CubitVector vertex_vector, normal_vector;
  CubitVector left_vector, right_vector;

    // Loop through face's loops.
  face_ptr->ref_edge_loops( ref_edge_loop );
  double return_val = CUBIT_DBL_MAX;
  
  for( int i = ref_edge_loop.size(); i>0; i-- )
  {
      // Look for this edge in the list.
    ref_edge_list = ref_edge_loop.get_and_step();
    if( !ref_edge_list->move_to( this ) )
      continue;
      // Look for other_edge in the list.
    if( !ref_edge_list->move_to( other_edge_ptr ) )
    {
      PRINT_ERROR("Attempted to get angle between two "
                  "RefEdges not on the same Loop.\n");
      return 0.;
    }

      // If other_edge is before this edge in the loop . . .
    if( ref_edge_list->next() == this )
    {   // (other = left, this = right)
        // Get the vector of the the common vertex.
      vertex_vector =
        other_edge_ptr->common_ref_vertex(this, face_ptr)->coordinates();
        // Get the normal vector of the face.
      normal_vector = face_ptr->normal_at( vertex_vector );
        // Get the vectors of both edges.
      other_edge_ptr->tangent( vertex_vector, left_vector,
                               this, face_ptr );
      this->tangent( vertex_vector, right_vector,
                     ref_edge_list->next(2), face_ptr);
    }

      // If this edge is before other_edge in the loop . . .
    else if( ref_edge_list->prev() == this )
    {   // (this = left, other = right)
        // Get the vector of the the common vertex.
      vertex_vector =
          common_ref_vertex(other_edge_ptr,face_ptr)->coordinates();
        // Get the normal vector of the face.
      normal_vector = face_ptr->normal_at( vertex_vector );
        // Get the vectors of both edges.
      this->tangent( vertex_vector, left_vector,
                     other_edge_ptr, face_ptr);
      other_edge_ptr->tangent( vertex_vector, right_vector,
                               ref_edge_list->next(), face_ptr );
    }

      // Otherwise, error.
    else
    {
      PRINT_ERROR("Attempted to get angle between two "
                  "non-consecutive edges.\n");
      return 0.;
    }

      // Return the angle between the vectors.
    return_val = normal_vector.vector_angle( right_vector, -left_vector );
    break;
  }
    // If this far, then error.
  if (CUBIT_DBL_MAX == return_val) {
    PRINT_ERROR("Attempted to get angle between two edges "
                "where one is not found on the given face.\n");
    return_val = 0.0;
  }
  
  return return_val;
}  

int RefEdge::dimension() const
{
  return 1;
}

double RefEdge::measure()
{
  return get_arc_length();
}

CubitString RefEdge::measure_label()
{
  return "length";
}

int RefEdge::num_of_common_ref_face( RefEdge *other_edge )
{
  DLIList<RefFace*> ref_faces, ref_faces_for_other_edge;
  this->ref_faces(ref_faces);
  other_edge->ref_faces(ref_faces_for_other_edge);
                                                                                
  int i, j;
  int num = 0;
  for (i = 0; i < ref_faces.size(); i++)
  {
     RefFace * ref_face = ref_faces.get_and_step();
     for (j = 0; j < ref_faces_for_other_edge.size(); j++)
        if (ref_face == ref_faces_for_other_edge.get_and_step())
           num++;
  }
                                                                                
  return num;
}
 
// Get one common RefFace between two edges.
RefFace *RefEdge::common_ref_face( RefEdge *other_edge )
{
  DLIList<RefFace*> ref_faces, ref_faces_for_other_edge;
  this->ref_faces(ref_faces);
  other_edge->ref_faces(ref_faces_for_other_edge);

  int i, j;
  for (i = 0; i < ref_faces.size(); i++)
  {
     RefFace * ref_face = ref_faces.get_and_step();
     for (j = 0; j < ref_faces_for_other_edge.size(); j++)
        if (ref_face == ref_faces_for_other_edge.get_and_step())
           return ref_face;
  }

  return NULL;
}

int RefEdge::common_ref_faces ( RefEdge* input_edge, DLIList<RefFace*> &common_face_list )
{
   int nedges = 0;
   DLIList<RefFace*> ref_faces, ref_faces_for_other_edge;
   this->ref_faces(ref_faces);
   input_edge->ref_faces(ref_faces_for_other_edge);                            

   int i, j;
   for (i = 0; i < ref_faces.size(); i++)
   {
      RefFace * ref_face = ref_faces.get_and_step();
      for (j = 0; j < ref_faces_for_other_edge.size(); j++)
         if (ref_face == ref_faces_for_other_edge.get_and_step())
         {
           common_face_list.append(ref_face);
           nedges++;
         }
   }

   return nedges;
}

RefVertex *RefEdge::common_ref_vertex( RefEdge *other_edge )
{
  RefVertex *this_start = start_vertex();
  RefVertex *this_end   = end_vertex();
  RefVertex *other_start = other_edge->start_vertex();
  RefVertex *other_end   = other_edge->end_vertex();
  
  if ( this_start == other_start || this_start == other_end )
  {
    return this_start;
  }
  else if ( this_end == other_start || this_end == other_end )
  {
    return this_end;
  }
  else
  {
    return NULL;
  }
}
CubitBoolean RefEdge::common_vertices( RefEdge *other_edge,
                                       DLIList<RefVertex*> &common_verts)
{
  CubitBoolean result = CUBIT_FALSE;
  RefVertex *this_start = start_vertex();
  RefVertex *this_end   = end_vertex();
  RefVertex *other_start = other_edge->start_vertex();
  RefVertex *other_end   = other_edge->end_vertex();
  
  if ( this_start == other_start || this_start == other_end )
  {
    common_verts.append(this_start);
    result = CUBIT_TRUE;
  }
  else if ( this_end == other_start || this_end == other_end )
  {
    common_verts.append(this_end);
    result = CUBIT_TRUE;
  }

  return result;
}
//-------------------------------------------------------------------------
// Purpose       : Return a pointer to the ref_vertex that is common
//                 between the two ref edges and is in the correct order
//                 going with respect to the face.
//
// Creator       : David White
//
// Creation Date : 03/15/97
//-------------------------------------------------------------------------

RefVertex *RefEdge::common_ref_vertex( RefEdge *next_ref_edge,
                                       RefFace *ref_face_ptr )
{
  CoEdge *co_edge_this = NULL;
  CoEdge *co_edge_next = NULL;
  
    //First get the two coedges that corrispond to
    //this ref_edge an the next one, with reference to
    //the ref_face_ptr.
  CubitStatus status = get_two_co_edges( next_ref_edge,
                                         ref_face_ptr,
                                         co_edge_this,
                                         co_edge_next );
  if (status == CUBIT_FAILURE )
      return NULL;
  assert(co_edge_this != NULL );
  assert(co_edge_next != NULL );
  RefVertex *common_vertex;
  
    //Now according to the sense get the vertex at the
    //end of this edge (start if reversed).
  if ( co_edge_this->get_sense() == CUBIT_FORWARD )
      common_vertex = end_vertex();
  else
      common_vertex = start_vertex();
  
    //Lets just do a sanitiy check...
  if (common_vertex == NULL ||
      !next_ref_edge->is_directly_related( common_vertex )) {
    
      // let's check for bad sense, then print warning and return
      // correct vertex
    common_vertex = other_vertex(common_vertex);
    if (common_vertex != NULL &&
        next_ref_edge->is_directly_related( common_vertex )) {
      
      PRINT_ERROR(" bad sense between curve %d and surface %d; please"
                  " report this.\n", id(), ref_face_ptr->id());
    }
    else {
      PRINT_ERROR("unable to find common vertex (curves %d, %d)",
                  this->id(), next_ref_edge->id());
      assert(CUBIT_TRUE);
      return (RefVertex*)NULL;
    }
  }
    //Hurray, Success...
  return common_vertex;
}

RefEdge* RefEdge::get_other_curve(RefVertex* common_vertex,
                         RefFace* ref_face_ptr)
{
  DLIList<RefEdge*> curves;
  ref_face_ptr->ref_edges(curves);
  for(int ii = curves.size(); ii>0; ii--)
  {
    RefEdge* temp_edge = curves.get_and_step();
    if((temp_edge->is_directly_related(common_vertex)) &&
       (this != temp_edge))
       return temp_edge;
  }
  return NULL;
}


//-------------------------------------------------------------------------
// Purpose       : Returns the co_edge that corresponds to 'this' ref_edge
//                 and the one that corresponds to next_ref_edge, with
//                 respect to the ref_face_ptr.
// Special Notes : next_ref_edge must follow 'this' ref_edge in a Loop
//                 on the ref_face_ptr, this is assumed so the function
//                 will assert if this is not done...
//
// Creator       : David White
//
// Creation Date : 03/15/97
//-------------------------------------------------------------------------
CubitStatus RefEdge::get_two_co_edges( RefEdge *next_ref_edge,
                                       RefFace *ref_face_ptr,
                                       CoEdge *&co_edge_this,
                                       CoEdge *&co_edge_next )
{
  DLIList<Loop*> loop_list;
  Loop *loop_ptr;
  CubitStatus status;
  DLIList<CoEdge*> co_edge_list;
  
    //First get the loops for this ref_face;
  ref_face_ptr->loops( loop_list );
    //Now we want to find the coedge list that
    //has 'this' refedge followed by the next one
  assert( loop_list.size() != 0 );
  
  for ( int i = loop_list.size(); i--; )
  {
      // get the ordered co-edges of this loop
    loop_ptr = loop_list.get_and_step();
    co_edge_list.clean_out();   
    status = loop_ptr->ordered_co_edges( co_edge_list );

      //Now find the coedges corresponding to this ref_edge
      //  and the next ref_edge.
    if ( status == CUBIT_SUCCESS ) {
      for ( int j = co_edge_list.size(); j--; ) 
      {
          // candidates
        co_edge_this = co_edge_list.get_and_step();
        co_edge_next = co_edge_list.get();
        
          // really correspond to this and next edge?
        if ( co_edge_this->get_ref_edge_ptr() == this &&
             co_edge_next->get_ref_edge_ptr() == next_ref_edge ) 
        {
          return CUBIT_SUCCESS;          
        }        
      }
    }
  }  
  co_edge_this = co_edge_next = NULL;
  PRINT_ERROR("in RefEdge::get_two_co_edges.\n"
              "Couldn't find CoEdgeList with this edge\n"
              "and the next edge passed in.\n");
  return CUBIT_FAILURE;
}


RefVertex *RefEdge::other_vertex( RefVertex *vertex )
{
  if ( vertex == start_vertex() )
  { 
    return end_vertex();
  }
  else if ( vertex == end_vertex() )
  {
    return start_vertex();
  }
  else
  {
    return NULL;
  }
}

RefFace *RefEdge::other_face(RefFace *not_face, RefVolume *ref_volume) 
{
    //- return the (an) other face sharing this edge, which also borders
    //- ref_volume if non-NULL
  DLIList<RefFace*> temp_faces;
  ref_faces(temp_faces);
  int i;
  for (i = temp_faces.size(); i > 0; i--) {
    RefFace *other_face = temp_faces.get_and_step();
    if (other_face != not_face &&
        (!ref_volume || other_face->is_directly_related(ref_volume)))
      return other_face;
  }
  
  return NULL;
}


CubitStatus RefEdge::relative_sense( RefEdge *ref_edge_ptr_2, 
                                     double tolerance_factor,
                                     CubitSense *sense,
                                     CubitBoolean &spatially_equal,
                                     CubitBoolean force_merge)
{
    // It is assumed that the endpoints have been determined to be
    // spatially equal before this function is called.
  
    // Algorithm: Get a point from 'this' (the 1/3 point).  See if it
    // lies on the second RefEdge also, within tolerance.  If so,
    // spatially_equal is set to true.  Next, compare the normal directions
    // at the common point.  If they point the same direction (dot product
    // of tangents is > 0), sense is forward.  Else, sense is reversed.  If
    // the tangent dot product is between -.5 and .5, a warning is printed.
  
  const double ONE_THIRD = 1.0/3.0;
  if (sense)
    *sense = CUBIT_FORWARD;
  CubitStatus result;
  
    // Find the point 1/3 along *this*
  CubitVector test_point_1, test_point_2;
  result = this->position_from_fraction( ONE_THIRD,
                                         test_point_1);
  if ( result != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Error in RefEdge::compare(refedges).\n"
                "Can't find position 1/3 along curve.\n");
    return CUBIT_FAILURE;
  }
  
    // See if the 1/3 point on *this* lies on the other curve
  if ( ref_edge_ptr_2->closest_point_trimmed(test_point_1, test_point_2)
       != CUBIT_SUCCESS )
  {
    return CUBIT_FAILURE;
  }
  if ( GeometryQueryTool::instance()->
       about_spatially_equal(test_point_1, test_point_2,tolerance_factor ))
  {
    spatially_equal = CUBIT_TRUE;
  }
  else if ( !force_merge ) 
  {
    return CUBIT_FAILURE;
  }
  
    // Now find the sense
  if (sense)
  {
    CubitVector tangent_1, tangent_2;
    result = this->closest_point(test_point_2,
                                 test_point_1,
                                 &tangent_1);
    if (result == CUBIT_SUCCESS)
      result = ref_edge_ptr_2->closest_point(test_point_1,
                                             test_point_2,
                                             &tangent_2);
    if ( result != CUBIT_SUCCESS )
    {
      PRINT_ERROR("Error in RefEdge::relative_sense(refedges).\n"
                  "Can't find Curve tangents.\n");
      return CUBIT_FAILURE;
    }
      
      // Find the sense
			
    //If one of the curves is zero-length, it will have a zero
    //tangent vector.
    double len_product = tangent_1.length() * tangent_2.length();
    if( len_product > CUBIT_DBL_MIN )
    {
		
      double dot_product = (tangent_1 % tangent_2) / len_product;
      if (dot_product < 0)
      *sense = CUBIT_REVERSED;
    	//if (dot_product > -.5 && dot_product < .5)
      //      tangent_warning = CUBIT_TRUE;
    }
    else
    {
      //If one of the tangents is zero-length, one of the curves had
      //better be as well.
      assert( (measure() * ref_edge_ptr_2->measure()) < CUBIT_RESABS );
    }
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Spatially compare two RefEdges.  Does not go down to EDGE
//                 level.  It does it at the ref_edge level so the parameter
//                 values are consistant.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 04/07/97
//-------------------------------------------------------------------------
CubitBoolean RefEdge::about_spatially_equal(
    RefEdge* ref_edge_ptr_2,
    double tol_factor,
    CubitSense* sensePtr,
    CubitBoolean notify_refEntity )

{
  if( this == ref_edge_ptr_2)
  {
    if (sensePtr)
      *sensePtr = CUBIT_FORWARD;
    if (notify_refEntity)
      remove_compare_data();
    return CUBIT_TRUE;
  }
  
  CubitBoolean spatially_equal = CUBIT_FALSE;
  CubitSense rel_sense = CUBIT_FORWARD;
  CubitStatus stat = CUBIT_SUCCESS;
  stat = relative_sense( ref_edge_ptr_2, tol_factor,
                         &rel_sense, spatially_equal);

  if (stat != CUBIT_SUCCESS || !spatially_equal)
  {
    return CUBIT_FALSE;
  }

  if (sensePtr)
    *sensePtr = rel_sense;

    //compare the start and end vertices to be spatially equal.
  RefVertex* this_start = start_vertex();
  RefVertex* this_end = end_vertex();
  RefVertex* edge2_start = ref_edge_ptr_2->start_vertex();
  RefVertex* edge2_end = ref_edge_ptr_2->end_vertex();
  
    // Swap vertices to simplify things later.
  if (rel_sense == CUBIT_REVERSED)
    std::swap(edge2_start, edge2_end);

    //compare vertex locations unless force_merge is true.
      // closed curve case
    if (this_start == this_end || edge2_start == edge2_end)
    {
      if ((this_start != this_end)   ||
          (edge2_start != edge2_end) ||
          !this_start->about_spatially_equal(edge2_start, tol_factor, CUBIT_FALSE))
        return CUBIT_FALSE;
    }
    else
    {
      if ((this_start == edge2_end) ||
          (this_end == edge2_start) ||
          !this_start->about_spatially_equal(edge2_start, tol_factor, CUBIT_FALSE) ||
          !this_end->about_spatially_equal(edge2_end, tol_factor, CUBIT_FALSE))
        return CUBIT_FALSE;
    }

    //Now if they match report it.
    //Do vertices explicitly here rather than in
    //RefVertex::about_spatially_equal(..) because we don't call 
    //RefVertex::about_spatially_equal(..) if force_merge is true.
  if (notify_refEntity)
  {
    this->notify(ref_edge_ptr_2, COMPARISON_FOUND);
    if (this_start != edge2_start)
      this_start->notify(edge2_start, COMPARISON_FOUND);
    else
      this_start->remove_compare_data();
    if (this_end != edge2_end)
      this_end->notify(edge2_end, COMPARISON_FOUND);
    else
      this_end->remove_compare_data();
  }
  
  return CUBIT_TRUE;
}

int RefEdge::validate()
{
    //- This function determines whether the entity is valid.
    //- Several types of checks can be done, 
  int error = 0;
  
    // Perform general RefEntity checks (measure > 0)
  error += RefEntity::validate();
  
    // Pass through to curve and add in its validation
  Curve *curve = get_curve_ptr();
  
    // check curve ptr
  if (curve != NULL) {
      // Check underlying curve
    DLIList <TopologyEntity*> bad_entities;
    error += curve->validate(entity_name(),bad_entities);
  } else {
    PRINT_WARNING("\tWARNING: Null underlying curve for %s, (%s %d)\n",
                  entity_name().c_str(), class_name(), id());
    error++;
  }
  return error;
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********

//-------------------------------------------------------------------------
// Purpose       : Initializes the member data of the RefEdge
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/07/96
//-------------------------------------------------------------------------
void RefEdge::initialize()
{
    // Initialize some member data
  refEdgeClone = 0;
  markedFlag = CUBIT_FALSE;

    // Make sure the arc length is not zero if there are start and end 
    // RefVertex'es already assigned to this RefEdge
  if ( get_arc_length() < CUBIT_DBL_MIN )  
  {
    if ( start_vertex() && end_vertex() )
    {
      CubitVector start_pt = start_vertex()->coordinates();
      CubitVector end_pt   = end_vertex()->coordinates();
      PRINT_WARNING (
          "WARNING (RefEdge::initialize): Edge has zero arclength.\n"
          "  Start vertex location is (%9.2f, %9.2f, %9.2f ).\n"
          "  End   vertex location is (%9.2f, %9.2f, %9.2f ).\n",
          start_pt.x(), start_pt.y(), start_pt.z(),
          end_pt.x(), end_pt.y(), end_pt.z() );
    }
    else if (!mSuppressEdgeLengthWarning)
    {
      PRINT_WARNING(
          "WARNING: Edge found with zero arclength\n"
          "  For cones, this may be normal.\n");
    }
    
  }
  
    // Set the Entity ID for this new RefEdge
   GeometryEntity* geom_ptr = get_geometry_entity_ptr();
   int saved_id = geom_ptr->get_saved_id();
   if ( !saved_id || RefEntityFactory::instance()->get_ref_edge(saved_id) )
   {
     saved_id =  RefEntityFactory::instance()->next_ref_edge_id();
     geom_ptr->set_saved_id(saved_id);
   }
   entityId = saved_id;
  
     // read and initialize attributes
   auto_read_cubit_attrib();
   auto_actuate_cubit_attrib();

     // Assign a default entity name
   assign_default_name();
}

//-------------------------------------------------------------------------
// Purpose       : Return a pointer to the start RefVertex of this RefEdge.
//
// Special Notes : The assumption is that there is only 1 Chain associated
//                 with each RefEdge. Also, the RefVertex associated with
//                 the first CoVertex in this Chain is the start RefVertex
//                 of this RefEdge and the RefVertex associated with the
//                 last CoVertex in this Chain is the end RefVertex of this
//                 RefEdge.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/21/96
//-------------------------------------------------------------------------
RefVertex* RefEdge::start_vertex() 
{
    // Get the first (and only) Chain associated with this RefEdge.
  Chain* chain_ptr = this->get_chain_ptr();
  
    // Ask the Chain for its first RefVertex
  if (chain_ptr == NULL)
  {
    return NULL;
  }
  
  else
  {
    return chain_ptr->start_vertex();
  }
}

//-------------------------------------------------------------------------
// Purpose       : Return a pointer to the end RefVertex of this RefEdge.
//
// Special Notes : The assumption is that there is only 1 Chain associated
//                 with each RefEdge. Also, the RefVertex associated with
//                 the first CoVertex in this Chain is the start RefVertex
//                 of this RefEdge and the RefVertex associated with the
//                 last CoVertex in this Chain is the end RefVertex of this
//                 RefEdge.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/21/96
//-------------------------------------------------------------------------
RefVertex* RefEdge::end_vertex() 
{
    // Get the first (and only) Chain associated with this RefEdge.
  Chain* chain_ptr = get_chain_ptr();
  
    // Ask the Chain for its end (last) RefVertex
  if (chain_ptr == NULL)
  {
    return NULL;
  }
  
  else
  {
    return chain_ptr->end_vertex();
  }
}

void RefEdge::reverse_topology() 
{
  Chain* chain_ptr = get_chain_ptr();
  chain_ptr->reverse_direction();

    // switch co_edge senses
  DLIList<SenseEntity*> co_edge_list;
  get_sense_entity_list(co_edge_list);
  for ( int i = co_edge_list.size(); i--; ) 
    co_edge_list.get_and_step()->reverse_sense();
}

CubitSense RefEdge::sense( RefFace *face )
{
  DLIList<CoEdge*> co_edge_list;
  get_co_edges(co_edge_list, face);
  CoEdge *coedge;
  CubitSense my_sense = CUBIT_UNKNOWN;
  for ( int i = co_edge_list.size(); i--; ) {
    coedge = co_edge_list.get_and_step();
    if(coedge->get_sense() == CUBIT_FORWARD)
      {
	if (my_sense == CUBIT_REVERSED)
	  return CUBIT_UNKNOWN;
	my_sense = CUBIT_FORWARD;
      }
    else
      {
	if (my_sense == CUBIT_FORWARD)
	  return CUBIT_UNKNOWN;
	my_sense = CUBIT_REVERSED;
      }
  }
  return my_sense;
}

CubitBoolean RefEdge::is_tolerant()
{
   Curve* curve_ptr = get_curve_ptr();
   if (curve_ptr == NULL) {
      PRINT_WARNING("\tWARNING: Null underlying curve for %s, (%s %d)\n",
         entity_name().c_str(), class_name(), id());
      return CUBIT_FALSE;
   }

   return curve_ptr->is_tolerant();
}

CubitVector RefEdge::curve_center()
{
  CubitVector p1 = start_vertex()->coordinates();
  CubitVector p2 = end_vertex()->coordinates();
  if ( start_vertex() == end_vertex() )
    {
      mid_point(p2);
    }
  p1 += p2;
  p1 /= 2.0;
  return p1;
}

CubitStatus RefEdge::get_graphics( GMem& polyline, double tolerance )
{
  Curve* curve_ptr = get_curve_ptr();
  if (!curve_ptr)
  {
    PRINT_ERROR("RefEdge %d is invalid -- no attached Curve.\n",id());
    return CUBIT_FAILURE;
  }
  
  int junk;
  return curve_ptr->get_geometry_query_engine()->get_graphics(curve_ptr,junk,
    &polyline, tolerance);
}

//-------------------------------------------------------------------------
// Purpose       : Reverse RefEdge sense relative to Curve(s)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 10/14/03
//-------------------------------------------------------------------------
void RefEdge::reverse_tangent()
{
  bridge_manager()->reverse_bridge_senses();
  reverse_topology();
}

void RefEdge::suppress_edge_length_warning(bool flag)
{
	mSuppressEdgeLengthWarning = flag;
}
    
