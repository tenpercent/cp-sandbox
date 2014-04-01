//-------------------------------------------------------------------------
// Filename      : RefFace.cpp
//
// Purpose       : This file contains the implementation of the class 
//                 RefFace. 
//
// Special Notes :
//
// Creator       : 
//
// Creation Date : 
//
// Owner         : 
//-------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include "CubitDefines.h"
#include "CubitVector.hpp"

#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include "CubitObserver.hpp"
#include "RefEntityFactory.hpp"
#include "GeometryQueryTool.hpp"
#include "GfxDebug.hpp"
#include "GeometryDefines.h"
#include "Loop.hpp"

#include "CoFace.hpp"
#include "CoEdge.hpp"

#include "Surface.hpp"

// lists
#include "DLIList.hpp"

#include "database.hpp"

#include "CastTo.hpp"
#include "CubitString.hpp"
#include "CubitUtil.hpp"
//for measuring/cacheing the area of the surface.
#include "GeomMeasureTool.hpp"

#include "GeometryUtil.hpp"

//static RefEdge* find_edge_to_adjust(DLIList<RefEdge*>& ref_edge_list);

//-------------------------------------------------------------------------
// Purpose       : Constructor with a pointer to a Surface 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------
RefFace::RefFace(Surface* surfacePtr)
{
     // Set the GeometryEntity pointer   if (surfacePtr != NULL)
   if (surfacePtr != NULL)
   {
      set_geometry_entity_ptr(surfacePtr) ;
   }
   else
   {
      PRINT_ERROR("In the RefFace(Surface*) constructor\n"
                  "       Input Surface pointer is NULL\n");
      assert(CUBIT_FALSE);
   }
   
     // Initialize the member data
   initialize();
}

//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes : Note that the GeometryEntity associated with this
//                 RefEntity is deleted in the destructor of the
//                 BasicTopologyEntity class.
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------
RefFace::~RefFace()
{
    // Delete the hardpoints associated with this RefFace
  
  
//     // Delete the contents of the SurfVertexList
//   for ( i = surfVertexList.size(); i > 0; i --)
//   {
//     delete surfVertexList.get_and_step();
//   }
  
    // Delete the contents of the HardPointList (these are RefVertices).
    // Calling the remove() function ensures that the DAG cleaned out 
    // appropriately.
  int i;
  for ( i = hardPointList.size(); i > 0; i --)
  {
    hardPointList.get_and_step()->remove_from_DAG();
  }
  
  remove_from_observers();	
}

double RefFace::get_crack_length()
{          
  double length = 0.0;
  if( is_parametric() )
  {
    double u_period_val, v_period_val;
//These just give the total period for the surface, I'm assuming that
//it will usually be 2PI
    double lower, upper;
    CubitBoolean u_period = is_periodic_in_U( u_period_val );
    CubitBoolean v_period = is_periodic_in_V( v_period_val );
      //The Booleans tell us if (duh) it's periodic in that direction
    if( ( u_period ) && ( v_period ) )//Uh-oh, periodic in both
//directions, now we have to do some magic to find out the crack
//length.
    {
        // find whether loops differ more in u or v direction.
        // the direction that they differ the most in is the direction
        // that we will step along to determine crack length.
      int m;
      double lowerv[2], upperv[2], loweru[2], upperu[2];
        //the eight variables above keep track of the lowest and
        //highest u- and v-parameter values for both boundaries of the
        //surface
      
      lowerv[0] = lowerv[1] = CUBIT_DBL_MAX;
      upperv[0] = upperv[1] = CUBIT_DBL_MAX;
      loweru[0] = loweru[1] = CUBIT_DBL_MAX;
      upperu[0] = upperu[1] = CUBIT_DBL_MAX;
      DLIList<Loop*> loops;
      this->loops( loops );
        // should be exactly 2 loops
      for( m = 0; m < loops.size() && m < 2; m++ )
      {
        int n;
        Loop* func_loop = loops.get_and_step();
        DLIList<RefEdge*> edges;
        func_loop->ordered_ref_edges( edges );
        for( n = edges.size(); n > 0; n-- )
        {
            // sample edges at 8 locations to approximate maximum and minimum
            // parameter extents of the line
          RefEdge* func_edge = edges.get_and_step();
          const int num_segments = 8;
          CubitVector location[num_segments];
          int i;
          for ( i = num_segments; i--; )
          {
            func_edge->position_from_fraction( ((double) i) / ((double) num_segments),
                                               location[i] );
            
          }
          CubitVector* closest_location = NULL;
          double temp_u, temp_v;

          for ( i = num_segments; i--; )
          {
              //now we fill temp_u and temp_v with the u and v values
              //of the closest point on the surface to the location[i]
              //vector--Note--the closest-point should be right where
              //the vector is, since that's where we got the vector
              //from
            this->u_v_from_position( location[i], temp_u, temp_v, closest_location );
            if( upperv[m] == CUBIT_DBL_MAX )
               upperv[m] = temp_v;
            if( lowerv[m] == CUBIT_DBL_MAX )
               lowerv[m] = temp_v;
            
            if( temp_v > upperv[m] )
               upperv[m] = temp_v;
            if( temp_v < lowerv[m] )
               lowerv[m] = temp_v;
            
            if( upperu[m] == CUBIT_DBL_MAX )
               upperu[m] = temp_u;
            if( loweru[m] == CUBIT_DBL_MAX )
               loweru[m] = temp_u;
            
            if( temp_u > upperu[m] )
               upperu[m] = temp_u;
            if( temp_u < loweru[m] )
               loweru[m] = temp_u;
          }
        }
      }
      if(DEBUG_FLAG(99))
      {
        int debug_int;
        
        for(debug_int = 0; debug_int <2; debug_int ++)
        {
          PRINT_INFO("Loweru[%d] = %f Upperu[%d] = %f\n",
                     debug_int,loweru[debug_int],debug_int,upperu[debug_int]);
        }
        
        for(debug_int = 0; debug_int <2; debug_int ++)
        {
          PRINT_INFO("Lowerv[%d] = %f Upperv[%d] = %f\n",
                     debug_int,lowerv[debug_int],debug_int,upperv[debug_int]);
        }
      }
 
//Okay, we look at our u-params.  The upperu of one loop should be
//lower than the loweru of the other loop.  If they cross, then we're
//most likely looking at the direction parallel to the loop, and the
//lower[1] and lower[0] values are probably very close to each other
      double u_start, u_end, v_start, v_end;
      u_start = u_end = v_start = v_end = 0.;
//u_diff is the gap between the loops.
        //u_start is the midpoint of the lower loop's u_space
        //u_end is the midpoint of the upper loop's u_space
      double u_diff = 0.;
                  
      if ( upperu[1] < loweru[0] )
      {
        u_diff = loweru[0] - upperu[1];
        u_start = (loweru[1] + upperu[1]) / 2.;
        u_end   = (loweru[0] + upperu[0]) / 2.;
      }
      else if ( upperu[0] < loweru[1] )
      {
        u_diff = loweru[1] - upperu[0];
        u_start = (loweru[0] + upperu[0]) / 2.;
        u_end   = (loweru[1] + upperu[1]) / 2.;
      }

        //Okay, we have u_start and u_end.  We have to decide how to
        //get from one to the other.  This comes up because u_start
        //might be smaller than u_end, but the trimmed surface might
        //not exist in between.  If this is the case, then we need to
        //go around the other way to get to u_end
      double low_u, mid_u, high_u;
      double low_v, mid_v, high_v;
      
      mid_u = (u_start + u_end )/2.0;
      
      get_param_range_V( low_v, high_v );
      PRINT_DEBUG_99("low_v = %f high_v = %f\n",low_v, high_v);
      PRINT_DEBUG_99("u_start = %f u_end = %f\n",u_start, u_end);
      mid_v = (low_v + high_v) / 2.;
      double delta_u = fabs(u_end - u_start);
            
      switch(point_containment(mid_u, mid_v))
      {
        case CUBIT_PNT_OUTSIDE:
           PRINT_DEBUG_99("Point outside\n");
           delta_u = -(u_period_val - delta_u);
           break;
        case CUBIT_PNT_INSIDE:
           PRINT_DEBUG_99("Point inside\n");
           break;
        case CUBIT_PNT_BOUNDARY:
           PRINT_DEBUG_99("Point on boundary\n");
           break;
        case CUBIT_PNT_UNKNOWN:
        default:
           PRINT_DEBUG_99("Point Unknown\n");
           break;
      }
      
      delta_u = delta_u / 10.;
         

      double v_diff = 0.;
      if ( upperv[1] < lowerv[0] )
      {
        v_diff = lowerv[0] - upperv[1];
        v_start = (lowerv[1] + upperv[1]) / 2.;
        v_end   = (lowerv[0] + upperv[0]) / 2.;
      }
      
      else if ( upperv[0] < lowerv[1] )
      {
        v_diff = lowerv[1] - upperv[0];
        v_start = (lowerv[0] + upperv[0]) / 2.;
        v_end   = (lowerv[1] + upperv[1]) / 2.;
      }

      mid_v = (v_start + v_end )/2.0;
      
      get_param_range_U( low_u, high_u );
      PRINT_DEBUG_99("low_u = %f high_u = %f\n",low_u, high_u);
      PRINT_DEBUG_99("v_start = %f v_end = %f\n",v_start, v_end);
      mid_u = (low_u + high_u) / 2.;
      double delta_v = fabs(v_end - v_start);
            
      switch(point_containment(mid_u, mid_v))
      {
        case CUBIT_PNT_OUTSIDE:
           PRINT_DEBUG_99("Point outside\n");
           delta_v = -(v_period_val - delta_v);
           break;
        case CUBIT_PNT_INSIDE:
           PRINT_DEBUG_99("Point inside\n");
           break;
        case CUBIT_PNT_BOUNDARY:
           PRINT_DEBUG_99("Point on boundary\n");
           break;
        case CUBIT_PNT_UNKNOWN:
        default:
           PRINT_DEBUG_99("Point Unknown\n");
           break;
      }
      delta_v = delta_v / 10.;
      
        //Okay, by now we have the delta_u and delta_v values.
        //Theoretically they are the distance between the two
        //parameters' midpoints, divided by 10

        //u_diff and v_diff are the gap between the lowest point on
        //the high one and the highest point on the low loop, they
        //will always be positive or zero

        //Note, if a loop is a complete circle, then the parameter
        //value going around that way will yield a delta_? and ?_diff
        //of zero
      int i;
      CubitVector start1, end1;
      CubitVector start2, end2;
      double length1 = 0, length2 = 0;
      PRINT_DEBUG_99("delta_v %e, delta_u %e\n",delta_v, delta_u); //zzyk
      PRINT_DEBUG_99("v_diff %e, u_diff %e\n",v_diff, u_diff); //zzyk

        //We want to get the longest crack length in the surface, so
        //if v_diff is smaller than u_diff, we'll step around the
        //u_direction, and vice versa
      if( v_diff < u_diff )
      { 
        get_param_range_V( low_v, high_v );
        PRINT_DEBUG_99("low_v = %f high_v = %f\n",low_v, high_v);
        PRINT_DEBUG_99("u_start = %f u_end = %f\n",u_start, u_end);
        mid_v = (low_v + high_v) / 2.;

        for( i = 0; i < 10; i++ )
        {
          start1 = this->position_from_u_v( (u_start+(delta_u*i)), low_v );
          end1 = this->position_from_u_v( (u_start+(delta_u*(i+1))), low_v );
          if(DEBUG_FLAG(99))
          {
            GfxDebug::draw_vector(start1,end1, CUBIT_RED);
            GfxDebug::flush();
            PRINT_INFO("At position %f in u\n",u_start + (delta_u*i));
          }
          length1 += start1.distance_between( end1 );
        }
        for( i = 0; i < 10; i++ )
        {
          start2 = this->position_from_u_v( (u_start+(delta_u*i)), mid_v );
          end2 = this->position_from_u_v( (u_start+(delta_u*(i+1))), mid_v );
          if(DEBUG_FLAG(99))
          {
            GfxDebug::draw_vector(start2,end2, CUBIT_BLUE);
            GfxDebug::flush();
            PRINT_INFO("At position %f in u\n",u_start + (delta_u*i));
          }
          length2 += start2.distance_between( end2 );
        }
      }
      else //u_diff is less than v_diff
      {
        get_param_range_U( low_u, high_u );
        PRINT_DEBUG_99("low_u = %f high_u = %f\n",low_u, high_u);
        PRINT_DEBUG_99("v_start = %f v_end = %f\n",v_start, v_end);
          //Remember, v_start is halfway between the low and high
          //points of the lower loop, v_end is halfway between the low
          //and high points of the upper loop
        mid_u = (low_u + high_u) / 2.;
        mid_v = (v_start + v_end )/2.0;
          //Now we step aroung the loop of edges in ten places, and
          //add up the chord lengths of those 10 chords to get our
          //crack length.  We do it twice and average just to make
          //sure everything is kosher.
        for( i = 0; i < 10; i++ )
        {
          start1 = this->position_from_u_v( low_u, (v_start+(delta_v*i)) );
          end1 = this->position_from_u_v( low_u, (v_start+(delta_v*(i+1))) );
          if(DEBUG_FLAG(99))
          {
            GfxDebug::draw_vector(start1,end1, CUBIT_RED);
            GfxDebug::flush();
            PRINT_INFO("At position %f in v\n",v_start + (delta_v*i));
          }
          length1 += start1.distance_between( end1 );
        }
        for( i = 0; i < 10; i++ )
        {
          start2 = this->position_from_u_v( mid_u, (v_start+(delta_v*i)) );
          end2 = this->position_from_u_v( mid_u, (v_start+(delta_v*(i+1))) );
          if(DEBUG_FLAG(99))
          {
            GfxDebug::draw_vector(start2,end2, CUBIT_BLUE);
            GfxDebug::flush();
            PRINT_INFO("At position %f in v\n",v_start + (delta_v*i));
          }
          length2 += start2.distance_between( end2 );
        }
      }
      
      length = ( (length1 + length2) / 2.0 );
    }
    else if( v_period )
    {
      int i;
      get_param_range_U( lower, upper );
      double low_v, high_v, mid_v;
      get_param_range_V( low_v, high_v );
      mid_v = (low_v + high_v) / 2.;
      double delta = ( ( upper - lower ) / 10.0 );
      CubitVector start, end;
      for( i = 0; i < 10; i++ )
      {
        start = this->position_from_u_v( (lower+(delta*i)), mid_v );
        end = this->position_from_u_v( (lower+(delta*(i+1))), mid_v );
      //  GfxDebug::draw_vector(start, end, CUBIT_RED);
      //  GfxDebug::flush();
        length += start.distance_between( end );
      }
    }
    else if( u_period )
    {
      int i;
      get_param_range_V( lower, upper );
      double low_u, high_u, mid_u;
      get_param_range_U( low_u, high_u );
      mid_u = (low_u + high_u) / 2.;
      double delta = ( ( upper - lower ) / 10.0 );
      CubitVector start, end;
      for( i = 0; i < 10; i++ )
      {
        start = this->position_from_u_v( mid_u, (lower+(delta*i)) );
        end = this->position_from_u_v( mid_u, (lower+(delta*(i+1))) );
        length += start.distance_between( end );
      }
    }
    else
    {
      PRINT_WARNING( "Surface %d is_periodic, but has neither a "
                     "u_period or v_period.\n", id() );
    }
  }
  if ( length == 0.0 )
  {
    PRINT_DEBUG_99("Crack length from uv was zero. Trying geometric method.\n");
    //This basically means that this function failed.  Lets try another
    //method that may be less generic but will work for map/submap type 
    //surfaces where getting this crack length is more important.
    length = find_crack_length_no_uv();

  }
  PRINT_DEBUG_99("Crack_length is %f\n",length);
  return length;
}
double RefFace::find_crack_length_no_uv()
{
  double length = 0.0;
  //assume that this surface has exactly 2 loops, again this is for
  //map/submap.  This function also assumes the surface is periodic.
  //ie, this basically assumes some sort of cylinder wall, torus wall or
  //something that looks like that.
  if ( this->num_loops() != 2 || !this->is_periodic() )
    return length;

  //Pick a point on one of the loops, find the closest point to that on the other
  //loops.  Then find the closest point to that point on the first loop.  Keep iterating
  //till we have the two closest points between the loops (may not be closest, but they
  //are mutually the closest to each other.).
  DLIList<Loop*> loops;
  this->loops(loops);
  Loop *first_loop = loops.get();
  Loop *other_loop = loops.next();
  RefEdge *tmp_edge = first_loop->co_edge()->get_ref_edge_ptr();
  CubitVector this_point, closest_point;
   // use the mid point of the first curve as the start point.
  if ( tmp_edge->mid_point(this_point) != CUBIT_SUCCESS )
    return length;
  CoEdge *other_co_edge = NULL;
  GeometryUtil *gu = GeometryUtil::instance();
  CoEdge *closest_co_edge = gu->closest_loop_coedge(other_loop, this_point,
                                                 other_co_edge, &closest_point);
  if ( closest_co_edge == NULL )
    return length;
  CubitVector this_closest;
  int counter = 0;
  const int MAX_LOOP_ITR = 6;
  for(;;)
  {
    if ( counter >= MAX_LOOP_ITR )
      break;
    counter++;
    //find a new this_closest from closest_point.
    closest_co_edge = gu->closest_loop_coedge(first_loop, closest_point,
                                           other_co_edge, &this_closest);
    if ( closest_co_edge == NULL )
      return length;
    if ( (this_point-this_closest).length_squared() <= GEOMETRY_RESABS*GEOMETRY_RESABS )
      break;
    //make this_closest our point.
    this_point = this_closest;
    //get a new closest_point.
    closest_co_edge = gu->closest_loop_coedge(other_loop, this_point,
                                           other_co_edge, &closest_point);
    if ( closest_co_edge == NULL )
      return length;
  }
  //if we are here we succeeded.  Our points are, this_point and closest_point.
  //Now, find some linear steps between these points, move them to the surface,
  //and calculate the length.
  int num_steps = 10;
  //Find the num_steps positions by treating this_point and closest_point as the
  //end points of a curve and parameterizing it...
  int ii;
  double param_step = 1/(double)num_steps;
  double curr_param = param_step;
  CubitVector start_point = this_point;
  CubitVector next_point;
  CubitVector diff_vec = closest_point - this_point;
  for ( ii = 0; ii < num_steps; ii++ )
  {
    next_point = this_point + curr_param*(diff_vec);
    this->move_to_surface(next_point);
    length += (next_point-start_point).length();
    start_point = next_point;
    curr_param += param_step;
  }
  PRINT_INFO("Crack length (no u-v) is = %f\n", length);

  return length;
}


CubitVector RefFace::normal_at ( const CubitVector& location, 
                                 RefVolume* volume, 
                                 double* u_guess, double* v_guess)
{
     // The Surface::normal_at function not only returns the normal, but 
     // also returns the point on the actual surface that is closest to 
     // the input location.  The normal is actually computed at *this*
     // (closest) point on the surface, not at the input location.
   CubitVector normal;
   CubitStatus result;
   Surface* surface_ptr = get_surface_ptr();
   
  if(!u_guess || !v_guess) // no guess was provided
  {
    if (u_guess || v_guess)
    {
      PRINT_ERROR("normal_at(): neither or both of u_guess and "
                  "v_guess must be specified.\n");
      assert(u_guess && v_guess);
    }
    result = surface_ptr->closest_point(location, NULL, &normal);
  }
  else
  {
    result = surface_ptr->closest_point_uv_guess(location, 
                                                       *u_guess, *v_guess,
                                                       NULL, &normal );
  }

   if (result == CUBIT_FAILURE)
   {
      PRINT_ERROR("In RefFace::normal_at\n"
                  "       Could not compute the requested normal at "
                  "location {%f %f %f} on %s (surface %d).\n",
                  location.x(), location.y(), location.z(),
                  this->entity_name().c_str(),
                  this->id());
      assert ( result == CUBIT_SUCCESS );
      return CubitVector(0.0, 0.0, 0.0);
   }
   
   if (surface_ptr->bridge_sense() == CUBIT_REVERSED)
     normal = -normal;
   
   if ( volume )
   {
      CubitSense s = sense( volume );
      if( s != CUBIT_FORWARD && s != CUBIT_REVERSED )
      {
        if(!(volume->is_sheet()))
        {
          PRINT_ERROR("Surface %d has bad sense information with respect to volume %d\n"
                      "       Probably is a 2-sided surface embedded in volume %d\n"  
                      "       Cannot handle this case.\n", id(), volume->id(), volume->id() );
          return CubitVector(0.0, 0.0, 0.0);
        }
      }
      
      if ( s == CUBIT_REVERSED )
          normal = -normal;
   }
   
   return normal;
}

void RefFace::find_closest_point_trimmed(CubitVector from_point,
                                         CubitVector& point_on_surface)
{
   get_surface_ptr()->closest_point_trimmed(from_point, point_on_surface);
}

void RefFace::move_to_surface ( CubitVector& location, 
                                double* u_guess, double* v_guess )
{
  CubitVector closest_location;
  CubitStatus result;

  if(!u_guess || !v_guess) // no guess was provided
  {
    if (u_guess || v_guess)
    {
      PRINT_ERROR("move_to_surface(): neither or both of u_guess and "
                  "v_guess must be specified.\n");
      assert(u_guess && v_guess);
    }
    result = get_surface_ptr()->closest_point(location, &closest_location);
  }
  else
  {
    result = get_surface_ptr()->closest_point_uv_guess(location, 
                                                       *u_guess, *v_guess,
                                                       &closest_location);
  }
  

  if (result == CUBIT_FAILURE)
  {
    PRINT_ERROR("In RefFace::move_to_surface\n"
                "       Could not compute the closest point to "
                "location {%f %f %f} on %s (surface %d).\n",
                location.x(), location.y(), location.z(),
                this->entity_name().c_str(),
                this->id());
    assert ( result == CUBIT_SUCCESS );
    return;
  }

  location.set ( closest_location.x(), 
                closest_location.y(), 
                closest_location.z() );
}

CubitPointContainment RefFace::point_containment( const CubitVector &point )
{
  Surface *surf = get_surface_ptr();
  return surf->point_containment(point);
}
  
CubitPointContainment RefFace::point_containment( double u, double v )
{
  Surface *surf = get_surface_ptr();
  return surf->point_containment(u, v);
}
/*
CubitPointContainment RefFace::point_containment( CubitVector &point, 
                                                  double u, double v )
{
  Surface *surf = get_surface_ptr();
  return surf->point_containment(point, u, v);
}
*/
CubitStatus RefFace::get_principal_curvatures( const CubitVector& point,
                                               double& curvature1,
                                               double& curvature2,
                                               RefVolume* ref_volume_ptr )
{
   Surface* surface_ptr = get_surface_ptr();
   
     // Call the relevant function to compute the curvatures
   CubitStatus status = surface_ptr->
       principal_curvatures( point, curvature1, curvature2 );
   if ( status != CUBIT_SUCCESS )
       return status;
       
   if (surface_ptr->bridge_sense() == CUBIT_REVERSED) {
     curvature1 = -curvature1;
     curvature2 = -curvature2;
   }
       
   if ( ref_volume_ptr ) {
      CubitSense s = sense( ref_volume_ptr );
      assert( s == CUBIT_FORWARD || s == CUBIT_REVERSED );
      if (  s == CUBIT_REVERSED ) {
         curvature1 = -curvature1;
         curvature2 = -curvature2;
      }
   }
   return CUBIT_SUCCESS;
}

int RefFace::genus() 
{
  int nloops = num_loops();
  if (nloops > 0) return (nloops - 1);
  else {
      // need to compute poles
    int num_poles = 0;
    double u_low, u_high, v_low, v_high;
    get_param_range_U(u_low, u_high);
    get_param_range_V(v_low, v_high);
    if (is_singular_in_U(u_low)) num_poles++;
    if (is_singular_in_U(u_high)) num_poles++;
    if (is_singular_in_V(v_low)) num_poles++;
    if (is_singular_in_V(v_high)) num_poles++;
    return -(num_poles+1);
  }
}

CubitVector RefFace::center_point ()
{
   CubitVector center_pt = bounding_box().center();
   move_to_surface(center_pt);
   return center_pt;
}

int RefFace::number_of_Loops ()
{
   int number_of_loops = 0;
   
     // Get the GroupingEntities (Loops) associated with this 
     // BasicTopologyEntity (RefFace)
   DLIList<GroupingEntity*> loopList;
   if ( this->get_grouping_entity_list(loopList) == CUBIT_SUCCESS)
   {
      number_of_loops = loopList.size();
   }
   
   else
   {
      PRINT_ERROR("In RefFace::number_of_Loops\n"
                  "       Unknown problem retrieving Loops "
                  "for %s (surface %d).\n",
                  entity_name().c_str(),
                  this->id());
      number_of_loops = 0;
   }
   
   return number_of_loops;
}

CubitSense RefFace::sense(RefVolume* volume)
{
   SenseEntity* co_face = find_sense_entity(volume);
   return co_face ? co_face->get_sense() : CUBIT_UNKNOWN;
}

CubitSense RefFace::sense( RefFace* face_ptr ) 
{
	DLIList<RefEdge*> edge_list, other_edge_list;
	ref_edges( edge_list );
	face_ptr->ref_edges( other_edge_list );
	edge_list.intersect( other_edge_list );
	
	CubitSense result = CUBIT_UNKNOWN;
	if( edge_list.size() > 0 )
	{
		RefEdge* edge = edge_list.get_and_step();
		if( edge->sense(this) == edge->sense(face_ptr) )
			result = CUBIT_REVERSED;
		else result = CUBIT_FORWARD;
	}
	for( int i = edge_list.size(); i > 1; i-- )
	{
		RefEdge* edge = edge_list.get_and_step();
		CubitSense temp = (edge->sense(this)==edge->sense(face_ptr)) 
			? CUBIT_FORWARD : CUBIT_REVERSED;
		if( temp != result )
		{
			result = CUBIT_UNKNOWN;
			break;
		}
	}
	return result;
}
			


//-------------------------------------------------------------------------
// Purpose       : Spatially compare two RefFaces.  Compare bounding boxes
//                 first, then compare each of the ref-edges.  This function
//                 works strictly off of the ref-entities. 
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 04/07/97
//-------------------------------------------------------------------------
CubitBoolean RefFace::about_spatially_equal( 
    RefFace* ref_face_ptr_2,
    double tolerance_factor,
    CubitBoolean notify_refEntity,
    CubitBoolean test_bbox,
    int test_internal )
{
     // Get rid of the trivial case...
   if( this == ref_face_ptr_2)
   {
      if (notify_refEntity)
        remove_compare_data();
      return CUBIT_TRUE;
   }

   const double tolerance = tolerance_factor * GEOMETRY_RESABS;
   CubitBox box_1 = this->bounding_box();
   CubitBox box_2 = ref_face_ptr_2->bounding_box();
   if (!box_1.overlap( tolerance, box_2) )
    return CUBIT_FALSE;


   GeometryQueryTool* gqt = GeometryQueryTool::instance();
   DLIList<RefEdge*> ref_edge_list_1, ref_edge_list_2;
   this->ref_edges( ref_edge_list_1 );
   ref_face_ptr_2->ref_edges( ref_edge_list_2 );
      
        //compare the size of the two lists.
   if ( ref_edge_list_1.size() != ref_edge_list_2.size() )
     return CUBIT_FALSE;
   
   if (test_internal == 2) // Do internal test for splines only
   {
     const GeometryType this_type = this->geometry_type();
     const GeometryType othr_type = ref_face_ptr_2->geometry_type();
     if (this_type != SPLINE_SURFACE_TYPE    &&
         this_type != BEST_FIT_SURFACE_TYPE  &&
         this_type != UNDEFINED_SURFACE_TYPE &&
         othr_type != SPLINE_SURFACE_TYPE    &&
         othr_type != BEST_FIT_SURFACE_TYPE  &&
         othr_type != UNDEFINED_SURFACE_TYPE   )
      test_internal = 0;
     else
      test_bbox = CUBIT_FALSE;
   }
      
   
     //This compare precedure does the following :
     //   1. Test the bounding boxes of the 2 faces for equality;
     //      If they are "equal" (within "resabs*tolerance_factor"):
     //   2. Compare the ref-edges.
     //   3. Test a point on the two surfaces. 
     //   4. When notify_refEntity is CUBIT_TRUE, whenever find an
     //      two ReEntity's are spatially equal, notify the RefEntity.
     //**** Reorderd by J.Kraftcheck, Sept 22, 2003 ****
     // - Check mergable curves first.  Then check boxes and finally
     //   the internal position.
     
   DLIList<Loop*> loop_list_1, loop_list_2;
   DLIList<CoEdge*> loop_1_coedges;
   this->loops( loop_list_1 );
   ref_face_ptr_2->loops( loop_list_2 );
   if( loop_list_1.size() != loop_list_2.size() )
     return CUBIT_FALSE;

   CubitSense relative_sense = compare_alignment( ref_face_ptr_2 );

   // match each loop in loop_list_1 with one in loop_list_2
   for( int i1 = loop_list_1.size(); i1 > 0; i1-- )
   {
     Loop* loop_1 = loop_list_1.get_and_step();
     loop_1_coedges.clean_out();
     loop_1->ordered_co_edges( loop_1_coedges );
     bool loop_match = false;

     // check every loop in loop_list_2 to see if it matches
     // the current loop from loop_list_1
     for( int i2 = loop_list_2.size(); (i2 > 0) && !loop_match; i2-- )
     {
       Loop* loop_2 = loop_list_2.step_and_get();
       loop_match = loop_2->about_spatially_equal( loop_1_coedges,
                                                   relative_sense,
                                                   tolerance_factor,
                                                   notify_refEntity );
     }
                                          
     // loop from loop_list_1 did not match any loop in loop_list_2
     if( ! loop_match )
       return CUBIT_FALSE;

     // found a match for the current one, so remove it
     loop_list_2.extract();

   } // for( loop_list_1 )


   if ( test_bbox )
   {
       // This test checks to see that the min and max vectors of the
       // bounding boxes are within 10% of the length of the bbox diagonal.
       // Note that this assumes the default values of resabs=1e-6 and
       // tolerance_factor=500

       // It has already been determined that the RefEdges of the
       // surfaces are mergeable, so the bounding boxes of the 
       // RefEdges should be equivalent.  Consider the bounding box
       // of each RefFace to be the RefFace's box united with the
       // box of all the curves.  This removes any potential issues
       // with non-tight bounding boxes for spline curves from 
       // consideration, while still comparing any extend of the boxes
       // that is the result of some internal feature of the surfaces.
     if (ref_edge_list_1.size())
     {
       int i;
       CubitBox edge_box = ref_edge_list_1.step_and_get()->bounding_box();
       for (i = ref_edge_list_1.size(); i > 1; i--)
        edge_box |= ref_edge_list_1.step_and_get()->bounding_box();
       for (i = ref_edge_list_2.size(); i > 0; i--)
        edge_box |= ref_edge_list_2.step_and_get()->bounding_box();
       box_1 |= edge_box;
       box_2 |= edge_box;
     }

     CubitVector tol_vect(
         CUBIT_MIN(box_1.x_range(), box_2.x_range()),
         CUBIT_MIN(box_1.y_range(), box_2.y_range()),
         CUBIT_MIN(box_1.z_range(), box_2.z_range()) );
     tol_vect *= 200.0 * tolerance;

     if( tol_vect.x() < tolerance ) tol_vect.x(tolerance);
     if( tol_vect.y() < tolerance ) tol_vect.y(tolerance);
     if( tol_vect.z() < tolerance ) tol_vect.z(tolerance);
     
     if( (fabs(box_1.minimum().x() - box_2.minimum().x()) > tol_vect.x()) ||
         (fabs(box_1.maximum().x() - box_2.maximum().x()) > tol_vect.x()) ||
         (fabs(box_1.minimum().y() - box_2.minimum().y()) > tol_vect.y()) ||
         (fabs(box_1.maximum().y() - box_2.maximum().y()) > tol_vect.y()) ||
         (fabs(box_1.minimum().z() - box_2.minimum().z()) > tol_vect.z()) ||
         (fabs(box_1.maximum().z() - box_2.maximum().z()) > tol_vect.z()) )
     {
        return CUBIT_FALSE;
     }
   }

   if ( test_internal != 0 )
   {
       //test a point in the middle.
     CubitVector center_1, center_2;
     CubitVector temp_1 = this->center_point();
     this->find_closest_point_trimmed( temp_1, center_1);
       //Okay, now we have the point.  See if this point is on the other 
       //surface.
     ref_face_ptr_2->find_closest_point_trimmed( center_1, center_2 );
     if ( !gqt->about_spatially_equal(center_1, center_2, tolerance_factor ) )
       return CUBIT_FALSE;
   }


     // If we have come this far, we have found matches for
     // every edge of the FACEs. Now notify the associated RefEntities
     // that a match was found.
   if (notify_refEntity == CUBIT_TRUE )
   {
      this->notify(ref_face_ptr_2,COMPARISON_FOUND);
   }

   return CUBIT_TRUE;

}    

//-NOTE: For this function it is assumed that the second_ref_face
//- is spacially equivalient to the first one.
//- This function could explode if this is not followed.
//- If you don't know about the closness of the two faces, check the
//- previous compare function first.
CubitSense RefFace::compare_alignment( RefFace* second_ref_face_ptr )
{
     //Get the sense by testing the two RefFace's at their common
     //center point.
   CubitVector center_point = this->center_point();
   CubitVector normal_this, normal_second;
   normal_this = this->normal_at( center_point );
   normal_second = second_ref_face_ptr->normal_at( center_point );
   
   double dot = normal_this % normal_second;
   
   CubitSense sense = CUBIT_FORWARD;
   if ( dot < 0 )
   {
      sense = CUBIT_REVERSED;
   }
// Moved this warning into merge tool because this function can
// be used by other code, for which this warning is misleading.
// j.k. - 10/11/01
//   else
//   {
//      PRINT_WARNING("Merging %s (surface %d) and %s (surface %d) "
//                    " which have the same sense.\n"
//                    "This may indicate bad geometry.\n",
//                    entity_name().c_str(), id(),
//                    second_ref_face_ptr->entity_name().c_str(),
//                    second_ref_face_ptr->id() );
//   }
   return sense;
}

class LoopAngles 
{
public:
   Loop* loopPtr;
   double angleMetric;
   
   LoopAngles(Loop *loop_ptr )
    { loopPtr = loop_ptr; }
   
   double angle_metric() 
    { return angleMetric; }
};
//create the sorted dl list for LoopAngles.
#include "SDLList.hpp"
SDLListdeclare(SDLLoopAngles, LoopAngles*, angle_metric, double)

CubitStatus RefFace::ordered_loops( DLIList<Loop*> &loop_list )
{
   CubitStatus status;
   SDLLoopAngles loop_angle_list;
   DLIList<Loop*> temp_loop_list;
     //Get all of the loops for this RefFace.
   loops( temp_loop_list );
   
   if ( temp_loop_list.size() < 2 )
   {
     loop_list += temp_loop_list;
   }
   else
   {
       //Now order them from inside to outside.
     Loop *loop_ptr;
     LoopAngles *loop_angles;
     
     for ( int ii = temp_loop_list.size(); ii > 0; ii-- )
     {
       loop_ptr = temp_loop_list.get_and_step();
       loop_angles = new LoopAngles( loop_ptr );
       status = loop_ptr->get_angle_metric( loop_angles->angleMetric );
       if ( status == CUBIT_FAILURE )
       {
         PRINT_ERROR("In RefFace::ordered_loops\n"
                     "       Unknown problem computing the angle metric"
                     " of the Loop.\n");
         delete loop_angles;
         
         while( loop_angle_list.size() != 0 )
           delete loop_angle_list.remove();
         
         return CUBIT_FAILURE;
       }
      loop_angle_list.append( loop_angles );
     }
     
     loop_angle_list.sort();
     loop_angle_list.reset();
     for ( int jj = loop_angle_list.size(); jj > 0; jj-- )
     {
       loop_list.append( loop_angle_list.get_and_step()->loopPtr );
     }
     
       //Delete loop_angles
     while( loop_angle_list.size() != 0 )
       delete loop_angle_list.remove();
   }
   
   return CUBIT_SUCCESS;
}

int RefFace::co_edge_loops ( DLIList<DLIList<CoEdge*>*>& co_edge_loops )
{
   DLIList<CoEdge*>* co_edge_list;
   DLIList<DLIList<CoEdge*>*> temp_loop_list;
   DLIList<Loop*> loop_list;
   CubitStatus status = CUBIT_FAILURE;
   Loop *loop_ptr;
   
     //Get the ordered loops (outside to inside);
   status = ordered_loops( loop_list );
   if ( status == CUBIT_FAILURE )
       return status;
   
     //Now get the co_edges associated with the loops.
   for ( int ii = loop_list.size(); ii > 0; ii-- )
   {
      loop_ptr = loop_list.get_and_step();
      
        // Get the CoEdges on this Loop (the "24" is just a memory allocation
        // chunking value and doesn't imply that we have a list of size 24!!)
      co_edge_list = new DLIList<CoEdge*>(24);
      
        //Get the ref_edges with respect to the loop.
      status = loop_ptr->ordered_co_edges( *co_edge_list );
      if ( status == CUBIT_FAILURE )
      {
         delete co_edge_list;
         while( temp_loop_list.size() != 0 )
             delete temp_loop_list.remove();
         return status;
      }
      temp_loop_list.append( co_edge_list );
   }
   co_edge_loops += temp_loop_list;
   
   return CUBIT_SUCCESS;
}

int RefFace::ref_edge_loops ( DLIList<DLIList<RefEdge*>*>& ref_edge_loops )
{
// NOTE: all of the ref_edge_list's will need to be deleted by 
//       the calling function...

   DLIList<RefEdge*>* ref_edge_list;
   DLIList<DLIList<RefEdge*>*> temp_loop_list;
   DLIList<Loop*> loop_list;
   CubitStatus status = CUBIT_FAILURE;
   Loop *loop_ptr;
   
     //Get the ordered loops (outside to inside);
   status = ordered_loops( loop_list );
   if ( status == CUBIT_FAILURE )
       return status;
   
     //Now get the ref-edges associated with the loops.
   for ( int ii = loop_list.size(); ii > 0; ii-- )
   {
      loop_ptr = loop_list.get_and_step();
      
        // Get the RefEdges on this Loop (the "24" is just a memory allocation
        // chunking value and doesn't imply that we have a list of size 24!!)
      ref_edge_list = new DLIList<RefEdge*>(24);
      
        //Get the ref_edges with respect to the loop.
      status = loop_ptr->ordered_ref_edges( *ref_edge_list );
      if ( status == CUBIT_FAILURE )
      {
         delete ref_edge_list;
         while( temp_loop_list.size() != 0 )
             delete temp_loop_list.remove();
         return status;
      }
      temp_loop_list.append( ref_edge_list );
   }
   ref_edge_loops += temp_loop_list;
   
   return CUBIT_SUCCESS;
}


RefVolume* RefFace::ref_volume()
{
   DLIList<RefEntity*> entity_list;
   DLIList<RefVolume*> vol_list;
   RefVolume *vol_ptr;
   
     // Get the list of RefVolumes that own this RefFace.
   get_parent_ref_entities( entity_list );
   CAST_LIST( entity_list, vol_list , RefVolume);
   
     // Return the first valid RefVolume from the list.
   vol_list.reset();
   for( int i = vol_list.size(); i>0; i-- )
   {
      vol_ptr = vol_list.get_and_step();
      if( vol_ptr )
          return vol_ptr;
   }
   
     // Print ERROR if no valid RefVolume was found.
   PRINT_ERROR("No RefVolume found for the RefEdge.\n");
   return NULL;
}

//NOTE: There could be more than one CoFace that is associated with this
// volume and ref-faces, (hard-surfaces).
CoFace* RefFace::get_matching_CoFace(RefVolume* ref_volume_ptr)
{
  return dynamic_cast<CoFace*>(find_sense_entity(ref_volume_ptr));
}


void RefFace::add_hard_point( RefVertex* ref_vertex_ptr)
{
   hardPointList.append( ref_vertex_ptr );
}

void RefFace::hard_points(  DLIList<RefVertex*>& new_hard_point_list )
{
   new_hard_point_list = hardPointList;
}

CubitVector RefFace::position_from_u_v (double u, double v)
{
     // Get the Surface this object points to
   Surface* surfacePtr = get_surface_ptr();
   
     // Make sure we get a valid Surface
   assert(surfacePtr != NULL) ;
   
     // Ask the Surface to do the real work
   return surfacePtr->position_from_u_v(u, v) ;
}

CubitStatus RefFace::u_v_from_position (CubitVector const& location,
                                        double& u, 
                                        double& v,
                                        CubitVector* closest_location )
{
     //- This function returns the {u, v} coordinates of the point 
     //- on the Surface closest to the input point (specified in global
     //- space). The closest_location is also returned.
   
     // pass call directly to surface
   Surface *surface = get_surface_ptr();
   
     // check surface ptr
   assert(surface != NULL);
   
   if(is_parametric() == CUBIT_TRUE)
   {
      return surface->u_v_from_position (location, u, v,
                                      closest_location);
   }
   else
   {
     return CUBIT_FAILURE;
   }
}


CubitBoolean RefFace::is_parametric()
{
     //- This function determines whether the underlying geometry of the
     //- Surface is parametrically defined or not.  Returns CUBIT_TRUE if 
     //- it is and CUBIT_FALSE if it is not.
   
     // pass call directly to surface
   Surface *surface = get_surface_ptr();
   
     // check surface ptr
   assert(surface != NULL);
   
   return surface->is_parametric();
}

CubitBoolean RefFace::get_param_range_U( double& lower_bound,
                                         double& upper_bound )
{
     //- Returns the lower and upper parametric bounds of the 
     //- surface in U, if it is parametric.  Otherwise, it returns
     //- CUBIT_FALSE and zeroes for the upper and lower parametric
     //- bounds.
   
     // pass call directly to surface
   Surface *surface = get_surface_ptr();
   
     // check surface ptr
   assert(surface != NULL);
   
   return surface->get_param_range_U(lower_bound, upper_bound);
}


CubitBoolean RefFace::get_param_range_V( double& lower_bound,
                                         double& upper_bound )
{
     //- Returns the lower and upper parametric bounds of the 
     //- surface in V, if it is parametric.  Otherwise, it returns
     //- CUBIT_FALSE and zeroes for the upper and lower parametric
     //- bounds.
   
     // pass call directly to surface
   Surface *surface = get_surface_ptr();
   
     // check surface ptr
   assert(surface != NULL);
   
   return surface->get_param_range_V(lower_bound, upper_bound);
}

CubitBoolean RefFace::is_periodic()
{
     //- This function determines whether the underlying geometry of the
     //- Surface is periodic or not.  Returns CUBIT_TRUE if it is and 
     //- CUBIT_FALSE if it is not.
   
     // pass call directly to surface
   Surface *surface = get_surface_ptr();
   
     // check surface ptr
   assert(surface != NULL);
   
   return surface->is_periodic();
}

CubitBoolean RefFace::is_periodic_in_U( double& period )
{
     //- Determines whether the surface object is 
     //- periodic in the U direction or not.  If it is, it
     //- returns CUBIT_TRUE and the value of the period. Otherwise,
     //- it returns CUBIT_FALSE and a value of 0.0 or the period.
   
     // pass call directly to surface
   Surface *surface = get_surface_ptr();
   
     // check surface ptr
   assert(surface != NULL);
   
   return surface->is_periodic_in_U(period);
}

CubitBoolean RefFace::is_periodic_in_V( double& period )
{
     //- Determines whether the surface object is 
     //- periodic in the V direction or not.  If it is, it
     //- returns CUBIT_TRUE and the value of the period. Otherwise,
     //- it returns CUBIT_FALSE and a value of 0.0 or the period.
   
     // pass call directly to surface
   Surface *surface = get_surface_ptr();
   
     // check surface ptr
   assert(surface != NULL);
   
   return surface->is_periodic_in_V(period);
}
CubitBoolean RefFace::is_singular_in_U( double u_param )
{
     //- Determines whether the surface object is 
     //- singular in the U direction or not.  
     // pass call directly to surface
   Surface *surface = get_surface_ptr();
   
     // check surface ptr
   assert(surface != NULL);
   
   return surface->is_singular_in_U(u_param);
}
CubitBoolean RefFace::is_singular_in_V( double v_param )
{
     //- Determines whether the surface object is 
     //- singular in the V direction or not.  
     // pass call directly to surface
   Surface *surface = get_surface_ptr();
   
     // check surface ptr
   assert(surface != NULL);
   
   return surface->is_singular_in_V(v_param);
}

RefEdge* RefFace::common_ref_edge ( RefFace* input_face_ptr )
{
   DLIList<RefEdge*> this_edge_list;
   ref_edges ( this_edge_list );
   
   for ( int i = this_edge_list.size(); i > 0; i--)
   {
      RefEdge* edge = this_edge_list.get_and_step();
      if (edge->find_sense_entity(input_face_ptr))
        return edge;
   }
   
   return NULL;
}

int RefFace::common_ref_edges ( RefFace* input_face_ptr, DLIList<RefEdge*> &common_edge_list )
{
   DLIList<RefEdge*> this_edge_list;
   ref_edges ( this_edge_list );
   int nedges = 0;
   
   for ( int i = this_edge_list.size(); i > 0; i--)
   {
      RefEdge* edge = this_edge_list.get_and_step();
      if (edge->find_sense_entity(input_face_ptr))
      {
        common_edge_list.append(edge);
        nedges++;
      }
   }
   
   return nedges;
}

RefVolume* RefFace::common_ref_volume ( RefFace* input_face_ptr )
{
   DLIList<RefVolume*> this_volume_list;
   DLIList<RefVolume*> input_volume_list;
   
   ref_volumes ( this_volume_list );
   input_face_ptr->ref_volumes ( input_volume_list );
   
   for ( int i = this_volume_list.size(); i > 0; i--)
   {
      if (input_volume_list.move_to (this_volume_list.get()))
      {
         return this_volume_list.get();
      }
      
      this_volume_list.step();
   }
   
   return NULL;
}

int RefFace::dimension() const
{
   return 2;
}
double RefFace::area()
{
  return GeomMeasureTool::measure_area(this);
}

double RefFace::measure()
{
   return this->area();
}

CubitString RefFace::measure_label()
{
   return "area";
}

//-------------------------------------------------------------------------
// Purpose       : Return a pointer to the surface associated with a face.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
Surface* RefFace::get_surface_ptr() 
{
  return CAST_TO(get_geometry_entity_ptr(), Surface);
}

const Surface* RefFace::get_surface_ptr() const
{
  return CAST_TO(get_geometry_entity_ptr(), Surface);
}


//-------------------------------------------------------------------------
// Purpose       : This function returns CUBIT_TRUE if the underlying 
//                 geometry of the face is planar. CUBIT_FALSE otherwise.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/17/96
//-------------------------------------------------------------------------
CubitBoolean RefFace::is_planar() 
{
     // Cast the generic GeometryEntity pointer to Surface pointer
   Surface* surfacePtr = this->get_surface_ptr() ;
   
     // Check if we have a valid Surface. If so, return the result of
     // querying the Surface if it is planar.
   if ( surfacePtr != NULL )
   {
      GeometryType geo_type;
      geo_type = surfacePtr->geometry_type();
      return geo_type == PLANE_SURFACE_TYPE ? CUBIT_TRUE : CUBIT_FALSE;
   }
   else
   {
      PRINT_WARNING("In RefFace::is_planar\n"
                    "         %s (surface %d) is not associated with a valid\n"
                    "         underlying geoemtric Surface\n",
                    entity_name().c_str(), id()) ;
      return CUBIT_FALSE ;
   }
}

//To get planar surface's origin and normal, both are pure output.
CubitStatus RefFace::get_point_normal( CubitVector& origin, CubitVector& normal )
{
   if( is_planar() == CUBIT_FALSE)
      return CUBIT_FAILURE;

   Surface* surface_ptr = get_surface_ptr();

   if( surface_ptr != NULL )
   {
      if( surface_ptr->get_point_normal( origin, normal ) == CUBIT_FAILURE )
         return CUBIT_FAILURE;
   }
   else 
   {
      PRINT_WARNING("In RefFace::get_point_normal\n"
                    "         %s (surface %d) is not associated with a valid\n"
                    "         underlying geoemtric Surface\n",
                    entity_name().c_str(), id()) ;
      return CUBIT_FAILURE;
   }
   
   if (surface_ptr->bridge_sense() == CUBIT_REVERSED)
     normal = -normal;
   
   return CUBIT_SUCCESS;
}

int RefFace::validate()
{
     //- This function determines whether the entity is valid.
     //- Several types of checks can be done, 
   int error = 0;
   
     // Perform general RefEntity checks (measure > 0)
   error += RefEntity::validate();
   
     // Pass through to surface and add in its validation
   Surface *surface = get_surface_ptr();
   
     // check surface ptr
   if (surface != NULL) {
        // Check underlying surface
	   DLIList <TopologyEntity*> bad_entities;
      error += surface->validate(entity_name(), bad_entities);
   } else {
      PRINT_WARNING("\tWARNING: Null underlying surface for %s, (%s %d)\n",
                    entity_name().c_str(), class_name(), id());
      error++;
   }
   return error;
}

//-------------------------------------------------------------------------
// Purpose       : Initializes all member data
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 09/25/96
//-------------------------------------------------------------------------
void RefFace::initialize()
{
     // Set the Entity ID for this new RefFace
   GeometryEntity* geom_ptr = get_geometry_entity_ptr();
   int saved_id = geom_ptr->get_saved_id();
   if ( !saved_id || RefEntityFactory::instance()->get_ref_face(saved_id) )
   {
     saved_id =  RefEntityFactory::instance()->next_ref_face_id();
     geom_ptr->set_saved_id(saved_id);
   }
   entityId = saved_id;
   
     // Default graphics attributes
   hardPointColor = 1;

     // initialize meshing data

     // initialize the bounding box
   CubitBox bound_box = bounding_box();
   maxPositionDeviation = 
       (bound_box.maximum() - bound_box.minimum()).length()/10.0;

     // read and initialize attributes
   auto_read_cubit_attrib();
   auto_actuate_cubit_attrib();

#ifdef ALPHA_TREADSWEEP
   if(entityId != saved_id)
     geom_ptr->set_saved_id(entityId);
#endif
   
     // Assign a default entity name
   assign_default_name();   
   
}

void RefFace::reverse_topology() 
{
   
   int i;
   
     // switch sense going up in dimension 
   DLIList<CoFace*> co_face_list;
   co_faces( co_face_list );
   for ( i = co_face_list.size(); i--; ) {
      CoFace *co_face = co_face_list.get_and_step();
      co_face->set_sense( CubitUtil::opposite_sense( co_face->get_sense() ) );    
   }
   
     // switch sense going down in dimension 
   DLIList<Loop*> loop_list;
   loops( loop_list );
   for ( i = loop_list.size(); i--; )
       loop_list.get_and_step()->reverse_direction();
}

void RefFace::reverse_normal()
{
   bridge_manager()->reverse_bridge_senses();
   reverse_topology();
}  

CubitBoolean RefFace::set_outward_normal( RefVolume *volume )
{
   CubitSense vol_sense = sense( volume );
   if ( vol_sense == CUBIT_UNKNOWN )
      return CUBIT_FALSE;
   assert( vol_sense == CUBIT_FORWARD || vol_sense == CUBIT_REVERSED );
   if ( vol_sense == CUBIT_REVERSED ) {
      reverse_normal();
      return CUBIT_TRUE;
   }
   // else already right pointing right way.
   return CUBIT_FALSE;
}

CubitStatus RefFace::get_graphics( GMem& facets,
                                   unsigned short normal_tolerance,
                                   double distance_tolerance,
                                   double longest_edge )
{
  int junk1, junk2, junk3;
  Surface* surf_ptr = get_surface_ptr();
  if (!surf_ptr)
  {
    PRINT_ERROR("RefFace %d is invalid -- no attached Surface.\n",id());
    return CUBIT_FAILURE;
  }
  
  return surf_ptr->get_geometry_query_engine()->
    get_graphics(surf_ptr, junk1, junk2, junk3, &facets, 
      normal_tolerance, distance_tolerance, longest_edge );
}

CubitBoolean RefFace::is_closed_in_U()
{
  Surface* surface_ptr = get_surface_ptr();
  return surface_ptr->is_closed_in_U();
}

CubitBoolean RefFace::is_closed_in_V()
{
  Surface* surface_ptr = get_surface_ptr();
  return surface_ptr->is_closed_in_V();
}
