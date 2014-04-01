//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//
// Filename      : RefFace.hpp
//
// Purpose       : A RefFace is a BasicTopologyEntity.  It is a topological
//                 entity that represents a bounded surface.  A RefFace is a
//                 contiguous point set, but can have multiple "holes"
//                 within it. It is a topologically 2-D entity.
//
//                 A RefFace can be meshed.                
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 07/11/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef REFFACE_HPP
#define REFFACE_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN MOTIF INCLUDES         **********
// ********** END MOTIF INCLUDES           **********

// ********** BEGIN OPEN INVENTOR INCLUDES **********
// ********** END OPEN INVENTOR INCLUDES   **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CastTo.hpp"
#include "BasicTopologyEntity.hpp"
#include "LocalStart.h"


// ********** END CUBIT INCLUDES           **********

// ********** BEGIN MACROS DEFINITIONS     **********
// ********** END MACROS DEFINITIONS       **********

// ********** BEGIN TYPEDEF DEFINITIONS    **********
// ********** END TYPEDEF DEFINITIONS      **********

// ********** BEGIN ENUM DECLARATIONS      **********
// ********** END ENUM DECLARATIONS        **********

// ********** BEGIN FORWARD DECLARATIONS   **********


class CoEdge;
template <class X> class DLIList;
class RefVertex;
class CoFace;
class Loop;
class CoFace;
class SurfMeshTool;
class SurfVertexType;
class TDUVSpace;
class Surface;
class GMem;

// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT RefFace : public BasicTopologyEntity
{
public :
  
  friend class RefEntityFactory;
    //- the factory is allowed to call the (private) constructors

    /* constructors/destructors */

  virtual ~RefFace() ;
    //- The destructor

  static const char* get_class_name()
     {
       return "Surface";
     }

  virtual const char* class_name() const
     {
       return get_class_name();
     }
  
    /* topology */
  DagType dag_type() const { return DagType::ref_face_type(); }
  const type_info& entity_type_info() const { return typeid(RefFace); }


  CubitStatus get_co_faces( DLIList<CoFace*> &co_faces_found_list,
                            RefVolume *input_ref_volume_ptr );
    //R CubitStatus
    //R- CubitSuccess/CubitFailure
    //O co_faces_found_list
    //O-Populates the co_faces_found_list with the CoFaces that are
    //-associated with this RefFace.  Because the input RefVolume is sent
    //-in it fills the list just with the CoFaces that are associated with
    //-this RefVolume.  Note that usually there will be only one CoFace with
    //-the RefVolume but with HardSurfaces there may be more (2).

  CubitStatus ordered_loops( DLIList<Loop*> &loop_list );
    //- Gets the loops in order from outside to inside.
    //- This function is used for all the loop extracting from the refface,
    //- i.e., nodes, ref-edges, ref-vertices.  This function orders the
    //- loops based on the angle metric calculation.  I believe this metric
    //- is actually currently calculated in AcisGeometryEngine but I see
    //- no reason that it should be done there.  I believe it could be
    //- done in the Loop class.

  int co_edge_loops ( DLIList<DLIList<CoEdge*>*>& co_edge_loops );
    //- Returns a list of lists.  Each of the included lists contains a list
    //- of CoEdges and represents an ordered list of CoEdges associated
    //- with each of the Loops of this RefFace.
  
  int ref_edge_loops ( DLIList<DLIList<RefEdge*>*>& ref_edge_loops );
    //- Returns a list of lists.  Each of the included lists contains a list
    //- of RefEdges and represents an ordered list of RefEdges associated
    //- with each of the Loops of this RefFace.
    //- NOTE: All of the ref_edge_lists in ref_edge_loops will 
    //-       need to be deleted by the *calling* function.  
  
  void ref_vertex_loops( DLIList<DLIList<RefVertex*>*>& ref_vert_loop_list );
    //- Returns a list of lists.  Each of the included lists contains a list
    //- of RefVertex'es and represents an ordered list of RefVertex'es 
    //- associated with each of the Loops of this RefFace.  
    //- NOTE: All of the ref_vertex_lists in the ref_vert_loop_list will 
    //-       need to be deleted by the *calling* function.
  
  int  number_of_Loops ();
    //- Returns the number of Loops associated with this RefFace
  
  
  RefVolume* ref_volume();
    //- Return the first RefVolume pointer to the volume which owns
    //- this RefFace
    //- Note: There may be more than one RefVolume that owns this RefFace.
    //-       This method just gets the first in the list.

  void hard_points( DLIList<RefVertex*>& new_hard_point_list );
    //- Populate the input DLIList<RefVertex*> with the list of hard points 
    //- that are defined for this RefFace
  
  int adjoins ( RefFace* input_face_ptr );
    //- Returns CUBIT_TRUE if this RefFace adjoins (is connected via a RefEdge)
    //- the input RefFace
  
  RefVolume* common_ref_volume ( RefFace* input_face_ptr );
    //- Returns a common RefVolume* if this RefFace shares one with the
    //- input RefFace

  RefEdge* common_ref_edge ( RefFace* input_face_ptr );
    //- Returns a common RefEdge* if this RefFace shares one with the
    //- input RefFace

  int common_ref_edges ( RefFace* input_face_ptr, DLIList<RefEdge*> &common_edge_list );
    //- Returns all common RefEdges that this face shares with the input face

  CoFace* get_matching_CoFace(RefVolume* ref_volume_ptr) ;
    //R CoFace*
    //R- Returned CoFace pointer
    //I ref_volume_ptr
    //I- The RefVolume to which matching is done.
    //- This function returns the CoFace that is associated with both "this"
    //- RefFace as well as with the input RefVolume.  The function is useful
    //- when a merge operation has resulted in a RefFace that is shared by
    //- more than 1 RefVolume (most often, two). In this case, the RefFace
    //- would be associated with more than one CoFace, each belonging to a 
    //- different RefVolume.
    //- If there is no match (i.e., this RefFace is not associated with the
    //- input RefVolume) then a NULL pointer is returned.

  int genus();
    //- return genus of this surfaces, which is defined as
    //- g = (L ? L-1 : -(P+1)) where L = # loops, P = # poles
  
    /*  geometry */

  virtual CubitVector center_point();
    //- Return the approximate (spatial) center of this RefFace
  
  CubitSense sense(RefVolume* volume);
    //-Determines the sense of "this" with respect to the passed-in volume.
		
	CubitSense sense(RefFace* face_ptr);
		//-Determine the relative sense of the passed face with respect to
		//-this face using the senses of common RefEdges.  i.e. if 
		//-CUBIT_REVERSED is returned, than the sense of the passed face
		//-should be the opposite of that of this face with respect to
		//-any volume.  CUBIT_UNKNOWN is returned if there are no
		//-common RefEdges between RefFaces or there is more than one
		//-common RefEdge.

  CubitSense get_geometry_sense();
    //- Gets the sense of the reface with respect to the underlying
    //- geometry engines representation of the surface.

  CubitBoolean about_spatially_equal ( RefFace* ref_face_ptr_2,
                                       double tolerance_factor = 1.0,
                                       CubitBoolean notify_refEntity =
                                       CUBIT_FALSE,
                                       CubitBoolean test_bbox = CUBIT_TRUE,
                                       int test_internal = 0 );
  
    //R CubitBoolean
    //R-CUBIT_TRUE/CUBIT_FALSE
    //I RefFace*, double, CubitBoolean
    //I- Second RefFace to compare, Tolerance factor to for GEOMETRY_RESABS,
    //I- and flag for notifying compared RefEntities.
    //O CubitBoolean
    //O- If the two RefFaces are spatially equal within the GEOMETRY_RESABS*
    //- the tolerance_factor, then CUBIT_TRUE will be returned.  Otherwise
    //- CUBIT_FALSE is returned.
    //- The comparison is done by first checking the bounding boxes of the
    //- RefFaces.  If this test is passed then the ref_edges of each face
    //- are looped through and compared.  A bounding box check for each
    //- edge is also done first before a comparison, for speed.
  
  CubitSense compare_alignment( RefFace* second_ref_face_ptr );
    //R CubitSense
    //R- Sense of this reface with respect to the second one passed in.
    //I RefFace *
    //I- pointer to second ref face with which the alignment is compared.
    //- This function will compare the sense of the two ref-faces, or
    //- rather their normals.
    //- NOTE: It is ASSUMED that BOTH reffaces are SPATIALLY EQUAL.
    //- If this is not followed this could explode.
  

  CubitVector normal_at(const CubitVector& location, RefVolume* volume=NULL, double* u_guess = NULL, double* v_guess = NULL);
    //- Calculate normal for input location (optional input RefVolume to 
    //- allow for feature consolidation). 
    //- Note that the input location is modified to the coordinates
    //- of the closest point on the surface.
    //-
    //- MJP NOTE:
    //- In the previous implementation, the result of this function call 
    //- would not only be the returned unit vector which is the normal
    //- at the location, but the function would also fill in the myPosition
    //- and myParametricPosition data members of RefFace.  These data
    //- members have been removed in this implementation of RefFace.
    //- However, the Surface::normal_at function that gets called returns
    //- an additional parameter which is the location on the underlying
    //- surface that is closest to the input location.
  void reverse_normal();
    //- switch the sense of this face with respect to the underlying
    //- geometry, so that all normals point in the opposite
    //- direction. The orientation of quads on the surface are switched
    //- to agree with this normal.
  virtual void reverse_topology();

//========  Change Code by DZ of Cat,  on 10/29/98 8:46:59 AM  ========
  CubitBoolean set_outward_normal( RefVolume *volume );
    //- Set the normal of this face to point outward wrt to given volume.
    //- Assumes there is only one coface of the volume for this RefFace.
    //- Uses the above "reverse_normal" function.
    //- return true only if reverse_normal function called.
//========  Change End by DZ of Cat,  on 10/29/98 8:46:59 AM  ========

  virtual void move_to_surface ( CubitVector& location, double* u_guess = NULL, double* v_guess = NULL);
    //- Moves the given node back onto its surface
  
  void find_closest_point_trimmed(CubitVector from_point, 
                                  CubitVector& point_on_surface);
    //R void
    //I CubitVector
    //I- point from which to find closest point on trimmed surface
    //O CubitVector
    //O- point on trimmed surface closest to passed-in point 
    //- This function finds the closest point on a TRIMMED surface to the
    //- passed-in point.  

  CubitPointContainment point_containment( const CubitVector &point );
  CubitPointContainment point_containment( double u, double v );
//  CubitPointContainment point_containment( CubitVector &point, double u, double v );
    //R CubitPointContainment - is the point outside, inside or on the boundary?
    //R- CUBIT_PNT_OUTSIDE, CUBIT_PNT_INSIDE, CUBIT_PNT_BOUNDARY, 
    //   CUBIT_PNT_UNKNOWN
    //I CubitVector
    //I- position to check, known to be on the Surface
    //I double
    //I- u coordinate, if known (significantly faster, if this is known - however
    //                           if not known let the function figure it out)
    //I double
    //I- v coordinate, if known (significantly faster, if this is known - however
    //                           if not known let the function figure it out)
    // NOTE: POINT MUST LIE ON THE SURFACE FOR THIS FUNCTION TO WORK PROPERLY.
  
  CubitStatus get_principal_curvatures( const CubitVector& point,
                                        double& curvature1,
                                        double& curvature2,
                                        RefVolume* ref_volume_ptr = NULL );
    //R CubitStatus 
    //R- CUBIT_SUCCESS/FAILURE
    //I point
    //I- Input location. The coordinates of this input point are 
    //I- modified to those of the point closest to this one, on the
    //I- surface of this RefFace.
    //O curvature1/curvature2
    //O- Output principal curvature values.
    //I ref_volume_ptr
    //I- Input RefVolume pointer
    //- This function first computes the point on the surface closest to the
    //- input point and sets the values of "point" to this closest
    //- location.  Then, the principal curvatures of the surface at this 
    //- new point are computed and returned. If the input RefVolume pointer
    //- is not NULL, it is used when computing the curvatures.
  
  CubitVector position_from_u_v (double u, double v);
    //- Return a CubitVector (representing a position vector corresponding 
    //- to the input point in {u,v} space
  
  CubitStatus u_v_from_position (CubitVector const& location,
                                 double& u, 
                                 double& v,
                                 CubitVector* closest_location = NULL );
    //R CubitStatus 
    //R- CUBIT_SUCCESS/FAILURE
    //I location
    //I- The input point in global space
    //O u, v
    //O- The returned u, v coordinate values (in local parametric space)
    //O- of the closest_point
    //O closest_location
    //O- The point on the Surface closest to the input location
    //I refvolume_ptr
    //I- The reference RefVolume with respect to which, the normal
    //I- is to be computed.  If the pointer is NULL, then the 
    //I- first underlying solid model entity is used to compute the
    //I- normal.
    //- This function returns the {u, v} coordinates of the point 
    //- on the Surface closest to the input point (specified in global
    //- space). The closest_location is also returned.
  
  CubitBoolean is_parametric();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- This function determines whether the underlying geometry of the
    //- Surface is parametrically defined or not.  Returns CUBIT_TRUE if 
    //- it is and CUBIT_FALSE if it is not.
  
  CubitBoolean get_param_range_U( double& lower_bound,
                                  double& upper_bound );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O lower_bound
    //O- The lower bound of the parametric range in the U direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //O upper_bound
    //O- The upper bound of the parametric range in the U direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //- Returns the lower and upper parametric bounds of the 
    //- surface in U, if it is parametric.  Otherwise, it returns
    //- CUBIT_FALSE and zeroes for the upper and lower parametric
    //- bounds.
  
  CubitBoolean get_param_range_V( double& lower_bound,
                                  double& upper_bound );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O lower_bound
    //O- The lower bound of the parametric range in the V direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //O upper_bound
    //O- The upper bound of the parametric range in the V direction.
    //O- This is set to 0.0 if the surface is not parametric.
    //- Returns the lower and upper parametric bounds of the 
    //- surface in V, if it is parametric.  Otherwise, it returns
    //- CUBIT_FALSE and zeroes for the upper and lower parametric
    //- bounds.
  
  CubitBoolean is_periodic();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- This function determines whether the underlying geometry of the
    //- Surface is periodic or not.  Returns CUBIT_TRUE if it is and 
    //- CUBIT_FALSE if it is not.
  
  CubitBoolean is_periodic_in_U( double& period );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the U direction.
    //- Determines whether the surface object is 
    //- periodic in the U direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.
  
  CubitBoolean is_periodic_in_V( double& period );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- The value of the period in the V direction.
    //- Determines whether the surface object is 
    //- periodic in the V direction or not.  If it is, it
    //- returns CUBIT_TRUE and the value of the period. Otherwise,
    //- it returns CUBIT_FALSE and a value of 0.0 or the period.

  CubitBoolean is_singular_in_U( double u_param );
  CubitBoolean is_singular_in_V( double v_param );
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //I double u/v parameter value.
    //- Determines if the surface is singular in a given direction
    //- at a given parameter value.

  CubitBoolean is_closed_in_U();  
  CubitBoolean is_closed_in_V();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- Determines if the surface is closed, smoothly or not in the
    //- given parameter direction.
    //- A periodic surface is always closed but a closed surface is
    //- is not always periodic.

  CubitStatus uv_derivitives( double u_param,
                              double v_param,
                              CubitVector &du,
                              CubitVector &dv );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //O- du, dv
    //- Determines the u and v derivitives from the given parameter
    //- values.

  virtual int dimension() const; 
    //- Returns the geometric dimension of RefFace entities. 

  double area();
    //- get the area of the underlying surface

  virtual double measure();
  virtual CubitString measure_label();
  
  CubitBoolean is_planar();
    //R CubitBoolean
    //R CUBIT_TRUE/CUBIT_FALSE
    //- This function returns CUBIT_TRUE if the underlying geometry
    //- of the face is planar. CUBIT_FALSE otherwise.

  CubitStatus get_point_normal( CubitVector& origin, CubitVector& normal );
    //- Only valid for planar surfaces
    //- Finds the underlying plane's origin and normal (unit) vector
    //- Returns CUBIT_FAILURE if surface is not a plane
  
  virtual int validate();
    //- Check that entity is valid. Returns number of problems detected.

  double get_crack_length();
    //- return the length of the periodic crack, or 0.0 if non-periodic.

    /* geometry modification */

  void add_hard_point(  RefVertex* ref_vertex_ptr );
    //- Add a hard point to this RefFace. 
    //-***************************************************************
    //- MJP Note:
    //- Currently, the new RefVertex that is created is NOT PART OF THE
    //- main DAG datastructure. The new RefVertex, however, has its
    //- own little mini-DAG which consists of a single DAGNode.
    //- Discuss this with the team before making it a part of the main DAG. 
    //- It is, however, deleted appropriately when the RefFace is deleted.
    //-***************************************************************

    /* other functions */

  Surface* get_surface_ptr() ;
  Surface const* get_surface_ptr() const ;
    //R Surface*
    //R- A pointer to the Surface to which the current 
    //R- face points. 
    //- This function returns a pointer to the Surface
    //- to which the current face points.
  
  CubitStatus get_graphics( GMem& results, 
                            unsigned short normal_tolerance = 15,
                            double distance_tolerance = 0.0,
                            double longest_edge = 0.0 );

protected :

  RefFace(Surface* surfacePtr) ;
    //- The constructor with a pointer to a Surface.
  
  DLIList<RefVertex*> hardPointList;

private:

  RefFace( const RefFace& );
  void operator=( const RefFace& );
  
  void initialize ();
    //- initialization method
  
#ifdef BOYD17 
  CubitBoolean amParametric;
#endif
  double maxPositionDeviation;
  
#ifdef BOYD17 
  int faceEdgeCount;
  int refFaceClone;
#endif
  int hardPointColor;

  double find_crack_length_no_uv();
  //- returns the crack length of the periodic surface (which can only have
  //- two loops, if there is some problem with the uv space.
     
};
#endif

