//-------------------------------------------------------------------------
//
// Filename      : RefEdge.hpp 
//
// Purpose       : This file contains the declarations of the class 
//                 RefEdge.
//
// Special Notes : 
//
// Creator       : Xuechen Liu 
//
// Creation Date : 07/11/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef REF_EDGE_HPP
#define REF_EDGE_HPP

#include "BasicTopologyEntity.hpp"
#include "CastTo.hpp"

const double maxSegmentLengthError = 0.005;

class CoEdge;

class CubitVector;
class FacetEvalTool;
class Curve;

class GMem;

// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT RefEdge : public BasicTopologyEntity
{
public :
   
// ********** BEGIN STATIC FUNCTION DECLARATIONS     **********
// ********** END STATIC FUNCTION DECLARATIONS     **********

  friend class RefEntityFactory;
    //- the factory is allowed to call the (private) constructors

  virtual ~RefEdge() ;
    //- The destructor

  static const char* get_class_name()
     {
       return "Curve";
     }

  virtual const char* class_name() const
     {
       return get_class_name();
     }
  
  DagType dag_type() const { return DagType::ref_edge_type(); }
  const type_info& entity_type_info() const { return typeid(RefEdge); }
  
  virtual CubitStatus get_point_direction( CubitVector& origin, 
                                           CubitVector& direction );
    //- Only valid for straight lines
    //- Finds the underlying line's origin and direction unit vector
    //- Returns CUBIT_FAILURE if curve is not a line

  virtual CubitStatus get_center_radius( CubitVector& center, double& radius );
    //- Only valid for arcs
    //- Finds the underlying arc's center point and radius
    //- Returns CUBIT_FAILURE if curve is not an arc
  
  // ********* TOPOLOGY ********************
  
  RefVertex* start_vertex();
  RefVertex* end_vertex();
    //R RefVertex*
    //R- Returned RefVertex pointer
    //- These functions get the start and end RefVertex'es of this
    //- RefEdge.
    //-
    //- MJP Notes:
    //- The start and end RefVertices are cached in the RefEdge (see
    //- private member data, start/end_RefVertex). These need to be
    //- updated if the end RefVertices could have changed during a
    //- geometric operation.  See the routine, RefEdge::update to
    //- see how this is done.
    
  virtual void reverse_topology();

  Chain* get_chain_ptr();
    //R Chain*
    //R- Pointer to the first Chain
    //- This function returns a pointer to the first Chain that this
    //- RefEdge points to.

  CubitStatus get_co_edges( DLIList<CoEdge*> &co_edges_found_list,
			    RefFace *input_ref_face_ptr = NULL );
    //R CubitStatus
    //R-CubitSucces/CubitFailure
    //O co_edges_found_list and optional input_ref_face_ptr
    //O- Gets the coedges that
    //- correspond to "this" RefEdge, they are populated into the
    //- co_edges_found_list. If an input_ref_face_ptr is sent it
    //- gets the co_edges that are just associated with that ref_face.
    //- Remember that usually there will be only one but for hard lines
    //- and sipes there may be two co-edges/ref-edge/ref-face.
    //- Note: This function will just blindly append the co-edges it
    //- finds into the co_edges_found_list.  A merge-unique might be
    //- more expensive than we want.

  double angle_between( RefEdge *other_edge_ptr, RefFace *face_ptr );
    //- Returns the "inside" angle between this edge and
    // other_edge_ptr.  Figures out which edge comes first.  Only
    // requirement is that both edges must be on face_ptr and on
    // the same loop.  Angle is in radians.
  
  CubitSense sense( RefFace *face );
    //- Sense of this edge wrt the given face.
    //- Return CUBIT_FORWARD or CUBIT_REVERSE.
    //- Returns CUBIT_UNKNOWN if there is more than one coedge, and
    //- senses don't agree, or if there are no coedges.
    
  virtual int dimension() const; 
    //- returns dimension of the actual entity.

  int num_of_common_ref_face( RefEdge *other_edge );
    // get the number of common refFace for two refEdges.

  RefFace *  common_ref_face (RefEdge* otherRefEdgePtr);
    //  - Find one common RefFace

  int common_ref_faces ( RefEdge* otherRefEdgePtr, DLIList<RefFace*> &common_face_list );
    //- Returns all common RefFaces that this edge shares with the input edge

  RefVertex* common_ref_vertex( RefEdge* otherRefEdgePtr );
    //- Find common RefVertex

  RefVertex* common_ref_vertex( RefEdge* next_ref_edge,
                                RefFace* ref_face_ptr );
    //- Finds the common vertex between the next_ref_edge and
    //- this edge, assuming that this edge comes before next_ref_edge
    //- in a the loop around the surface. This is needed because
    //- two edges may share two vertices...

  CubitBoolean common_vertices( RefEdge* otherRefEdgePtr, 
                                DLIList<RefVertex*> &common_verts);
    //-Populates the common_verts with the vertices that are common
    //-between the two curves.
    //-Returns CUBIT_TRUE if there are common vertices/ false otherwise.

  RefEdge* get_other_curve(RefVertex* common_vertex,
                           RefFace* ref_face_ptr);
    //- Finds the other curve on the ref_face_ptr that shares
    //- the common vertex
  
  CubitStatus get_two_co_edges( RefEdge *next_ref_edge,
                                RefFace *ref_face_ptr,
                                CoEdge *&co_edge_this,
                                CoEdge *&co_edge_next );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //O next_ref_edge, ref_face_ptr, co_edge_this, co_edge_next
    //O-Returns the co_edge that corrisponds to 'this' ref_edge
    //- and the one that corrisponds to next_ref_edge, with
    //- respect to the ref_face_ptr.
    //- Special Note : next_ref_edge must follow 'this' ref_edge in a Loop
    //- on the ref_face_ptr, this is assumed so the function
    //- will assert if this is not done...  

  RefFace *other_face(RefFace *not_face, RefVolume *ref_volume = NULL);
    //- return the (an) other face sharing this edge, which also borders
    //- ref_volume if non-NULL

  RefVertex* other_vertex( RefVertex* refVertexPtr );
    //- Returns the vertex at the other end of the edge.

  RefVertex* closest_vertex(const CubitVector &point3);
  
  // ********* GEOMETRY ********************
  
  Curve* get_curve_ptr() ;
  Curve const* get_curve_ptr() const ;
    //R Curve*
    //R- A pointer to the Curve to which this RefEdge points. 
    //- This function returns a pointer to the Curve
    //- to which the current RefEdge points.
      
  CubitVector start_coordinates();
  CubitVector end_coordinates();
    //R CubitVector
    //R- Returned location.
    //- These functions return the start and end global coordinate
    //- locations of the RefEdge.
    //-
    //- NOTE:
    //- These coordinates remain consistent throughout the life of the
    //- RefEdge, regardless of the fact that the actual start and
    //- end RefVerex'es may change (e.g., after a merge operation).
      
  virtual void move_to_curve ( CubitVector& vector );
    //- Moves the given location (CubitVector or CubitNode) to the closest 
    //- point on the Curve

  CubitStatus get_interior_extrema(DLIList<CubitVector*>& interior_points,
                                   CubitSense& return_sense) const;

  CubitStatus closest_point_trimmed( CubitVector const& location, 
                                     CubitVector& closest_location);
    //R void    
    //I location
    //I- The point to which the closest point on the Curve is desired.
    //O closest_location
    //O- The point on the Curve, closest to the input location which
    //O- will be on the Curve.  This is input as a reference 
    //O- so that the function can modify its contents.

  
  CubitStatus closest_point( CubitVector const& location, 
                             CubitVector& closest_location,
                             CubitVector* tangent_ptr = NULL,
                             CubitVector* curvature_ptr = NULL);
    //R void    
    //I location
    //I- The point to which the closest point on the Curve is desired.
    //O closest_location
    //O- The point on the Curve, closest to the input location which
    //O- might not be on the Curve.  This is input as a reference 
    //O- so that the function can modify its contents.
    //O tangent_ptr
    //O- The tangent to the Curve (output as a unit vector) at the 
    //O- closest_location.
    //O curvature_ptr
    //O- The curvature of the Curve at the closest_location.
    //- This function computes the point on the Curve closest to the input 
    //- location.
    //-
    //- If the tangent and/or curvature is required, then the calling code
    //- is responsible for allocating space for the CubitVector(s) and
    //- sending in the relevant non-NULL pointers.  If either of these
    //- pointers is NULL, the related quantity is not computed.
    //-
    //- Notes:
    //- The tangent direction is always in the positive direction of the 
    //- *owning RefEdge*, regardless of the positive direction of the
    //- underlying solid model entities, if any.

  CubitPointContainment point_containment( const CubitVector &point );
    //R CubitPointContainment - is the point on bounds of the curve?
    //R- CUBIT_PNT_OFF, CUBIT_PNT_ON, CUBIT_PNT_UNKNOWN
    //I CubitVector
    //I- position to check, known to be on the Curve
    // NOTE: POINT MUST LIE ON THE CURVE FOR THIS FUNCTION TO WORK PROPERLY.
      
  void tangent ( const CubitVector &point, CubitVector& tangent_vec );
    //- Return the tangent for the point on this RefEdge that is closest
    //- to the input "point".  The tangent direction is always in the
    //- positive direction of the RefEdge.  The positive direction of
    //- the RefEdge is an invariant through its lifecycle.


  CubitStatus tangent( const CubitVector &point,
                       CubitVector& tangent_vec,
                       RefFace *ref_face_ptr );
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE  
    //O point, tangent_vec, ref_face_ptr
    //O- Get the correct tangent with respect to
    //- the ref_face_ptr.
    //- This tangent function is not the safest method for getting
    //- the tangent on a surface.  In face if there could be
    //- another direction, this function could fail...(assert).
    //- This function assumes that there is only 1 co-edge per
    //- this ref_face for this REfEdge.
  
  CubitStatus tangent( const CubitVector &point, CubitVector& tangent_vec,
                       RefEdge *next_ref_edge, RefFace *ref_face_ptr );
    //- Retruns the tangent for the point on this RefEdge that is closest
    //- to the input "point".  Also, because the next_ref_edge and
    //- the ref_face_ptr are given, the tangent is oriented correctly.
    //- NOTE: It is assumed that 'this' RefEdge and next_ref_edge are
    //- in order in a "Loop" sense on the given ref_face_ptr.  The next
    //- ref_edge obviously must follow 'this' edge in this loop.
  
  virtual double measure();
    //- A generic geometric extent function.
    //- Returns volume for RefVolumes, area for RefFaces, length for RefEdge,
    //- and 1.0 for RefVertices
    //- A RefGroup calculates the maximum dimension of its contained
    //- entities and returns the sum of the measures() of all entities
    //- of that dimension.
    //- Default return value is 0.0 for all other RefEntities.

  virtual CubitString measure_label();
    //- Returns the type of measure: (volume, area, length, or N/A)

  double get_arc_length();
  double get_arc_length ( const CubitVector &point1,
                          const CubitVector &point2 );
  double get_arc_length ( const CubitVector &point1, int whichEnd );
    //- Various arc length calculations, some are redundant

  double get_chord_length();
    //- Calculates and returns the straight-line distance
    //- between the startRefVertex and the endRefVertex
  
  CubitVector curve_center();
  //- return a guess at the "centroid" of the curve, not the midpoint.
  //- For a non-periodic curves, returns the average of vertices
  //- For periodic curves, returns the average of the midpoint and the vertex.

  virtual CubitVector center_point();
    //- Returns location at the "actual" center of this RefEdge (along the 
    //- arc of the RefEdge) 

  CubitStatus mid_point ( const CubitVector &point1, const CubitVector &point2,
                          CubitVector& midPoint );
    //- Calculate midpoint between the 2 points on this RefEdge

  CubitStatus mid_point ( CubitVector &mid_point);
    //- Calculate midpoint on this RefEdge

  CubitStatus position_from_fraction( double fraction_along_curve,
                                      CubitVector &ouput_position );
    //R CubitStatus
    //I fraction in parameter space along refedge. (1/3,2/3,4/5...)
    //O- position where percent in parameter space lies.
    //-This function takes the given fraction, finds the parameter value,
    //-and calculates this position. This is based off the vgi curve.

  double start_param();
    //R double
    //- this function returns the starting parameter of the underlying curve
  
  double end_param();
    //R double
    //- this function returns the ending parameter of the underlying curve
  
  CubitBoolean get_param_range( double& start_param,
                                double& end_param );
    //R CubitBoolean
    //R- CUBIT_TRUE/FALSE
    //O start_param, end_param
    //O- The "lower" and "upper" parameter values of the RefEdge.
    //- This function returns the parameter bounds for the RefEdge.
    //- The start_param represents the parameter value of the 
    //- start location of the RefEdge and the end_param that of the
    //- end location.
    //-
    //- CUBIT_TRUE is returned if the RefEdge is defined parametrically
    //- and CUBIT_FALSE if it is not.  In the latter case, the 
    //- output values of start_ and end_param are undetermined.
    //-
    //- MJP Note:
    //- The numercial value of the start_param could be higher
    //- than that of the end_param.

  double u_from_position (const CubitVector& input_position);
    //R double
    //R- The returned "u" parameter value in local parametric space
    //I input_position
    //I- The input position for which "u" is to be computed.
    //- This function returns the coordinate of a point in the local
    //- parametric (u) space that corresponds to the input position in
    //- global (world) space.  The input point is first moved to the
    //- closest point on the Curve and the parameter value of that
    //- point is determined.

  CubitStatus position_from_u (double u_value,
                               CubitVector& output_position);
    //R CubitStatus
    //R- CUBIT_SUCCESS/FAILURE
    //I u_value
    //I- The input u parameter value
    //O output_position
    //O- The output position
    //- This function returns the coordinates of a point in the global
    //- (world) space that corresponds to the input parametric position 
    //- in the local space.
    //-
    //- If the input parameter value is not defined for the Curve, then 
    //- the input CubitVector is not modified and CUBIT_FAILURE is
    //- returned. Otherwise, position is appropriately modified and
    //- CUBIT_SUCCESS is returned.
    //-
    //- If the curve is periodic, the input u_value is first "normalized"
    //- to the fundamental period of the Curve before its position
    //- in global space is determined.

  double u_from_arc_length ( double root_param,
                             double arc_length );
    //R double
    //R- Returned parameter value
    //I root_param
    //I- The parameter value of the "root point"
    //I arc_length
    //I- A distance along the Curve
    //- This function returns the parameter value of the point that is
    //- "arc_length" away from the root point in the
    //- positive sense direction of the owning RefEdge.
    //- 
    //- A negative value for distance would force the search to go in the 
    //- negative (sense) direction of the RefEdge.
    //-
    //- NOTE:
    //- The important assumption that is made in this routine is that
    //- the end points of the RefEdge that owns this CurveACIS are the same
    //- as the end points of the first solid model entity in the list of 
    //- solid model entities associated with this Curve.
      
  double fraction_from_arc_length(RefVertex *root_vertex,
                                  double     length);

  CubitStatus point_from_arc_length ( double root_param, 
                                      double arc_length,
                                      CubitVector& new_point );
    //- See below. This allows you to directly specify the starting parameter,
    //- instead of it being determined from the root_point.
 
  CubitStatus point_from_arc_length ( const CubitVector& root_point, 
                                      double arc_length,
                                      CubitVector& new_point );
    //- Return a point arc_length distance from root_point
    //- on this RefEdge. 
    //-
    //- If arc_length is negative, the new point
    //- is in the negative sense direction (along the RefEdge) from
    //- the root point.
    //- 
    //- If the curve is not periodic and the point arc_length away
    //- from root_point in the appropriate direction goes beyond
    //- the end point of the RefEdge, that end point is returned
    //- as new_point.
    //-
    //- If the curve is periodic and the point arc_length away
    //- from root_point in the appropriate direction goes beyond
    //- the end point of the RefEdge, wrap around is done.
    //-
    //- NOTE: I have had problems with this function if the root point
    //- is the start vertex, for some reason by a factor of 1e-7 the
    //- parameter value is different between the start_param and the
    //- parameter found from the start_vertex location.  So I will switch
    //- EdgeMeshTool to go off the start parameter instead of this function.
  

  double length_from_u( double parameter1,
                        double parameter2 );

  CubitBoolean is_periodic();
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //- This function determines whether the underlying geometry of the
    //- Curve is periodic or not.  Returns CUBIT_TRUE if it is and 
    //- CUBIT_FALSE if it is not.
  
  CubitBoolean is_periodic( double& period);
    //R CubitBoolean
    //R- CUBIT_TRUE/CUBIT_FALSE
    //O period
    //O- Returned period value
    //- This function determines whether the underlying geometry of the
    //- Curve is periodic or not.  Returns CUBIT_TRUE if it is and 
    //- CUBIT_FALSE if it is not.
    //-
    //- If it is periodic, then it returns the period in the input
    //- reference variable, "period". This value is set to 0.0 if
    //- the Curve is not periodic.

  int get_mark();
  void set_mark( int set_value );
    //- get/set generic 2-bit mark.

  CubitStatus relative_sense( RefEdge *ref_edge_ptr_2,
                              double tolerance_factor,
                              CubitSense *sense,
                              CubitBoolean &spatially_equal,
                              CubitBoolean force_merge = CUBIT_FALSE);
    //- calculates the relative sense of two refedges and tries to
    //- see if there are two points that are close on each edge (1/3 and 2/3
    //- param value).


  CubitBoolean about_spatially_equal( RefEdge* ref_edge_ptr_2,
                                      double tolerance_factor = 1.0,
                                      CubitSense* sensePtr = NULL,
                                      CubitBoolean notify_refEntity =
                                      CUBIT_FALSE );

    //R CubitBoolean
    //R-CUBIT_TRUE/CUBIT_FALSE
    //I RefEdge, double, CubitSense*, CubitBoolean
    //I- Second RefEdge to compare, Tolerance factor to for GEOMETRY_RESABS,
    //I- and flag for notifying compared RefEntities.
    //O CubitSense*, and returned CubitBoolean.
    //O- if the two refEdges are spatially equal within the GEOMETRY_RESABS*
    //- the tolerance_factor, then CUBIT_TRUE will be returned.  Otherwise
    //- CUBIT_FALSE is returned.  If there is a non-null sesnePtr sent in,
    //- the sense pointer value will be assinged the relative sense of the
    //- RefEdges, ie, CUBIT_FORWARD/CUBIT_REVERSED.
    //- To do the comparison, the end points, 1/3 and 2/3 world
    //- points are spatially compared.

  virtual int validate();
    //- Check that entity is valid. Returns number of problems detected.

//=========  Add Code by SRS of Cat,  3/3/99 2:28:31 PM  =========
  CubitBoolean is_tolerant();
    //- Returns CUBIT_TRUE if refedge is a tolerant edge.  Tolerant
    //- edges can be created in the ACIS healer if it is unable to
    //- heal the edge.
//=========  Code End by SRS of Cat,  3/3/99 2:28:31 PM  =========

  CubitStatus get_graphics( GMem& polyline, double tolerance = 0.0 );
    
  void reverse_tangent();

  static void suppress_edge_length_warning(bool flag);


protected:

  RefEdge(Curve* curvePtr) ;
    //- The constructor with a pointer to a Curve.

private:
  void initialize();
    //- Initializes all member data

  double get_lower_param_bound();
  double get_upper_param_bound();
    //- Get the parameter values at the start and end locations of
    //- the RefEdge.  If it is periodic, this is taken into account
    //- and the lower and upper bounds of the period are returned.

  int refEdgeClone;
   
  cBit markedFlag : 2; //- generic flag, one or two bits?
  
  RefEdge( const RefEdge& );
  void operator=( const RefEdge& );
  static bool mSuppressEdgeLengthWarning;
};

// ********** BEGIN HELPER CLASSES         **********
// ********** END   HELPER CLASSES         **********

// ********** BEGIN INLINE FUNCTIONS       **********

inline int
RefEdge::get_mark()
{ return markedFlag; }

inline void
RefEdge::set_mark( int set_value )
{ markedFlag = set_value; }
  

// ********** END INLINE FUNCTIONS         **********
 
// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********
 
// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********
 
#endif

