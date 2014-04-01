//-------------------------------------------------------------------------
// Filename      : GeoemtricModifyEngine.hpp
//
// Purpose       : Define the interface for all solid model modify
//                 engines.
//
// Special Notes : This is an abstract base class.
//
// Creator       : Tim Tautges
//
// Creation Date : 2/01
//
// Owner         : Tim Tautges
//-------------------------------------------------------------------------

#ifndef GEOMETRY_MODIFY_ENGINE_HPP
#define GEOMETRY_MODIFY_ENGINE_HPP

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitGeomConfigure.h"

class CubitVector;
class CubitPlane;
template <class X> class DLIList;


class TopologyBridge;
class Point;
class Curve;
class Surface;
class Lump;
class BodySM;
class LoopSM;
class GeometryEntity;
class GeometryQueryEngine;

class CUBIT_GEOM_EXPORT GeometryModifyEngine
{

   public:

  virtual ~GeometryModifyEngine() {}
    //- virtual destructor

  virtual Point* make_Point( CubitVector const& position ) const = 0;
    //R Point*
    //R- Returned pointer to a Point object.
    //I position
    //I- Input coordinates of the point to be created.
    //- This function creates a Point, given coordinates.  The particular
    //- type of Point object that is created depends on the specific
    //- modeling engine.  For example, if the engine
    //- is AcisGeometryEngine, then a PointACIS is created and returned.

  virtual Curve* make_Curve(Curve *curve_ptr) const = 0;
    //- creates a curve from an existing curve.  This creates totally
    //- new topology.  This function is useful for constructing geometry
    //- from existing geometry.
  
  virtual Curve* make_Curve( Point const* point1_ptr,
                             Point const* point2_ptr,
                             Surface* ref_face_ptr,
                             const CubitVector *third_point = NULL) const = 0;
    //- Create a curve exactly on the give ref_face.
    //- Make sure the points are on the underlying surface.
  
  virtual Curve* make_Curve( GeometryType curve_type,
                             Point const* point1_ptr,
                             Point const* point2_ptr,
                             DLIList<CubitVector*>& vector_list,
                             Surface* ref_face_ptr = NULL) const = 0;
  
  virtual Curve* make_Curve( GeometryType curve_type,
                             Point const* point1_ptr,
                             Point const* point2_ptr,
                             CubitVector const* intermediate_point_ptr,
                             CubitSense sense) const = 0;
    //R Curve*
    //R- Returned pointer to a Curve object.
    //I curve_type 
    //I- The type of curve to be created.
    //I point1_ptr, point2_ptr
    //I- Input end points of the curve to be created.
    //I intermediate_point_ptr
    //I- The coordinates of an intermediate point required for the
    //I- generation of the curve.
    //I vector_list
    //I- Input list of CubitVectors. The new Curve interpolates these
    //I- positions.
    //I ref_face_ptr
    //I- If the optional Surface pointer is provided, then the input
    //I- locations (CubitVectors) are moved to the surface of the 
    //I- Surface before interpolation is done.
    //I- NOTE: The end Points are *not* moved to the surface. Only the
    //I-       (intermediate) points are.
    //- This function creates a Curve, of type, curve_type, given 
    //- the end points of the curve and an intermediate point or a set
    //- of points to interpolate through.
    //- The particular type of Curve object 
    //- that is created depends on the specific modeling engine.  For 
    //- example, if the engine is AcisGeometryEngine, then a CurveACIS 
    //- is created and returned.
    //-
    //- NOTE on some hideous programming
    //-   This function not only creates the Curve, but also
    //-   creates the VGI entities Chain and CoVertexes that
    //-   will be connected up with the RefEdge. For the sake
    //-   of expediency we are assuming that this function is
    //-   called when RefEdge needs to be made from RefVertices.
    //-
    //- This function creates non-linear curves (e.g., elliptic and 
    //- parabolic curves. curve_type indicates what type of curve is
    //- required.
    //- 
    //- Elliptical:
    //- The portion of the ellipse (circle) that is generated goes from
    //- point1 to point2.
    //-
    //- Parabolic:
    //- Construct a parabolic arc from 3 points. intermediate_point is the 
    //- peak of the parabola - in this case, the point which is equidistant 
    //- from the start and end points of the parabola. The 3 points must form
    //- an isosceles triangle. This definition limits the user to generation 
    //- of the tip of parabolic shapes only.
  
  virtual Surface* make_Surface( Surface *old_surface_ptr,
                                 bool extended_from = false) const = 0;
    //R Surface*
    //R- Pointer to a newly created Surface object.
    //I Surface*
    //I- The surface from which we want to create a new one.
    //- This function creates a new surface from an existing one.
    //- The new surface is attached to ACIS geometry.  The acis
    //- goemetry is attached to a full data structure, loop, lump, bodies..
  
  virtual Surface* make_Surface( GeometryType surface_type, 
                                 DLIList<Curve*>& curve_list,
                                 Surface *old_surface_ptr = NULL,
                                 bool check_edges = true ) const = 0;
    //R Surface*
    //R- Pointer to a newly created Surface object.
    //I surface_type
    //I- The type of surface to be created.
    //I curve_list
    //I- list of curves to be used as the bounds of the Surface
    //- This function creates a Surface, of type, surface_type, given 
    //- the list of curves.
    //- The particular type of Surface object created depends on the 
    //- specific modeling engine.  For example, if the engine is 
    //- AcisGeometryEngine, then a SurfaceACIS is created.
  
  virtual CubitStatus stitch_surfs( 
			DLIList<BodySM*> &surf_bodies,
                        BodySM *& stitched_Body )const= 0;
   //I List of surface_bodys you want stitched together
   //IO If stitching of FACEs is successful, stitched_Body is resultant BodySM.
   //   Means that all surfaces could be stitched together, forming a single BODY,
   //   but not a single FACE.
   //returns CUBIT_SUCCESS if stitching was successful, other return FAILURE

  virtual Lump* make_Lump( DLIList<Surface*>& surface_list ) const = 0;
    //R Lump*
    //R- Pointer to a newly created Lump object.
    //I lump_type
    //I- The type of surface to be created.
    //I surface_list
    //I- list of surfaces to be used as the bounds of the Lump
    //- This function creates a Lump, of type, lump_type, given 
    //- the list of surfaces.
    //- The particular type of Lump object created depends on the 
    //- specific modeling engine.  For example, if the engine is 
    //- AcisGeometryEngine, then a LumpACIS is created.
    //-
    //- NOTE on some hideous programming
    //-   This function not only creates the Lump, but also
    //-   creates the VGI entities Shell and CoFaces that
    //-   will be connected up with the RefVolume. For the sake
    //-   of expediency we are assuming that this function is
    //-   called when RefVolume needs to be made from Surfaces.
    //-

  virtual BodySM* make_BodySM( Surface * ) const = 0;
  
  virtual BodySM* make_BodySM( DLIList<Lump*>& /*lump_list*/ ) const = 0;
    //R BodySM*
    //R- Pointer to a newly created BodySM object.
    //I lump_list
    //I- list of lumps to be used to create the BodySM
    //- This function creates a BodySM given the list of lumps.
    //- The particular type of BodySM object created depends on the 
    //- specific modeling engine.  For example, if the engine is 
    //- AcisGeometryEngine, then a BodyACIS is created. Non-solid
    //- model based engines can use the default implementation 
    //- to return a NULL pointer.
    //-
    //- NOTE on some hideous programming
    //-   This function not only creates the Lump, but also
    //-   creates the VGI entities CoVolumes that
    //-   will be connected up with the Body. For the sake
    //-   of expediency we are assuming that this function is
    //-   called when Body needs to be made from RefVolumes.
    //-
  
    //HEADER- Functions for creation of primitives
      virtual BodySM* sphere(double radius) const = 0 ;
      //R Body*
      //R- A pointer to a newly created Body
      //I radius
      //I- Radius of the sphere
      //- Creates a sphere and assigns it to a Body.
      //- Returns pointer to the newly created body.

      virtual BodySM* brick ( double wid, double dep, double hi ) const = 0 ;
      //R Body*
      //R- A pointer to a newly created Body
      //I wid
      //I- Width of the brick
      //I dep
      //I- Depth of the brick
      //I hi
      //I- Height of the brick
      //- Creates a cuboid and assigns it to a Body.
      //- Returns pointer to the newly created body.

      virtual BodySM* brick( const CubitVector &center, 
                             const CubitVector axes[3],
                             const CubitVector &extension) const = 0 ;
      //R Body*
      //R- A pointer to a newly created Body
      //I center
      //I- Center location of the brick
      //I axes
      //I- XYZ axes of brick
      //I extension
      //I- Size of brick, equivalent to 1/2 width, height, depth
      //- Creates a cuboid and assigns it to a Body.
      //- Returns pointer to the newly created body.

      virtual BodySM* prism( double height, int sides, double major, 
                           double minor) const = 0 ;
      //- Creates an ACIS prism and assigns it to a Body $
      //- {height, major, minor} input height, major and minor radii. $
      //- {sides} input number of sides. Must be >= 3.
      //- Returns the ID of the new Body or CUBIT_FAILURE

      virtual BodySM* pyramid( double height, int sides, double major, 
                             double minor, double top=0.0) const = 0 ;
      //- Creates an ACIS pyramid and assigns it to a Body $
      //- {height, major, minor} input height, major and minor radii. $
      //- {sides} input number of sides. Must be >= 3.
      //- {top} radius at top of pyramid.
      //- Returns the ID of the new Body or CUBIT_FAILURE

      virtual BodySM* cylinder( double hi, double r1, double r2, 
                                double r3) const = 0;
      //- Creates an ACIS frustum and assigns it to a Body $
      //- {hi} input height $
      //- {r1} input radius in x-direction at base $
      //- {r2} input radius in y-direction at base $
      //- {r3} input radius in x-direction at top
      //- Returns the ID of the new Body or CUBIT_FAILURE

      virtual BodySM* torus( double r1, double r2 ) const = 0 ;
      //- Creates an ACIS torus and assigns it to a Body $
      //- {r1} input major_radius $
      //- {r2} input minor_radius
      //- Returns the ID of the new Body or CUBIT_FAILURE 

      virtual BodySM* planar_sheet ( const CubitVector& p1,
                                     const CubitVector& p2,
                                     const CubitVector& p3,
                                     const CubitVector& p4 ) const = 0;
      //- Creates a solid body consisting of a planar sheet (no volume)
      //- {p1} - 1st corner of the sheet
      //- {p2} - 2nd corner of the sheet
      //- {p3} - 3rd corner of the sheet
      //- {p4} - 4th corner of the sheet

      virtual BodySM* copy_body ( BodySM* bodyPtr) const = 0 ;
      //R Body*
      //R- A pointer to the newly created body
      //I bodyPtr
      //I- A pointer to the body to be copied
      //- This function makes a copy of the input Body and returns a
      //- pointer to the newly created copy. The input Body and the newly
      //- created Body are geometrically identical.

      //HEADER- Functions for boolean operations

      virtual CubitStatus subtract(DLIList<BodySM*> &tool_body_list,
                                   DLIList<BodySM*> &from_bodies,
                                   DLIList<BodySM*> &new_bodies,
                                   bool imprint = false, 
                                   bool keep_old = false) const = 0;
      /********************************************************************************* 
      virtual CubitStatus imprint(BodySM* BodyPtr1, BodySM* BodyPtr2,
                                  BodySM*& newBody1, BodySM*& newBody2,
                                  bool keep_old) const = 0;
      *********************************************************************************/
      virtual CubitStatus imprint(DLIList<BodySM*> &from_body_list,
                                  DLIList<BodySM*> &new_from_body_list,
                                  bool keep_old,
                                  DLIList<TopologyBridge*> *new_tbs = NULL,
                                  DLIList<TopologyBridge*> *att_tbs = NULL) const = 0;

      virtual CubitStatus imprint( DLIList<BodySM*> &body_list,
                                   DLIList<Curve*> &ref_edge_list,
                                   DLIList<BodySM*>& new_body_list,
                                   bool keep_old_body,
                                   bool show_messages= true) const = 0;
      //- Imprints a list of Bodies with a list of RefEdges.  All
      //- entities must be ACIS entities.  Useful for splitting
      //- surfaces.  If edge pierces a surface a hardpoint will
      //- result at the pierce location.

      virtual CubitStatus imprint( DLIList<Surface*> &ref_face_list,
                                   DLIList<Curve*> &ref_edge_list,
                                   DLIList<BodySM*>& new_body_list,
                                   bool keep_old_body ) const = 0;
      //- Imprints a list of Surfaces with a list of RefEdges.  This is
      //- useful if the user has a curve which spans several surfaces on 
      //- a body and only wants to imprint to selected surfaces.  Algorithm 
      //- does not support imprinting to free surfaces.

      virtual CubitStatus imprint( DLIList<Surface*> &surface_list,
                                   DLIList<DLIList<Curve*>*> &curve_lists_list,
                                   BodySM*& new_body,
                                   bool keep_old_body ) const = 0;
      //- Imprints a list of Surfaces with list of Curves, sorted per
      //- Surface (ie., curve_lists_list is same length as surface_list).
      //- All input surfaces must be from the same body.

      virtual CubitStatus imprint( DLIList<BodySM*> &body_list,
                                   DLIList<CubitVector*> &vector_list,
                                   DLIList<BodySM*>& new_body_list,
                                   bool keep_old_body,
                                   DLIList<TopologyBridge*> *new_tbs = NULL,
                                   DLIList<TopologyBridge*> *att_tbs = NULL ) const = 0;
      //- Imprints a list of bodies with a list of vectors.  Useful for
      //- splitting curves and creating hardpoints on surfaces.

      virtual CubitStatus tolerant_imprint( DLIList<BodySM*> &bodies_in,
                                            DLIList<BodySM*> &new_bodies,
                                   DLIList<TopologyBridge*> *new_tbs = NULL,
                                   DLIList<TopologyBridge*> *att_tbs = NULL ) const = 0;

      virtual CubitStatus imprint_projected_edges( 
                                       DLIList<Surface*> &ref_face_list,
                                       DLIList<Curve*> &ref_edge_list,
                                       DLIList<BodySM*>& new_body_list,
                                       bool keep_old_body,
                                       bool keep_free_edges) const= 0  ;
      //- Imprints a list of Surfaces with a list of projected RefEdges.  

      virtual CubitStatus imprint_projected_edges(
                                       DLIList<Surface*> &ref_face_list,
                                       DLIList<BodySM*> &body_list,
                                       DLIList<Curve*> &ref_edge_list,
                                       DLIList<BodySM*>& new_body_list,
                                       bool keep_old_body,
                                       bool keep_free_edges) const= 0  ;
      //- Imprints a list of Bodies with a list of RefEdges which are projected
      //- to a list of Surfaces

      virtual CubitStatus project_edges( 
                             DLIList<Surface*> &ref_face_list,
                             DLIList<Curve*> &ref_edge_list_in,
                             DLIList<Curve*> &ref_edge_list_new,
                             bool print_error = true ) const = 0;
      //- Projects list RefEdges to a list of Surfaces
  
      virtual CubitStatus intersect(BodySM* tool_body_ptr,
                                    DLIList<BodySM*> &from_bodies,
                                    DLIList<BodySM*> &new_bodies,
                                    bool keep_old = false) const = 0;
      //R CubitStatus
      //R- CUBIT_SUCCESS/CUBIT_FAILURE
      //I BodyPtr1, BodyPtr2
      //I- Two Body pointers that will be Booleaned
      //O newBody (or newBody1, newBody2 for imprint)
      //O- The new Body build by boolean operation on two Bodys.
      //O- for imprint, the output is two Bodies.
      //- These functions perform boolean operations on two Bodys and 
      //- return the result through the output arguments. If the
      //- boolean operations fails at any stage, a NULL value is assigned 
      //- to the output argument(s) and the function returns 
      //- CUBIT_FAILURE. If everything goes well, the function returns 
      //- CUBIT_SUCCESS. The original Bodys are left untouched during the
      //- boolean operation. The boolean operation is carried out on 
      //- copies of the original Bodys.

	  virtual CubitStatus chop( DLIList<BodySM*> &bodies, 
                              DLIList<BodySM*> &intersectBodies,
                              DLIList<BodySM*> &outsideBodies, 
                              BodySM*& leftoversBody,
                              bool keep_old = false,
                              bool nonreg = false) const = 0;
    //R CubitStatus
    //R-the result of the chop operation: Success or Failure
    //I bodies
    //I-DLIList<Body*>: a list of Body pointers that will be united
    //O intersectBody, outsideBody, leftoversBody
    //O- new Bodies build by chop operation on the list of  Body pointers.
    //- This function performs a chop of a tool body on a blank body  and returns 
    //- the result through the output argument intersectBody, outsideBody, leftoversBody. If the chop
    //- operation went through OK, the function returns CUBIT_SUCCESS. If,
    //- for some reason, the chop operation did not go well, the output
    //- argument is assigned a NULL value and the function returns 
    //- CUBIT_FAILURE. In either case, the original Bodys are left 
    //- untouched. 
	  

      virtual CubitStatus thicken( DLIList<BodySM*>& bodies, 
                                   DLIList<BodySM*>& new_body_list,
                                   double depth,
                                   bool both = false) const = 0;
    //R CubitStatus
    //R-the result of the thicken operation: Success or Failure
    //I bodies
    //I-DLIList<Body*>: a list of Body pointers that will be thicken
    //O- new Bodies build by thicken operation on the list of  Body pointers.
    //- This function performs a thicken of sheet bodies and returns 
    //- the result through the output argument in_out_body. If the thicken
    //- operation went through OK, the function returns CUBIT_SUCCESS. If,
    //- for some reason, the thicken operation did not go well, the output
    //- argument is assigned a NULL value and the function returns 
    //- CUBIT_FAILURE.

     virtual CubitStatus hollow(DLIList<BodySM*>& bodies,
                                DLIList<Surface*>& surfs_to_remove,
                                DLIList<BodySM*>& new_bodies,
                                double depth) const = 0;
     //R CubitStatus
    //R-the result of the hollow operation: Success or Failure
    //I body
    //I- for OCC: a list of solid BodySM's that will be hollowed into a list of
    //I- thick solids.
    //I- surfs_to_remove: the faces to be removed from the original solid (OCC).
    //O- new Bodies build by hollow operation on the  body pointers.
    //- This function performs a hollow of solid body and returns
    //- the result through the output argument in_out_body. If the hollow
    //- operation went through OK, the function returns CUBIT_SUCCESS. If,
    //- for some reason, the hollow operation did not go well, the output
    //- argument is assigned a NULL value and the function returns
    //- CUBIT_FAILURE.

     virtual CubitStatus flip_normals( DLIList<Surface*>& face_list ) const = 0;
            //R CubitStatus
    //R-the result of the flip_normal operation: Success or Failure
    //I bodies
    //I-DLIList<Surface*>: a list of Face pointers that will be fliped
    //O- If the flip operation went through OK, the function returns CUBIT_SUCCESS. If,
    //- for some reason, the flip operation did not go well, the output
    //- returns CUBIT_FAILURE.
	  
	  
	  virtual CubitStatus unite( DLIList<BodySM*> &bodies, 
                               DLIList<BodySM*> &newBodies,
				                       bool keep_old = false) const = 0;
      //R CubitStatus
      //R- CUBIT_SUCCESS/CUBIT_FAILURE
      //I bodies
      //I- A list of Bodys that will be united
      //O newBody 
      //O- The new Body built by the unite operation on the list of Bodys.
      //- This function performs a unite of a list of Bodys and returns 
      //- the result through the output argument, "newBody". If the unite
      //- operation went through OK, the function returns CUBIT_SUCCESS. If,
      //- for some reason, the unite operation did not go well, the output
      //- argument is assigned a NULL value and the function returns 
      //- CUBIT_FAILURE. In either case, the original Bodys are left 
      //- untouched.

      //HEADER- Sweep-related functions. All of these are implemented only for 
      //HEADER- RefEntities whose underlying geometry is represented by a solid
      //HEADER- model such as ACIS. 
      
      virtual CubitStatus  sweep_translational(
                                      DLIList<GeometryEntity*>& ref_ent_list,
                                      DLIList<BodySM*>& result_body_list,
                                      const CubitVector& sweep_vector,
                                      double draft_angle,
                                      int draft_type,
                                      bool switchside,
                                      bool rigid,
                                      Surface* stop_surf = NULL,
                                      BodySM*  stop_body = NULL) const = 0;
      
      virtual CubitStatus  sweep_perpendicular(
                                      DLIList<GeometryEntity*>& ref_ent_list,
                                      DLIList<BodySM*>& result_body_list,
                                      double distance,
                                      double draft_angle,
                                      int draft_type,
                                      bool switchside,
                                      bool rigid,
                                      Surface* stop_surf = NULL,
                                      BodySM*  stop_body = NULL ) const = 0;
      
      virtual CubitStatus  sweep_rotational(
                                   DLIList<GeometryEntity*>& ref_ent_list,
                                   DLIList<BodySM*>& result_body_list,
                                   const CubitVector& point,
                                   const CubitVector& direction,
                                   double angle,
                                   int steps = 0,
                                   double draft_angle = 0.0,
                                   int draft_type = 0,
                                   bool switchside = false,
                                   bool make_solid = false,
                                   bool rigid = false,
                                   Surface *stop_surf = NULL,
                                   BodySM*  stop_body = NULL ) const = 0;
        
      virtual CubitStatus sweep_along_curve( 
                                   DLIList<GeometryEntity*>& ref_ent_list,
                                   DLIList<BodySM*>& result_body_list,
                                   DLIList<Curve*>& ref_edge_list,
                                   double draft_angle = 0.0,
                                   int draft_type = 0,
                                   bool rigid = false,
                                   Surface *stop_surf = NULL,
                                   BodySM*  stop_body = NULL  ) const= 0;
     

      virtual CubitStatus scale( BodySM *&body, const CubitVector& factors ) = 0;

      //HEADER- Webcut-related functions
      virtual CubitStatus webcut( 
              DLIList<BodySM*>& webcut_body_list, 
              const CubitVector &v1,
              const CubitVector &v2,
              const CubitVector &v3,
              DLIList<BodySM*>& results_list,
              bool imprint = false             
 )  = 0;
      //R int
      //R- Number of bodies that were webcut ( >= 0 )
      //I webcut_body_list
      //I- The list of bodies to be webcut
      //I plane
      //I- The plane to be used for webcutting.
      //I merge
      //I- A flag to decide whether the new bodies created by the
      //I- webcutting process should be merged or not.
      //I imprint
      //I- A flag to decide whether the new bodies created by the
      //I- webcutting process should be imprinted or not.
      //- This function webcuts a list of bodies through a plane.
      //- The newly created bodies are merged and imprinted depeding on 
      //- the respective flags.

      virtual CubitStatus webcut(
                    DLIList<BodySM*>& webcut_body_list, 
                    BodySM const* tool_body,
                    DLIList<BodySM*>& results_list,
                    bool imprint = false
                    ) = 0 ;
      //R int       
      //R- Number of bodies that were webcut ( >= 0 )
      //I webcut_body_list
      //I- The list of bodies to be webcut
      //I tool_body
      //I- The body to be used for webcutting.
      //I merge
      //I- A flag to decide whether the new bodies created by the
      //I- webcutting process should be merged or not.
      //I imprint
      //I- A flag to decide whether the new bodies created by the
      //I- webcutting process should be imprinted or not.
      //- This function webcuts a list of bodies using another body
      //- as the webcutting tool. The newly created bodies are 
      //- merged and imprinted depeding on the respective flags.
      virtual CubitStatus webcut_across_translate( 
                                              DLIList<BodySM*>& body_list, 
                                              Surface* plane_surf1,
                                              Surface* plane_surf2,
                                              DLIList<BodySM*>& results_list, 
                                              bool imprint = false) const = 0;
     
	  	  
	  //- R status
      //- R-results_list of bodies affected, or created from webcut.
      //- I- Bodies to be webcut, plane to define cuts, and imprint merge flags.
      //- This is an experimental function, hooked to the GUI for making
      //- bodies one to one sweeps.
  
      virtual CubitStatus section( DLIList<BodySM*> &section_body_list,
                                   const CubitVector &point_1,
                                   const CubitVector &point_2,
                                   const CubitVector &point_3,
                                   DLIList<BodySM*>& new_body_list,
                                   bool keep_normal_side,
                                   bool keep_old = false,
                                   bool keep_both_sides = false) = 0;
      //- Section will cut a list a bodies and keep a side of the bodies.
      //- The bodies are cut with a planar surface (surface will be extended).
  
      virtual CubitStatus split_body( BodySM *body_ptr,
                                      DLIList<BodySM*> &new_bodies ) = 0;
      //- Splits a body with multiple volumes into multiple bodies
      //- each having only one volume.
  
      virtual CubitStatus reverse_body( BodySM *body_to_reverse ) = 0;
      //- Reverse body (turn it inside-out).
  
      virtual CubitStatus split_periodic( BodySM *body_ptr,
                                          BodySM *&new_body ) = 0;
      //- Splits the periodic surfaces along the 0 and PI/2's periods.

      virtual CubitStatus regularize_body( BodySM *body_ptr,
                                           BodySM *&new_body_ptr ) = 0;
      //- Removes all unnecessary faces, edges and vertices from the body.   
      
	    virtual CubitStatus regularize_entity( GeometryEntity *old_entity_ptr,  
                                             BodySM *&new_body_ptr ) = 0;
      //- Removes all all unnessesary faces, curves, vertices and associated
      //- data from a refentity.
   

      // Added by CAT
      virtual CubitStatus offset_curves( DLIList<Curve*>& ref_edge_list, 
                                         DLIList<Curve*>& result_curve_list,
                                         double offset_distance, 
                                         const CubitVector& offset_direction, 
                                         int gap_type = 1 ) = 0;
      //- Creates curves offset from a chain of curves.  The offset direction is
      //- only used if there is one linear curve.  Otherwise, the offset direction
      //- is calculated by ACIS (the cross product of the wires tangent and the 
      //- planar normal).  The gap type is 0 - rounded, 1 - extended, 2 - natural.

      virtual Curve* trim_curve( Curve* trim_curve, 
                                      const CubitVector& trim_vector, 
                                      const CubitVector& keep_vector,
                                      bool keep_old = false ) = 0;
      //- Trims or extends a curve, up to the trim_vector.  If trimming, the 
      //- keep_vector determines which side of the curve to keep.  If the curve 
      //- is not free, the curve is automatically copied before trimming (so
      //- a new curve results).

      virtual CubitStatus create_solid_bodies_from_surfs( 
                                      DLIList<Surface*> &ref_face_list, 
                                      DLIList<BodySM*> &new_bodies,
                                      bool keep_old = false,
                                      bool heal = true) const = 0;
      //- Creates a single body from a set of faces.  The faces can only be attached
      //- to bodies if they are sheet bodies.  It is assumed that the calling code will
      //- check for this, ie GeometryTool.

      virtual Curve* create_arc_three( Point* pt1, 
                                       Point* pt2,
                                       Point* pt3,
                                       bool full = false ) = 0;
      virtual Curve* create_arc_three( Curve* curve1,
                                       Curve* curve2,
                                       Curve* curve3,
                                       bool full = false ) = 0;
  
      virtual Curve* create_arc_center_edge( Point* point1,
                                             Point* point2,
                                             Point* point3,
                                             const CubitVector &normal,
                                             double radius = CUBIT_DBL_MAX,
                                             bool full = false ) = 0;

      //- Methods to create arcs.  First uses 3 points on arc, next creates arc
      //- tangent to 3 curves, last creates arc using center and two points on arc.
      //- If full option is specified, a full circle is created.

      virtual CubitStatus create_curve_combine( DLIList<Curve*>& curve_list, 
                                    Curve *&new_curve_ptr ) = 0;
      //-  Uses the solid modeller to create a new RefEdge that is a combination 
      //- of the input chain of edges.  
      //-
  
      virtual GeometryQueryEngine *get_gqe() = 0;
        /**< all gme's should be able to return a gqe
         */

      virtual CubitBoolean is_modify_engine(const TopologyBridge *) const 
        {return CUBIT_FALSE;};
        /**< return CUBIT_TRUE if the tb_ptr belongs to this modify engine
         */

      virtual CubitStatus get_offset_intersections( 
                               Curve* ref_edge1, 
                               Curve* ref_edge2,
                               DLIList<CubitVector*>& intersection_list,
                               double offset = 0.0, 
                               bool ext_first = true ) = 0;
    //- Finds the intersections of a certain distance (offset) between two
    //- curves.  The two curves must lie in a plane.  The first curve is offset 
    //- the offset distance in both directions, and the bounded intersections with 
    //- the second curve are found.  The first curve and optionally be extended
    //- to infinity for the intersection calculation.  The intent of the function
    //- is so that the user can create a point on a curve a certain distance
    //- from another curve, as in specifying a reference location for a gage 
    //- diameter on an arc in an engineering drawing.

  virtual CubitStatus get_offset_intersections( 
                               Curve* ref_edge_ptr, 
                               Surface* ref_face_ptr,
                               DLIList<CubitVector*> &intersection_list,
                               double offset = 0.0, 
                               bool ext_surf = true ) = 0;
    //- Finds intersections (points) of the curve and surface.  The surface can 
    //- be offset - it is offset to each side and intersections are found.  By
    //- default the surface is extended to infinity (if possible) and the 
    //- intersections are found.  The function allocates the CubitVectors in 
    //- the returned list, so be sure to free them.
  
  virtual CubitStatus surface_intersection( Surface *surface1_ptr,
                                            Surface *surface2_ptr,
                                            DLIList<Curve*> &inter_graph,
                                            const double tol) const = 0;
    //- Finds the intersection of two bodies defined by the curves of the
    //- intersection graph.

  virtual CubitStatus get_mid_plane( const CubitVector &point_1,
		                     const CubitVector &point_2,
                                     const CubitVector &point_3,
		                     BodySM *body_to_trim_to,
		                     DLIList<BodySM*>& midplane_bodies ) const = 0;
    //- Finds the mid plane described by the 3 points and trims
    //- it to the body.  It returns the mid planes as part of the 
    //- midplane_body

  virtual CubitStatus get_spheric_mid_surface( Surface *surface_ptr1,
					       Surface *surface_ptr2,
					       BodySM *body_to_trim_to,
					       DLIList<BodySM*>& midsurface_bodies ) const = 0;

  virtual CubitStatus get_conic_mid_surface( Surface *surface_ptr1,
					       Surface *surface_ptr2,
					       BodySM *body_to_trim_to,
					       DLIList<BodySM*>& midsurface_bodies ) const = 0;

  virtual CubitStatus get_toric_mid_surface( Surface *surface_ptr1,
					       Surface *surface_ptr2,
					       BodySM *body_to_trim_to,
					       DLIList<BodySM*>& midsurface_bodies ) const = 0;

  virtual CubitStatus tweak_chamfer( DLIList<Curve*> &curve_list, 
                                     double left_offset,
                                     DLIList<BodySM*> &new_bodysm_list,
                                     double right_offset = -1.0,
                                     CubitBoolean keep_old_body = CUBIT_FALSE,
                                     CubitBoolean preview = CUBIT_FALSE )  = 0;
  /**<  Chamfer curves on solid bodies.  The left and right offsets are with 
    *   respect to the curve direction.  If the given right offset is negative,
    *   the left offset is used.  Users can preview to clarify the meaning of
    *   left and right.
    */

  virtual CubitStatus tweak_chamfer( DLIList<Point*> &point_list, 
                                     double offset1, 
                                     DLIList<BodySM*> &new_bodysm_list,
                                     Curve *edge1 = NULL,
                                     double offset2 = -1.0,
                                     Curve *edge2 = NULL,
                                     double offset3 = -1.0,
                                     Curve *edge3 = NULL,
                                     CubitBoolean keep_old_body = CUBIT_FALSE,
                                     CubitBoolean preview = CUBIT_FALSE )  = 0;
  /**<  Chamfer vertices on solid or sheet bodies.  On a solid body there can
    *   be up to 3 offsets; on a sheet body up to 2 offsets.  The offsets are
    *   in the direction of the supplied edges (for a solid) or the loop
    *   direction (for a sheet).  If multiple vertices are supplied, only one
    *   offset value is allowed and the edges are not used.
    */

  virtual CubitStatus tweak_fillet( DLIList<Curve*> &curve_list, 
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE )  = 0;
  /**<  Create a round fillet (or blend) at the given curves on solid bodies.
    */

  virtual CubitStatus tweak_fillet( Curve *curve_ptr, 
                                    double start_radius,
                                    double end_radius,
                                    BodySM *&new_bodysm_ptr,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) = 0;
  /**<  Create a round fillet (or blend) at the given curve on a solid body.
    *   The fillet has a variable radius from the start to the end of the curve.
    */

  virtual CubitStatus tweak_fillet( DLIList<Point*> &point_list, 
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) = 0;
  /**<  Create a round fillet (or blend) at the given vertices on sheet bodies.
    */

  virtual CubitStatus tweak_move( DLIList<Surface*> &surface_list,
                                  const CubitVector &delta,
                                  DLIList<BodySM*> &new_bodysm_list,
                                  CubitBoolean keep_old_body = CUBIT_FALSE,
                                  CubitBoolean preview = CUBIT_FALSE ) const = 0;
  /**<  Tweak specified faces of a volume or volumes along a vector.
    */

  virtual CubitStatus tweak_move( DLIList<Curve*> &curve_list,
                                  const CubitVector &delta,
                                  DLIList<BodySM*> &new_bodysm_list,
                                  CubitBoolean keep_old_body = CUBIT_FALSE,
                                  CubitBoolean preview = CUBIT_FALSE ) const = 0;
  /**<  Tweak specified curves of a sheet body along a vector.
    */

  virtual CubitStatus tweak_offset( DLIList<Surface*> &surface_list,
                                    double offset_distance,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const = 0;
  /**<  Tweak specified faces of a volume or volumes by offsetting those faces
    *   by the offset distance.
    */

  virtual CubitStatus tweak_offset( DLIList<Curve*> &curve_list,
                                    double offset_distance,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const = 0;
  /**<  Tweak specified curves of a sheet body or bodies by offsetting those
    *   curves by the offset distance.
    */

  virtual CubitStatus tweak_remove( DLIList<Surface*> &surface_list,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean extend_adjoining = CUBIT_TRUE,
                                    CubitBoolean keep_surface = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const = 0;
  /**<  Remove surfaces from a body or bodies and then extend the adjoining
    *   surfaces to fill the gap or remove the hole.
    */

  virtual CubitStatus tweak_remove( DLIList<Curve*> &curve_list,
                                    DLIList<BodySM*> &new_bodysm_list, 
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const = 0;
  /**<  Remove curves from a sheet body or bodies and then extend the remaining
    *   curves to fill the gap.  If an internal loop of curves is removed the
    *   hole is removed.
    */

  virtual CubitStatus tweak_target( DLIList<Surface*> &surface_list,
                                    DLIList<Surface*> &target_surf_list,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const = 0;
  /**<  Tweak specified faces of a volume or volumes up to target surfaces.
    */

  virtual CubitStatus tweak_target( DLIList<Curve*> &curve_list,
                                    DLIList<Surface*> &target_surf_list, 
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const = 0;
  /**<  Tweak specified edges of a surface or set of surfaces (in sheet
    *   bodies) up to a set of connected target surfaces.  This essentially
    *   extends or trims the attached surfaces of the sheet body.
    */

  virtual CubitStatus tweak_target( DLIList<Curve*> &curve_list,
                                    DLIList<Curve*> &target_curve_list, 
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg = CUBIT_FALSE,
                                    CubitBoolean keep_old_body = CUBIT_FALSE,
                                    CubitBoolean preview = CUBIT_FALSE ) const = 0;
  /**<  Tweak specified edges of a sheet body or bodies up to a set of target
    *   curves that are part of a sheet body.  The target is a surface created
    *   by thickening the owning surface of the target curve.
    */

  virtual CubitStatus remove_curve_slivers( BodySM *bodies, double lengthlimit ) const = 0;


  virtual CubitStatus create_net_surface( DLIList<Surface*>& ref_face_list, BodySM *& new_body,
                                          DLIList<DLIList<CubitVector*>*> &vec_lists_u, 
                                          DLIList<DLIList<CubitVector*>*> &vec_lists_v, 
                                          double net_tol = 1e-3, CubitBoolean heal = CUBIT_TRUE ) const = 0;

  virtual CubitStatus create_net_surface( DLIList<Curve*>& u_curves, DLIList<Curve*>& v_curves,
                                          BodySM *& new_body, double net_tol = 1e-3,
                                          CubitBoolean heal = CUBIT_TRUE ) const = 0;

  virtual CubitStatus create_offset_surface( Surface* ref_face_ptr, BodySM*& new_body, double offset_distance ) const = 0;

  virtual CubitStatus create_offset_body( BodySM* body_ptr, BodySM*& new_body, double offset_distance ) const = 0;

  virtual CubitStatus create_skin_surface( DLIList<Curve*>& curves, BodySM*& new_body ) const = 0;

  virtual CubitStatus loft_surfaces( Surface *face1, const double &takeoff1,
                                     Surface *face2, const double &takeoff2,
                                     BodySM*& new_body,
                                     CubitBoolean arc_length_option = CUBIT_FALSE,
                                     CubitBoolean twist_option = CUBIT_FALSE,
                                     CubitBoolean align_direction = CUBIT_TRUE,
                                     CubitBoolean perpendicular = CUBIT_TRUE,
                                     CubitBoolean simplify_option = CUBIT_FALSE) const = 0;

  virtual CubitStatus loft_surfaces_to_body( Surface *face1, const double &takeoff1,
									                  Surface *face2, const double &takeoff2,
									                  BodySM*& new_body,
									                  CubitBoolean arc_length_option,
									                  CubitBoolean twist_option,
									                  CubitBoolean align_direction,
									                  CubitBoolean perpendicular,
									                  CubitBoolean simplify_option) const = 0;
    
  virtual CubitStatus create_surface( DLIList<CubitVector*>& vec_list,
                                      BodySM *&new_body, Surface *ref_face_ptr, 
			                             CubitBoolean project_points ) const = 0;

  virtual CubitStatus create_weld_surface( CubitVector &root,
                                           Surface *ref_face1, double leg1, Surface *ref_face2, double leg2,
                                           BodySM *&new_body ) const = 0;

  protected: 
  
  private:

};

#endif

