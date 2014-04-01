/**
 * \file GeometryQueryTool.hpp
 *
 * \brief This class provides the interface for all the query-only geometry
 *        related operations to the outside.
 *
 * This class is implemented using the Singleton pattern.
 *
 * \author Tim Tautges
 *
 * \date 2/2001
 *
 */

#ifndef GEOMETRYQUERYTOOL_HPP
#define GEOMETRYQUERYTOOL_HPP

#include <stdio.h>
#include <typeinfo>
#include <list>
#include <set>
#if !defined(NT)
using std::type_info;
#endif

#include "GeometryDefines.h"
#include "DLIList.hpp"
#include "GeometryQueryEngine.hpp"
#include "IntermediateGeomEngine.hpp"

class RefGroup;
class Body;
class BodySM;
class Lump;
class Shell;
class ShellSM;
class Loop;
class LoopSM;
class Chain;
class Surface;
class Curve;
class TopologyBridge;
class CoEdgeSM;

class GeometryEntity;


class TopologyEntity;
class TopologyBridge;

class RefEntity ;
class RefVolume ;
class RefFace ;
class RefEdge ;
class CoEdge;
class RefVertex ;

class SenseEntity;
class BasicTopologyEntity;
class RefVolume ;
class RefFace ;
class RefEdge ;
class RefVertex ;
class Point;
class CubitEntity;
class CoEdgeSM;

class GeometryQueryEngine;

class CubitPlane ;
class CubitBox;
class CubitVector;
class CubitTransformMatrix;
class TBOwner;

class CUBIT_GEOM_EXPORT GeometryQueryTool
{
public :

  friend class GeometryModifyTool;
   
  bool ige_is_composite(TBOwner *bridge_owner);
  bool ige_is_partition(TBOwner *bridge_owner);


  void ige_remove_modified(DLIList<TopologyBridge*>& geometry_list);

  static GeometryQueryTool* instance( GeometryQueryEngine* gqePtr = NULL);
  /**<
    *  \return GeometryQueryTool* - Pointer to the singleton GeometryQueryTool object
    *  \arg SMEPtr
    *   Pointer to a GeometryQueryEngine object. The default value
    *   is set to NULL.
    *
    *  Return a pointer to the only instance of the class with the default
    *  geometric query engine set to the argument, if it is not NULL.
    *  If the argument is NULL, return the pointer to the existing
    *  instance without modifying the default SME. A valid GQEPtr *must*
    *  be supplied at the time of the first call to this function.
    *  Hence, this instance function should specifically be called at
    *  startup with a valid non-NULL input GQEPtr.
    */

  ~GeometryQueryTool();
  ///<  Destructor.


  CubitStatus register_intermediate_engine( IntermediateGeomEngine* engine );

  void ige_remove_imprint_attributes_after_modify(DLIList<BodySM*> &old_sms,
                                              DLIList<BodySM*> &new_sms);
  void ige_push_imprint_attributes_before_modify
                                ( DLIList<BodySM*> &geom_list );
  void ige_export_geom( DLIList<TopologyBridge*> &geom_list );
  void ige_import_geom( DLIList<TopologyBridge*> &geom_list );
  void ige_remove_attributes( DLIList<TopologyBridge*> &geom_list );
  void ige_attribute_after_imprinting( DLIList<TopologyBridge*> &new_tbs,
                                                    DLIList<TopologyBridge*> &att_tbs,
                                                    DLIList<BodySM*> &new_sms,
                                                        DLIList<Body*> &old_bodies);
  void ige_remove_attributes_from_unmodifed_virtual(DLIList<TopologyBridge*> &bridges);

  /*! <HR><H1> global-list-functions Global entity list functions </H1>*/

  CubitStatus cubit_entity_list( const char* keyword,
                                 DLIList<CubitEntity*> &entity_list,
                                 const CubitBoolean print_errors = CUBIT_TRUE);
  /**<  return ref entities in a generic cubit entity list (overwrites list).
    *  returns CUBIT_FAILURE if keyword is not a ref entity name, and
    *  optionally prints error message
    */

  CubitStatus ref_entity_list( char const* keyword,
                               DLIList<RefEntity*> &entity_list,
                               const CubitBoolean print_errors = CUBIT_TRUE);
  /**<  return ref entities in entity_list (overwrites list).
    *  returns CUBIT_FAILURE if keyword is not a ref entity name,
    *  and optionally prints error message
    */

  CubitBox model_bounding_box();
  void bodies (DLIList<Body*> &bodies);
  void ref_volumes (DLIList<RefVolume*> &ref_volumes);
  void ref_groups (DLIList<RefGroup*> &ref_groups);
  void ref_faces (DLIList<RefFace*> &ref_faces);
  void ref_edges (DLIList<RefEdge*> &ref_edges);
  void ref_vertices (DLIList<RefVertex*> &ref_vertices);
  ///< Append global lists to arguments

  int num_bodies() const;
  int num_ref_volumes() const;
  int num_ref_groups() const;
  int num_ref_faces() const;
  int num_ref_edges() const;
  int num_ref_vertices() const;
  ///< Return number of entities of the specified type

  RefEntity *get_ref_entity (const char *type, int id);
  RefEntity *get_ref_entity (const type_info& type, int id);
  ///< Get a ref entity of the specified type and id

  Body *get_body ( int id );
  RefVolume *get_ref_volume ( int id );
  RefGroup *get_ref_group ( int id );
  RefFace *get_ref_face ( int id );
  RefEdge *get_ref_edge ( int id );
  RefVertex *get_ref_vertex ( int id );
  ///< Get an entity of the specified type and id

  Body *get_first_body ();
  RefVolume *get_first_ref_volume ();
  RefGroup *get_first_ref_group ();
  RefFace *get_first_ref_face ();
  RefEdge *get_first_ref_edge ();
  RefVertex *get_first_ref_vertex ();
  ///< Get the first entity in the global list of the specified type

  Body *get_next_body ();
  RefVolume *get_next_ref_volume ();
  RefGroup *get_next_ref_group ();
  RefFace *get_next_ref_face ();
  RefEdge *get_next_ref_edge ();
  RefVertex *get_next_ref_vertex ();
  ///< Get the next entity in the global list of the specified type

  Body *get_last_body ();
  RefVolume *get_last_ref_volume ();
  RefGroup *get_last_ref_group ();
  RefFace *get_last_ref_face ();
  RefEdge *get_last_ref_edge ();
  RefVertex *get_last_ref_vertex ();
  /**< Get the last entity in the global list of the specified type

    <HR><H3>File operations (import, export)</H3>
  */

  CubitStatus get_free_ref_entities(DLIList<RefEntity*> &free_entities);
    //- returns the list of free entities in the model

  void get_connected_free_ref_entities(
    RefEntity *entity,
    const int merge_option,
    DLIList<Body*> &body_list,
    DLIList<RefFace*> &ref_face_list,
    DLIList<RefEdge*> &ref_edge_list,
    DLIList<RefVertex*> &ref_vertex_list );
    //- get the free entities connected to, but not necessarily
    //- vertically related to, the entity.  - If merge_option, then take
    //- into account the fact that the model may be merged (it's slower
    //- that way).

  CubitStatus save_temp_geom_files( DLIList<RefEntity*>& ref_entity_list,
                                  const char* filename,
                                  const CubitString &cubit_version,
                                  std::list<CubitString> &files_written,
                                  std::list<CubitString> &types_written);

  CubitStatus export_solid_model( DLIList<RefEntity*>& ref_entity_list,
                                  const char* filename,
                                  const char * filetype,
                                  int &num_ents_exported,
                                  const CubitString &cubit_version,
                                  const char* logfile_name = NULL );
  /**<
   * Export entities to a solid model file.
    *  \arg ref_entity_list
    *  A list of RefEntities to be exported or saved to a file.
    *  \arg file_name
    *  The name of the file to write to.
    *  \arg file_type
    *  An optional type of file.
    *  \arg logfile_name
    *  Optional - name of logfile.
    *   \return  CubitStatus - success/failure
    *
    *  Export the current CUBIT geometry (everything in the Model) to a
    *  solid model format. Valid file types are:
    *     "ACIS_SAT"    --  ACIS ASCII (SAT) file format
    *     "ACIS_SAB"    --  ACIS BINARY (SAB) file format
    *     "ACIS_DEBUG"  --  ACIS ASCII debug format
    *     "IGES"        --  IGES file
    *     "STEP"        --  STEP file
    *  No logfile gets created for SAB/SAT files, but for IGES and
    *  STEP file a logfile always gets created.  Default filenames
    *  are assigned if one is not given (iges_export.log, step_export.log).
    *
    *  The function returns CUBIT_FAILURE if anything goes wrong with
    *  export - improper file type, inaccessible file, mismatch between
    *  the underlying representation and file type. It returns
    *  CUBIT_SUCCESS if everything goes well.
    *
    *  NOTE: if the ref_entity_list is empty, GeometryQueryTool gets the list of
    *  all entities in the current model, including free ref entities
    *
    */

  CubitStatus export_solid_model(DLIList<RefEntity*>& ref_entity_list,
				 char*& p_buffer,
				 int& n_buffer_size,
				 bool b_export_buffer);
  
  CubitStatus import_solid_model(const char* file_name,
                                 const char* file_type,
                                 const char* logfile_name = NULL,
                                 CubitBoolean heal_step = CUBIT_TRUE,
                                 CubitBoolean import_bodies = CUBIT_TRUE,
                                 CubitBoolean import_surfaces = CUBIT_TRUE,
                                 CubitBoolean import_curves = CUBIT_TRUE,
                                 CubitBoolean import_vertices = CUBIT_TRUE,
                                 CubitBoolean free_surfaces = CUBIT_TRUE,
				 DLIList<RefEntity*> *imported_entities = NULL);

  /**<
   * Import all or specified entities in a solid model file.
    *  \arg file_ptr
    *  A pointer to the file to read (can be NULL for IGES and STEP files).
    *  \arg file_type
    *  Type of file.
    *  \arg heal_step - auto-healing of step bodies on import.  This is recommended
    *     because they always need it.
    *  \arg import_bodies (etc...)
    *  Should bodies be import.
    *   \return  CubitStatus - success/failure
    *
    *  Reads in geometry and creates the necessary Reference entities
    *  associated with the input geometry.
    *  Valid file types are:
    *     "ACIS_SAT"    --  ACIS ASCII (SAT) file format
    *     "ACIS_SAB"    --  ACIS BINARY (SAB) file format
    *     "IGES"        --  IGES file
    *     "STEP"        --  STEP file
    *
    *  Function can selectively import solid bodies, free surfaces, free
    *  curves, or free vertices.  For example, the user may not want
    *  to import any free entities.
    *
    *  The function returns CUBIT_FAILURE if anything goes wrong with
    *  import - improper file type, inaccessible file, mismatch between
    *  the underlying representation and file type. It returns
    *  CUBIT_SUCCESS if everything goes well.
    */

  // import entities in a solid model buffer
  CubitStatus import_solid_model(DLIList<RefEntity*> *imported_entities,
				 const char* pBuffer,
				 const int n_buffer_size);

  CubitStatus construct_refentities(DLIList<TopologyBridge*> &topology_bridges,
                                    DLIList<RefEntity*> *imported_entities = NULL);
    //- given a list of TB's, construct ref entities for them; if the 2nd list pointer is
    //- non-NULL, pass back the list of ref entities in that list

  CubitStatus read_geometry_file(char const* fileName,
                                 char const* includePath = NULL,
                                 char const* type = "ACIS_SAT");
  /**<
   * Read geometry from the specified file.
    *  \arg filename
    *  The name of the file to read.
    *  \arg includePath
    *  An optional path to be used to search for the file if the file
    *  is not in the current directory.
    *  \arg type
    *  An optional type of file.
    *   \return  CubitStatus - success/failure
    *
    *  Reads a geometry file (named fileName).  If the file cannot be found
    *  in the current directory, and if includePath is not NULL, the include
    *  path specified by includePath is added to the file name, and it is
    *  searched in the given order. The file type specified in the "type"
    *  argument determines what type of geometry file it is -- the default
    *  is an ACIS SAT file.
    *  Supported geometry file types are:
    *     "ACIS_SAT"  --  ACIS SAT file format
    *     "ACIS_SAB"  --  ACIS SAB file format (binary)
    *  The function returns CUBIT_FAILURE if it cannot open the file, or
    *  the "type" is not known. Otherwise it returns CUBIT_SUCCESS.
    */

  int fire_ray(Body* body,
               const CubitVector ray_point,
               const CubitVector unit,
               DLIList<double>& ray_params,
               DLIList<RefEntity*> *entity_list);
  /**<
   * Fire a ray at a body, passing back number and distances of hits.
    *  \arg body
    *  Body at which you are firing rays
    *  \arg ray_point
    *  beginning point of ray
    *  \arg unit
    *  direction of ray (unit vector)
    * num_hit
    *-number of bodies hit
    * ray_params
    * array of parameters along ray at which bodies were hit
    * entity_list
    * entities hit by ray
    *  \return - error flag
    *
    *  fire a ray at the specified body, returning the entities
    *  hit and the parameters along the ray; return CUBIT_FAILURE
    *  if error
    */

  int fire_ray(RefFace *face,
               const CubitVector ray_point,
               const CubitVector unit,
               DLIList<double>& ray_params);
  
  static void geom_debug( DLIList<TopologyEntity*> );
  /**<  Temporary debug tool -- Steve Jankovich
    */
  static void set_facet_bbox( CubitBoolean pass_flag )
      {useFacetBBox = pass_flag;}
  static CubitBoolean get_facet_bbox()
      {return useFacetBBox;}
  /**< get and set the useFacetBBox flag. (compute our own bbox
    * instead of relying on bad ACIS geometry)
    */


  CubitString get_engine_version_string();
  /**<  Calls get_engine_version_string() function of underlying ModelingEngine.
    */

  CubitStatus set_export_allint_version(int version);
  /**<  set the major/minor version of the geometry engine
    */

  int get_allint_version();
  /**<  get the major/minor version of the geometry engine
    */

  CubitStatus list_engine_versions(CubitString &versions);
  /**<  returns a string with the engine version choices
    */

  double get_sme_resabs_tolerance();
  double set_sme_resabs_tolerance( double new_resabs );
  /**<  Gets/Sets solid modeler's resolution absolute tolerance
    */

  CubitStatus set_sme_int_option( const char* opt_name, int val );
  CubitStatus set_sme_dbl_option( const char* opt_name, double val );
  CubitStatus set_sme_str_option( const char* opt_name, const char* val );
  /**<  Set solid modeler options
    */

  ///< <HR><H3>Topology/geometry creation functions</H3>

  Body* make_Body(BodySM *bodysm_ptr) const;
  RefFace* make_RefFace(Surface* surface_ptr ) const;
  RefEdge* make_RefEdge(Curve* curve_ptr) const;
  RefVertex* make_RefVertex(Point* point_ptr) const;

  static CubitSense relative_sense( Surface* surface1, Surface* surface2 );

  RefFace* make_free_RefFace(Surface *surface_ptr ) const;
  RefEdge* make_free_RefEdge(Curve *curve_ptr ) const;
  RefVertex* make_free_RefVertex(Point *point_ptr) const;
  /**<  These functions can be used to create free ref-entities
    *  from a geometry engine.  Just call populate_topology_
    *  bridges to create the sm_ptr of the desired type, then pass
    *  that to the appropriate function here.
    */


  ///<HR><H3>Topology and geometry deletion</H3>

  void delete_Body( DLIList<Body*>& body_list );
  /**<  Deletes all Bodies in the input list from the model. Their
    *  associated Solid Model entities are deleted as well, if they exist,
    *  and if remove_solid_model_entities is CUBIT_TRUE.
    */

  CubitStatus delete_Body( Body* body_ptr );
  /**<  Deletes the input Body from the model.
    *  Its associated Solid Model entities are deleted as well, if
    *  they exist, and if remove_solid_model_entities is CUBIT_TRUE.
    *  If all went well, the input Body* reference will be set to NULL
    *  as the Body itself has been deleted.
    *  Returns CUBIT_SUCCESS if all went well, otherwise, CUBIT_FAILURE.
    */

  CubitStatus delete_single_Body( Body* body_ptr );
  /**< Behaves exactly as delete_Body, but in addition checks to see if
    *  children of Body are merged.  In some cases, 2 entities can be
    *  forced-merged, where they are not spatially equal.  This regenerates
    *  the graphics on merged entities so they look correct after a partner
    *  has been deleted.
    */

  CubitStatus delete_RefEntity( RefEntity* ref_entity_ptr );
  CubitStatus delete_RefFace( RefFace* ref_face_ptr );
  CubitStatus delete_RefEdge( RefEdge* ref_edge_ptr );
  CubitStatus delete_RefVertex( RefVertex* ref_vertex_ptr );
  /**<  This function is used to delete free-floating RefFaces, RefEdges
    *  or RefVertex'es.  All underlying VGI and solid model entities (if
    *  any) are also deleted. The entities will *not* be deleted is
    *  the input RefEntity has an owner (e.g., if the input RefEdge
    *  is owned by a RefFace).
    */

  void cleanout_deactivated_geometry();
  void cleanout_temporary_geometry ();

  void delete_geometry();

  ///< <HR><H3>Miscellaneous geometry evaluation functions</H3>

  CubitStatus interpolate_along_surface( CubitVector *vector_1,
                                         CubitVector *vector_2,
                                         DLIList<CubitVector*> &vector_list,
                                         RefFace* ref_face_ptr,
                                         int number_points ) const;
  /**<  Creates a list of vectors the length of the number_points that
    *  interpolate between vector_1 and vector_2.  All of the points
    *  will lie on the underlying equation of the refface.
    */

  CubitBoolean about_spatially_equal (const CubitVector& Vec1,
                                      const CubitVector& Vec2,
                                      double tolerance_factor = 1.0);
  /**<  \return  CubitBoolean
    *  \return - CUBIT_TRUE/CUBIT_FALSE
    *  \arg Vec1
    *  A reference to the first vector.
    *  \arg Vec2
    *  A reference to the second vector.
    *  \arg tolerance_factor
    *  Factor by which the absolute internal tolerance shall be
    *  multiplied.
    *
    *  Returns CUBIT_TRUE if the input Vec1 and Vec2 are spatially
    *  equivalent within a tolerance.  The internal spatial tolerance
    *  value is multiplied by tolerance_factor before the (spatial)
    *  test is done. Else, returns CUBIT_FALSE.
    */

  CubitBoolean about_spatially_equal (RefVertex* refVertex1,
                                      RefVertex* refVertex2,
                                      double tolerance_factor = 1.0);

  double geometric_angle(RefEdge* ref_edge_1,
                         RefEdge* ref_edge_2,
                         RefFace* ref_face );
  double geometric_angle(CoEdge* co_edge_1,
                         CoEdge* co_edge_2 );
  /**< calculates internal surface angles given 2 refedges on the
   * surface.  CoEdge version correctly handles curves in a surface
   * twice.
   */

  double surface_angle( RefFace *ref_face_1, RefFace *ref_face_2,
                        RefEdge *ref_edge = NULL,
                        RefVolume *ref_volume = NULL,
                        double frac = 0.5);
    /**< Calculate dihedral angle at ref_edge between two faces of the
     *   volume.
     */

  CubitStatus get_intersections( RefEdge* ref_edge1,
                                 CubitVector& point1,
                                 CubitVector& point2,
                                 DLIList<CubitVector*>& intersection_list,
                                 CubitBoolean bounded = CUBIT_FALSE,
                                 CubitBoolean closest = CUBIT_FALSE);

  CubitStatus get_intersections( RefEdge* ref_edge1, RefEdge* ref_edge2,
                                 DLIList<CubitVector*>& intersection_list,
                                 CubitBoolean bounded = CUBIT_FALSE,
                                 CubitBoolean closest = CUBIT_FALSE );
  /**<  Finds the intersections of the two curves.  If the bounded flag is
    *  true, it finds only those intersections within the parameter range
    *  of both curves; otherwise it uses the extensions of the curves.  The
    *  closest option is currently valid only if both curves are straight,
    *  in which case it will return the 2 closest intersection locations,
    *  if the straight lines don't actually intersect.  So far, other than
    *  for straight lines, this function only works if both curves are ACIS
    *  curves, unless both curves are linear.  The function allocates the
    *  CubitVectors in the returned list, so be sure to free them.
    */

  CubitStatus get_intersections( RefEdge* ref_edge, RefFace* ref_face,
                                 DLIList<CubitVector*>& intersection_list,
                                 CubitBoolean bounded = CUBIT_FALSE );
  /**< Finds the intersections of the curve and surface.  The curve is extended
    * if the bounded flag is not false.  The function allocates the CubitVectors
    * in the returned list, so be sure to free them.
    */

  CubitStatus entity_extrema( RefEntity *ref_entity_ptr,
                              const CubitVector *dir1,
                              const CubitVector *dir2,
                              const CubitVector *dir3,
                              CubitVector &extrema,
                              RefEntity *&extrema_entity_ptr );
  CubitStatus entity_extrema( DLIList<RefEntity*> &ref_entity_list,
                              const CubitVector *dir1,
                              const CubitVector *dir2,
                              const CubitVector *dir3,
                              CubitVector &extrema,
                              RefEntity *&extrema_entity_ptr );
  /** Gets the extrema position along the first given direction. If there
    * is more than one extrema position, the other directions will be used
    * to determine a unique position.  Directions 2 and 3 can be NULL.
    * Entities supported include bodies, volumes, surfaces, curves and
    * vertices.  The entity the extrema is found on is also returned.
    */

  CubitStatus entity_entity_distance( GeometryEntity *ge1,
                                           GeometryEntity *ge2,
                                           CubitVector &pos1, CubitVector &pos2,
                                           double &distance );
  CubitStatus entity_entity_distance( RefEntity *ref_entity_ptr1,
                                      RefEntity *ref_entity_ptr2,
                                      CubitVector &pos1, CubitVector &pos2,
                                      double &distance );
  /** Gets the minimum distance between two entities and the closest positions
    * on those entities. Supports vertices, curves, surfaces, volumes and bodies.
    */

  CubitBox bounding_box_of_bodies();
  /**<  returns the bounding box of all bodies in the model
    */

  ///< <HR><H3>Merging functions.</H3>

  static void set_geometry_factor( double fac );
  static double get_geometry_factor();
  static void set_merge_test_bbox(CubitBoolean tof);
  static CubitBoolean get_merge_test_bbox();
  static void set_merge_test_internal(int tof);
  static int get_merge_test_internal();
  static void set_sliver_curve_cleanup_tolerance( double tol );
  static void set_sliver_surface_cleanup_tolerance( double tol );
  static double get_sliver_curve_cleanup_tolerance(); 
  static double get_sliver_surface_cleanup_tolerance(); 

  //Initializes all settings of this class
  static void initialize_settings();

  static CubitStatus import_actuate( DLIList<RefEntity*> &entity_list );

  CubitBoolean same_query_engine(DLIList<TopologyEntity*> &topo_list) const;
  /**<  Returns CUBIT_TRUE if all the entities have the same geometric query engine and
    *  if that is the same one as the default.
    */

  GeometryQueryEngine* common_query_engine(
    DLIList<TopologyEntity*>& topology_list,
    DLIList<TopologyBridge*>& engine_bridges,
    CubitBoolean allow_default_engine
    = CUBIT_FALSE ) const;
  /**<  \return GeometryQueryEngine*
    *   A GeometryQueryEngine common at least one
    *   TopologyBridge of each of the passed TopologyEntities, or
    *   NULL if no common geometry engine is found.
    *   \arg topology_list
    *   The input list of TopologyEntities
    *   \arg engine_bridges
    *   Pass back the list of TopolgyBridges associated with each
    *   of the passed TopologyEntities (topology_list) and owned
    *   by the returned geometry engine.
    *   \arg allow_virtual_engine
    *   Return VirtualGeometryEngine::instance() if no common
    *   geometry enginge can be found.
    *
    *   Look for a common geometry engine other than the
    *   VirtualGeometryEngine.  If no common geometry engine other
    *   than VGE can be found and allow_virtual_engine is FALSE,
    *   NULL is returned.  If allow_virtual_engine is TRUE, and no
    *   common geometry engine is found, VGE will be returned, and
    *   engine_bridges will be populated with any virtual geometry
    *   if possible, otherwise with the first topology bridge attached
    *   to each of the passed TopologyEntities.
    */

  CubitBoolean does_geom_contain_query_engine(DLIList<TopologyEntity*> &topo_list,
                                              GeometryQueryEngine *engine) const;
  /**<  \return CubitBoolean
    *   Determine if any of the input entities contain the given query engine
    */

  CubitBoolean does_geom_contain_query_engine(DLIList<RefEntity*> &ref_entity_list,
                                              GeometryQueryEngine *engine,
                                              CubitBoolean children_too = CUBIT_FALSE) const;
  /**<  \return CubitBoolean
    *   Determine if any of the input entities contain the given query engine
    */

  TopologyEntity* entity_from_bridge( TopologyBridge* bridge_ptr ) const;

  void add_gqe(GeometryQueryEngine *gqe_ptr);
    /**< add a geometry query engine to the list
     */

  CubitStatus remove_gqe(GeometryQueryEngine *gqe_ptr);
    /**< remove a geometry query engine from the list; returns CUBIT_FAILURE
     *   if it wasn't on the list
     */

  void get_gqe_list(DLIList<GeometryQueryEngine*> &gqe_list);
    /**< return the list of gqe's
     */

  GeometryQueryEngine *get_gqe();
    /**< return the first gqe on the list
     */

  CubitStatus set_default_gqe(GeometryQueryEngine* gqe);
  /**< set the default GeometryQueryEngine
   */

  bool contains_intermediate_geometry(RefEntity*) const;
  bool contains_intermediate_geometry(DLIList<RefEntity*>& ref_entitylist) const;
  bool is_intermediate_geometry(RefEntity*) const;
  bool is_intermediate_geometry(TopologyBridge*) const;

//  GeometryQueryEngine *get_gqe(const EntityType gqe_type);
    /**< return the gqe of the specified type
     */

  CubitStatus destroy_dead_entity( TopologyEntity* topo_ent, bool top = true ) const;
    //- Remove this entity and any dead children where
    //- a dead child a) has no parent entities and b)
    //- has no topology_bridges.

  //Translate
  CubitStatus translate( Body* entity, const CubitVector& delta,
                         bool check_to_transform = true);
  CubitStatus translate( BasicTopologyEntity* entity, const CubitVector& delta,
                         bool check_to_transform = true);

  //Rotate
  CubitStatus rotate   ( Body* entity, const CubitVector& axis, double degrees,
                         bool check_to_transform = true);
  CubitStatus rotate   ( Body* entity,
                         const CubitVector& point,
                         const CubitVector& direction,
                         double degrees,
                         bool check_to_transform = true);
  CubitStatus rotate   ( BasicTopologyEntity* entity, const CubitVector& axis, double degrees,
                         bool check_to_transform = true);

  //Scale
  CubitStatus scale    ( Body* entity, double factor, bool check_to_transform = true);
  CubitStatus scale    ( Body* entity, const CubitVector& factors, bool check_to_transform = true);
  CubitStatus scale    ( BasicTopologyEntity* entity, double factor, bool check_to_transform = true);
  CubitStatus scale    ( BasicTopologyEntity* entity, const CubitVector& factors,
                         bool check_to_transform = true);

  //Reflect
  CubitStatus reflect  ( DLIList<Body*> bodies, const CubitVector& axis );
  CubitStatus reflect  ( BasicTopologyEntity* entity, const CubitVector& axis,
                         bool check_to_transform = true);
  CubitStatus restore_transform( Body* body );

  CubitBoolean volumes_overlap( RefVolume *volume_1, RefVolume *volume_2);
  CubitBoolean bodies_overlap( Body *body_ptr_1, Body *body_ptr_2 );

    //R CubitBoolean
    //R- CUBIT_TRUE if the two bodies overlap, CUBIT_FALSE if they don't
    //R- overlap.  If the bodies are touching the function
    //R- should return CUBIT_FALSE.
    //I body_ptr_1, body_ptr_2
    //I- The two body pointers that are being tested for overlap.
    //-  The function uses the intersect call to test if the bodies
    //-  are overlaping.  The full intersect Boolean is needed to see if
    //-  the bodies actually overlap and don't just touch.

protected :

  GeometryQueryTool(GeometryQueryEngine* GQEPtr);
  /**<  Constructor for the (singleton) GeometryQueryTool object
    */

private :

  CubitBoolean okay_to_transform( Body* body ) const;
  CubitBoolean okay_to_transform( BasicTopologyEntity* bte ) const;
  CubitStatus notify_observers_of_transform( RefEntity* ref_entity ) const;
  CubitStatus notify_intermediate_of_transform( TopologyEntity* bte,
                                             const CubitTransformMatrix& xform
                                             ) const;

  Shell* make_Shell(ShellSM *shellsm_ptr, bool& shell_modified ) const;
  RefVolume* make_RefVolume( Lump* lump_ptr, bool& vol_modified ) const;

  CubitStatus make_merged_RefFace( Surface* surface_ptr ) const;
    // Helper function for make_RefFace(Surface*)
    // For a merged RefFace, check if the passed surface topology
    // is still compatible.  If so, remerge any unmerged child
    // entities and return CUBIT_SUCCESS.  If not, return CUBIT_FAILURE.

  RefEntity *check_mergeable_refentity(GeometryEntity *bridge) const;
    // check for mergeable ref entity, indicated by a merge attribute on the
    // bridge

  CoEdgeSM *find_merged_coedgesm( Surface* on_this_surface,
                                  CoEdgeSM* merged_with_this ) const;
    // check for mergeable sense entity, which is a child of the entity pointed
    // to by a merge attribute on bridge, which corresponds to re_ptr

  CubitStatus straightline_intersections( RefEdge* ref_edge1,
                                 CubitVector & origin2,
                                 CubitVector & dir2,
                                 DLIList<CubitVector*>& intersection_list,
                                 CubitBoolean bounded = CUBIT_FALSE,
                                 CubitBoolean closest = CUBIT_FALSE);

  static GeometryQueryTool* instance_;
  /**<  static pointer to the unique instance of this class.
    */

  static CubitBoolean useFacetBBox;
  /**< For use in calculating a bounding box, you can do it based on
    *  a set of facets, rather than what default modeling engine  uses.
    */

  DLIList<GeometryQueryEngine*> gqeList;
  /**<  The list of geometry query engines
    */

  struct IGEComp : public std::binary_function<IntermediateGeomEngine*,
                                              IntermediateGeomEngine*,
                                              bool> {
    bool operator() ( const IntermediateGeomEngine* ptr1,
                      const IntermediateGeomEngine* ptr2 ) const
      { return ptr1->level() < ptr2->level(); }
  };
  typedef std::set<IntermediateGeomEngine*,IGEComp> IGESet;
  IGESet igeSet;

  GeometryQueryEngine *default_gqe;
  /**<  The default geometry query engine
    */

  static double curveSliverCleanUpTolerance;
  // After imprinting, an attempt at removing sliver curves is made.
  // Curves less than this tolerance will be removed.

  static double surfaceSliverCleanUpTolerance;
  // After imprinting, an attempt at removing sliver surfaces is made.
  // Removes sliver surfaces whose maximum gap distance among the long 
  // edges is smaller than the tolerance and who have at most three long edges. 

  static double geometryToleranceFactor;
  /**<  This factor is the the multiplier for the ACIS resabs
    *  when comparingcurves.  ALWAYS when using a multiplier
    *  use this factor for consistancy.
    */

  static CubitBoolean bboxMergeTest;
  static int internalSurfaceMergeTest; //0=off, 1=all, 2=splines only
  ///<  Options for refface merging.

  int maxPersistentBodyId;
  int maxPersistentRefVolumeId;
  int maxPersistentRefGroupId;
  int maxPersistentRefFaceId;
  int maxPersistentRefEdgeId;
  int maxPersistentRefVertexId;

  static const int CGM_MAJOR_VERSION;
  static const int CGM_MINOR_VERSION;
};


inline void GeometryQueryTool::set_geometry_factor( double factor )
{
  if ( factor < .0099999999999 )
    return;
  else
    geometryToleranceFactor = factor;
}
inline double GeometryQueryTool::get_geometry_factor()
{
  return geometryToleranceFactor;
}

inline void GeometryQueryTool::set_sliver_surface_cleanup_tolerance( double tol )
{
  surfaceSliverCleanUpTolerance = tol;
}

inline void GeometryQueryTool::set_sliver_curve_cleanup_tolerance( double tol )
{
  curveSliverCleanUpTolerance = tol;
}

inline double GeometryQueryTool::get_sliver_curve_cleanup_tolerance()
{
  return curveSliverCleanUpTolerance; 
}

inline double GeometryQueryTool::get_sliver_surface_cleanup_tolerance()
{
  return surfaceSliverCleanUpTolerance; 
}

inline CubitBoolean GeometryQueryTool::get_merge_test_bbox()
{return bboxMergeTest;}
inline int GeometryQueryTool::get_merge_test_internal()
{return internalSurfaceMergeTest;}
inline void GeometryQueryTool::set_merge_test_bbox(CubitBoolean tof)
{bboxMergeTest = tof;}
inline void GeometryQueryTool::set_merge_test_internal(int tof)
{internalSurfaceMergeTest = tof;}

inline void GeometryQueryTool::add_gqe(GeometryQueryEngine *gqe_ptr)
{
  assert(gqe_ptr != 0);
  if (!gqeList.move_to(gqe_ptr)) gqeList.append(gqe_ptr);
}
  /**< add a geometry query engine to the list
   */

inline CubitStatus GeometryQueryTool::remove_gqe(GeometryQueryEngine *gqe_ptr)
{
  assert(gqe_ptr != 0);
  CubitStatus status = CUBIT_FAILURE;
  if (gqeList.move_to(gqe_ptr)) {
    gqeList.remove();
    status = CUBIT_SUCCESS;
  }
  return status;
}

inline void GeometryQueryTool::get_gqe_list(DLIList<GeometryQueryEngine*> &gqe_list)
{
  gqe_list += gqeList;
}
  /**< return the list of gqe's
   */

inline GeometryQueryEngine *GeometryQueryTool::get_gqe()
{
  GeometryQueryEngine *gqe = NULL;
  if (gqeList.size()) {
    gqeList.reset();
    gqe = gqeList.get();
  }
  return gqe;
}

  /**< Return the VirtualGeometryEngine if it exists.  Otherwise return null.
   */

#endif



