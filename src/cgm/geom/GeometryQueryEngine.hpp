//-------------------------------------------------------------------------
// Filename      : GeoemtricQueryEngine.hpp
//
// Purpose       : Define the interface for all solid model query
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

#ifndef GEOMETRIC_QUERY_ENGINE_HPP
#define GEOMETRIC_QUERY_ENGINE_HPP

#include "CubitObserver.hpp"
#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitGeomConfigure.h"

class CubitString;
class CubitVector;
template <class X> class DLIList;
class Body;
class GMem;

class TopologyBridge;
class GeometryEntity;
class BodySM;
class Lump; 
class Surface;
class Curve;
class Point;

class CUBIT_GEOM_EXPORT GeometryQueryEngine
{

public:

  virtual ~GeometryQueryEngine() {}
    //- virtual destructor

      virtual CubitStatus export_solid_model(
                                   DLIList<TopologyBridge*>& bridge_list,
                                   const char* file_name,
                                   const char* file_type,
                                   const CubitString &cubit_version,
                                   const char* logfile_name = NULL ) = 0;

      virtual CubitStatus export_solid_model(
                                   DLIList<TopologyBridge*>& bridge_list,
				   char*& p_buffer,
				   int& n_buffer_size,
				   bool b_export_buffer) = 0;

     //R CubitStatus
     //R- CUBIT_SUCCESS/CUBIT_FAILURE
     //I ref_entity_list
     //I- A list of RefEntities to be exported or saved to a file.
     //I file_name
     //I- The name of the file to write to.
     //I file_type
     //I- An optional type of file.
     //- Export the current CUBIT geometry (everything in the Model) to a
     //- solid model format. Valid file types are:
     //-    "ACIS_SAT"    --  ACIS ASCII (SAT) file format
     //-    "ACIS_SAB"    --  ACIS BINARY (SAB) file format
     //-    "ACIS_DEBUG"  --  ACIS ASCII debug format
     //-    "IGES"        --  IGES file
     //-    "STEP"        --  STEP file
     //I logfile_name
     //I- Optional - name of logfile.
     //- No logfile gets created for SAB/SAT files, but for IGES and
     //- STEP file a logfile always gets created.  Default filenames
     //- are assigned if one is not given (iges_export.log, step_export.log).
     //-
     //- The function returns CUBIT_FAILURE if anything goes wrong with
     //- export - improper file type, inaccessible file, mismatch between
     //- the underlying representation and file type. It returns
     //- CUBIT_SUCCESS if everything goes well.

     virtual CubitStatus save_temp_geom_file(
                                 DLIList<TopologyBridge*> &ref_entity_list,
                                 const char *filename,
                                 const CubitString &cubit_version,
                                 CubitString &created_file,
                                 CubitString &created_file_type ) = 0;

     virtual CubitStatus import_temp_geom_file(
                                 FILE* file_ptr,
                                 const char* file_name,
                                 const char* file_type,
                                 DLIList<TopologyBridge*> &bridge_list  ) = 0;

     virtual CubitStatus import_solid_model(
                                const char* file_name,
                                const char* file_type,
                                DLIList<TopologyBridge*>& imported_entities,
                                bool print_results = true,
                                const char* logfile_name = NULL,
                                bool heal_step = true,
                                bool import_bodies = true,
                                bool import_surfaces = true,
                                bool import_curves = true,
                                bool import_vertices = true,
                                bool free_surfaces = true
					                      ) = 0;

     virtual CubitStatus import_solid_model(DLIList<TopologyBridge*> &imported_entities,
					    const char* pBuffer,
					    const int n_buffer_size) = 0;
       
     //R CubitStatus
     //R- CUBIT_SUCCESS/CUBIT_FAILURE
     //I file_ptr
     //I- A pointer to the file to read (can be NULL for IGES and STEP files).
     //I file_type
     //I- Type of file.
     //- Reads in geometry and creates the necessary Reference entities
     //- associated with the input geometry.
     //- Valid file types are:
     //-    "ACIS_SAT"    --  ACIS ASCII (SAT) file format
     //-    "ACIS_SAB"    --  ACIS BINARY (SAB) file format
     //-    "IGES"        --  IGES file
     //-    "STEP"        --  STEP file
     //O imported_entities
     //O- List of top-level entities read from file
     //I print_results
     //I- If false, fail silently (don't write error messages to stdout or stderr)
     //I heal_step - auto-healing of step bodies on import.  This is recommended
     //              because they always need it.
     //I import_bodies (etc...)
     //I- Should bodies be import.
     //- Function can selectively import solid bodies, free surfaces, free
     //- curves, or free vertices.  For example, the user may not want
     //- to import any free entities.
     //-
     //- The function returns CUBIT_FAILURE if anything goes wrong with
     //- import - improper file type, inaccessible file, mismatch between
     //- the underlying representation and file type. It returns
     //- CUBIT_SUCCESS if everything goes well.

      virtual CubitStatus get_underlying_curves(Curve * curve_ptr,
                                 DLIList<TopologyBridge*>& curve_list)  ;
      // -currently used only by VirtualGeometryEngine to get the underlying
      // - curves from virtual curves.
      virtual CubitStatus get_underlying_surfaces(Surface * surf_ptr,
                                 DLIList<TopologyBridge*>& surf_list)  ;

      virtual CubitStatus get_intersections(
                                    Curve* ref_edge1, CubitVector& point1,
                                    CubitVector& point2,
                                    DLIList<CubitVector*>& intersection_list,
                                    bool bounded = false,
                                    bool closest = false ) = 0;

      virtual CubitStatus get_intersections(
                                    Curve* ref_edge1, Curve* ref_edge2,
                                    DLIList<CubitVector*>& intersection_list,
                                    bool bounded = false,
                                    bool closest = false ) = 0;
      //- Finds the intersections of the two curves.  If the bounded flag is
      //- true, it finds only those intersections within the parameter range
      //- of both curves; otherwise it uses the extensions of the curves.  The
      //- closest option is currently valid only if both curves are straight,
      //- in which case it will return the 2 closest intersection locations,
      //- if the straight lines don't actually intersect. The function allocates
      //- the CubitVectors in the returned list, so be sure to free them.

      virtual CubitStatus get_intersections(
                                   Curve* curve,
                                   Surface* surface,
                                   DLIList<CubitVector*>& intersection_list,
                                   bool bounded = false ) = 0;
      //- Finds the intersections of the curve and surface.  The curve is extended
      //- if the bounded flag is not false.  The function allocates the CubitVectors
      //- in the returned list, so be sure to free them.

      virtual CubitStatus entity_extrema(
                                DLIList<GeometryEntity*> &ref_entity_list,
                                const CubitVector *dir1,
                                const CubitVector *dir2,
                                const CubitVector *dir3,
                                CubitVector &extrema,
                                GeometryEntity *&extrema_entity_ptr ) = 0;
      //- Gets the extrema position along the first given direction. If there
      //- is more than one extrema position, the other directions will be used
      //- to determine a unique position.  Directions 2 and 3 can be NULL.
      //- Entities supported include bodies, volumes, surfaces, curves and
      //- vertices.  The entity the extrema is found on is also returned.

      virtual CubitStatus entity_entity_distance(
                                         GeometryEntity *ref_entity_ptr1,
                                         GeometryEntity *ref_entity_ptr2,
                                         CubitVector &pos1,
                                         CubitVector &pos2,
                                         double &distance ) = 0;
      //- Gets the minimum distance between two entities and the closest positions
      //- on those entities. Supports vertices, curves, surfaces, volumes and bodies.

      virtual void delete_solid_model_entities(DLIList<BodySM*>& body_list) const = 0;
  virtual CubitStatus delete_solid_model_entities( BodySM* body_ptr ) const = 0;
  virtual CubitStatus delete_solid_model_entities(Surface* surf_ptr ) const = 0;
  virtual CubitStatus delete_solid_model_entities( Curve* curve_ptr ) const = 0;
  virtual CubitStatus delete_solid_model_entities( Point* point_ptr ) const = 0;
      //- Deletes the solid model entities associcated with the input
      //- free-floating RefEntity. If the input RefEntity is not free-floating,
      //- then its underlying ACIS ENTITYs are not deleted and CUBIT_FAILURE
      //- is returned.
      //- Returns CUBIT_SUCCESS if all went well, otherwise, CUBIT_FAILURE.

      virtual CubitStatus fire_ray(BodySM *body,
                                   const CubitVector &ray_start_point,
                                   const CubitVector &unit_direction,
                                   DLIList<double> &ray_params,
                                   DLIList<GeometryEntity*> *entity_list = NULL
                                   ) const = 0;
      //- fire a ray at the specified body, returning the entities hit and
      //- the parameters along the ray; return CUBIT_FAILURE if error

  virtual CubitStatus get_isoparametric_points(Surface* ref_face_ptr,
                                               int &nu, int &nv,
                                               GMem *&gMem) const = 0;

  virtual CubitStatus get_u_isoparametric_points(Surface* ref_face_ptr,
                                                 double v, int& n,
                                                 GMem *&gMem) const = 0;

  virtual CubitStatus get_v_isoparametric_points(Surface* ref_face_ptr,
                                                 double u, int&n,
                                                 GMem *&gMem) const = 0;

  virtual CubitStatus transform_vec_position(
    CubitVector const& position_vector,
    BodySM *OSME_ptr,
    CubitVector &transform_vector ) const = 0;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I position_vector on body, body that was transformed
    //O vector showning transform.
    //O-Computes the transform from the transformed body.

  virtual const char* modeler_type() = 0;
  virtual int get_major_version() = 0;
  virtual int get_minor_version() = 0;
  virtual int get_subminor_version() = 0;
  virtual int get_allint_version();
  virtual CubitString get_engine_version_string() = 0;

  virtual bool is_intermediate_engine() {return FALSE;}

    //- pass a string back identifying this query engine
  virtual CubitStatus set_export_version(const int major,
                                         const int minor);
  virtual CubitStatus set_export_allint_version(const int version);
  virtual CubitStatus list_engine_versions(CubitString &versions);
    //- gets, sets, and lists choices for solid modeling engine version(s)

  virtual double get_sme_resabs_tolerance() const = 0; // Gets solid modeler's resolution absolute tolerance
  virtual double set_sme_resabs_tolerance( double new_resabs ) = 0;

  virtual CubitStatus set_int_option( const char* opt_name, int val ) = 0;
  virtual CubitStatus set_dbl_option( const char* opt_name, double val ) = 0;
  virtual CubitStatus set_str_option( const char* opt_name, const char* val ) = 0;
   //- Set solid modeler options

  virtual CubitStatus get_graphics( Surface* surface_ptr,
                                    int& number_of_triangles,
                                    int& number_of_points,
                                    int& number_of_facets,
                                    GMem* gMem,
                                    unsigned short normal_tolerance=15,
                                    double distance_tolerance=0,
                                    double longest_edge = 0 ) const = 0;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I ref_face_ptr
    //I- The RefFAce for which hoops facetting information will be
    //I- gathered.
    //O- The number of polygons (triangles) needed for facetting.
    //O number_points
    //O- The number of points needed for facetting
    //O number_Facets
    //O- The number of facets needed for facetting.
    //O gMem
    //O- The sturage place for facets (and points).
    //= This function gathersw and outputs ACIS facet (and point)
    //- information for hoops involved in facetting an RefFace.  If
    //- all goes well, CUBIT_Success is retuned.

  virtual CubitStatus get_graphics( Curve* curve_ptr,
                                    int& num_points,
                                    GMem* gMem = NULL,
                                    double tolerance = 0.0 ) const = 0;
    //R CubitStatus
    //R- CUBIT_SUCCESS/CUBIT_FAILURE
    //I ref_edge_ptr
    //I- The RefEdge for which hoops facetting information will be
    //I- gathered/
    //O numSteps
    //O- The number of edges gathered.
    //O gMem
    //O- The sturage place for edges involved in facetting.
    //I tolerance
    //I- The tolerance deviation used when facetting the curve.
    //- This function tessellates the edge which is typically needed
    //- for graphically displaying the edge. If all goes well,
    //- CUBIT_SUCCESS is returned. Otherwise, CUBIT_FAILURE is returned.


  virtual CubitStatus translate( BodySM* body, const CubitVector& offset ) = 0;
  virtual CubitStatus rotate   ( BodySM* body, const CubitVector& axis, double angle ) = 0;
  virtual CubitStatus scale    ( BodySM* body, double factor ) = 0;
  virtual CubitStatus scale    ( BodySM* body, const CubitVector& factors ) = 0;
  virtual CubitStatus reflect  ( BodySM* body, const CubitVector& axis) = 0;
  virtual CubitStatus restore_transform( BodySM* body ) = 0;

  virtual CubitStatus translate( GeometryEntity* ent, const CubitVector& offset ) = 0;
  virtual CubitStatus rotate   ( GeometryEntity* ent, const CubitVector& axis, double degrees ) = 0;
  virtual CubitStatus scale    ( GeometryEntity* ent, double factor ) = 0;
  virtual CubitStatus scale    ( GeometryEntity* ent, const CubitVector& factors ) = 0;
  virtual CubitStatus reflect  ( GeometryEntity* ent, const CubitVector& axis) = 0;

  virtual CubitBoolean bodies_overlap (BodySM *body_ptr_1,
                                       BodySM *body_ptr_2 ) const = 0;

  virtual CubitBoolean volumes_overlap (Lump *lump1, Lump *lump2 ) const = 0; 

   protected:

   private:

} ;

#endif
