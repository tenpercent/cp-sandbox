//-------------------------------------------------------------
//Class: GeomMeasureTool
//Description:  Calculates measurements and statistics for a given entity.
//Author: David R. White
//Date: May 28, 2002
//-------------------------------------------------------------
#ifndef GEOM_MEASURE_TOOL_HPP
#define GEOM_MEASURE_TOOL_HPP
class RefEntity;
class RefVertex;
class RefEdge;
class RefFace;
class RefVolume;
class Body;
class RefGroup;
class TopologyEntity;
class CubitVector;
#include "DLIList.hpp"
template <class X> class DLIList;
class GeomPoint;
class GeomSeg;
#include "CubitGeomConfigure.h"

typedef DLIList <GeomPoint*> PointList;
typedef DLIList <PointList*> PointLoopList;
typedef DLIList <GeomSeg*> SegList;
typedef DLIList <SegList*> SegLoopList;

const double GEOM_END_LOWER = 87.0;
const double GEOM_END_UPPER = 93.0;
const double GEOM_SIDE_LOWER = 177.0;
const double GEOM_SIDE_UPPER = 183.0;
const double GEOM_CORNER_LOWER = 267.0;
const double GEOM_CORNER_UPPER = 273.0;
const double GEOM_REVERSAL_LOWER = 355.0;

class AreaHashTuple 
{
public:
  AreaHashTuple()
    {}
  ~AreaHashTuple()
    {}
  RefFace *myFace;
  double myArea;
};

class CUBIT_GEOM_EXPORT GeomMeasureTool
{
private:
  static int is_narrow_region_at_point(RefEdge *e1,
                                      const CubitVector &pt_on_e1,
                                      RefEdge *e2,
                                      const double &tol_sq);
  static int has_narrow_region(RefFace* face, const double &tol,
                                       const double &tol_sq);
  static void get_edges_from_list(DLIList <RefEntity*> &entity_list,
                                  DLIList <RefEdge*> &ref_edges );
  static void get_faces_from_list(DLIList <RefEntity*> &entity_list,
                                  DLIList <RefFace*> &ref_faces );
  static void get_volumes_from_list(DLIList <RefEntity*> &entity_list,
                                    DLIList <RefVolume*> &ref_volumes );
  static void get_bodies_from_list(DLIList <RefVolume*> &entity_list,
                                   DLIList <Body*> &ref_bodies );
    ///
    ///converts the entity list, and expands the list to get the edges,
    ///faces or volumes of the entities.
    ///

  static double dist_sq_point_data( CubitVector &curr_point,
                                    GeomSeg *&curr_seg );
    ///
    /// finds the distance between the point and the segment.  Function
    /// used in the AbstractTree Nearest Neighbor calculation.

  static CubitStatus interior_angle(GeomSeg *first_seg,
                                    GeomSeg *next_seg,
                                    RefFace *face,
                                    double &angle );
    ///
    /// measures the interior angle between the two connected segments.
    /// Returns an error if the segments are not connected or if the
    /// surface normal length is zero.
    ///

   static void find_index(DLIList<AreaHashTuple*> *hash_table,
                          int table_size,
                          RefFace *ref_face,
                          AreaHashTuple *&curr_tuple );
    ///
    /// Finds the tuple for the given ref_face in the hash_table.
    /// This function is used in find_adjacent_face_ratios.
    ///

   static void get_adjacent_faces(RefFace* curr_face,
                                  DLIList<RefFace*> &adjacent_faces);
    ///
    /// Finds the adjacent faces to this face by traversing the DAG.
    ///

   
   static CubitBoolean find_opposite_edge( RefEdge* ref_edge,
                                           RefFace* ref_face,
                                           RefEdge *&other_edge,
                                           CubitVector &closest_point);
    ///
    /// finds the other edge which is about opposite to ref_edge.  If
    /// there are more than two loops, we return CUBIT_FALSE.  If
    /// there are two loops it finds the loop and curve on that loops,
    /// that has the closest point  to the mid point of ref_edge.
    /// If there is just one loop, it turns on that loop to find the
    /// opposite.  There can be no convex vertex angles on the loop.
    ///
   
   static CubitBoolean is_equal(double v1, double v2);
    ///
    /// Tests to see if the two double values are within CUBIT_RESABS.
    /// Returns TRUE if they are, false if they aren't.
    ///

 public:
  GeomMeasureTool()
    {}
  
  ~GeomMeasureTool()
    {}

  static void measure_curve_length(DLIList <RefEntity*> &entity_list,
                                   double &smallest,
                                   RefEdge *&smallest_edge,
                                   double &largest,
                                   RefEdge *&largest_edge,
                                   double &average,
                                   double &sum);
    ///
    ///Find the smallest, longest and average size of the curves.  Also
    ///compute the sum of all the curve.  Do that for a list of RefEntity pointers.
    ///

  static void measure_curve_length(DLIList <RefEdge*> &ref_edges,
                                   double &smallest,
                                   RefEdge *&smallest_edge,
                                   double &largest,
                                   RefEdge *&largest_edge,
                                   double &average,
                                   double &sum);
    ///
    ///Find the smallest, longest and average size of the curves.  Also
    ///compute the sum of all the curve.  Do that for a list of RefEdge pointers.
    ///
    

  static void measure_face_curves( RefFace *ref_face,
                                   double &min_curve_length,
                                   double &max_curve_ratio,
                                   RefEdge *&min_ref_edge);
    ///
    /// Finds the curve with the smallest length on the surface.
    /// Also find the biggest change in curve lengths where
    /// the ratio is the largest curve divided by the smallest curve.
    ///

  static CubitStatus measure_face_loops(RefFace *face,
                                        double &min_distance_between_loops,
                                        double &min_distance_in_one_loop,
                                        double &min_angle, double &max_angle,
                                        double tolerance);
    ///
    ///Find the smallest distance between multiple loops on the surface.
    ///Find the smallest distance between two curves.  This could be the length
    ///of the smallest curve...
    ///Also calculate the angles between curves on the faces.
    ///The tolerance is used to specify the faceting to approximate the boundaries.
    ///

  static void  measure_face_area (DLIList <RefEntity*> &entity_list,
                                  double &smallest,
                                  RefFace *&smallest_face,
                                  double &largest,
                                  RefFace *&largest_face,
                                  double &average,
                                  double &sum);
    ///
    ///Find the smallest, largest and average area of the faces contained
    ///in the entity list..  Also compute the sum of all the face area.
    ///

  static void  measure_face_area (DLIList <RefFace*> &ref_faces,
                                  double &smallest,
                                  RefFace *&smallest_face,
                                  double &largest,
                                  RefFace *&largest_face,
                                  double &average,
                                  double &sum);
    ///
    ///Find the smallest, largest and average area of the faces 
    ///in the face list..  Also compute the sum of all the face area.
    ///
  
  static void  measure_volume_volume (DLIList <RefEntity*> &entity_list,
                                      double &smallest,
                                      RefVolume *&smallest_volume,
                                      double &largest,
                                      RefVolume *&largest_volume,
                                      double &average,
                                      double &sum);
    ///
    ///Find the smallest, largest, and average volume.  Compute the sum of all
    ///the volumes for the volumes contained by the entities in the entity list.
    ///

  static void  measure_volume_volume (DLIList <RefVolume*> &ref_volumes,
                                      double &smallest,
                                      RefVolume *&smallest_volume,
                                      double &largest,
                                      RefVolume *&largest_volume,
                                      double &average,
                                      double &sum);
    ///
    ///Find the smallest, largest, and average volume.  Compute the sum of all
    ///the volumes in the list.
    ///

#ifdef BOYD14
  static void  measure_face_hydraulic_radius( DLIList <RefEntity*> &entity_list,
                                              double &smallest,
                                              RefFace *&smallest_face,
                                              double &largest,
                                              RefFace *&largest_face,
                                              double &average,
                                              double &sum );

  static void  measure_face_hydraulic_radius( DLIList <RefFace*> &ref_faces,
                                              double &smallest,
                                              RefFace *&smallest_face,
                                              double &largest,
                                              RefFace *&largest_face,
                                              double &average,
                                              double &sum );
#endif

#ifdef BOYD14
  static void measure_volume_hydraulic_radius(DLIList <RefEntity*> &entity_list,
                                              double &smallest,
                                              RefVolume *&smallest_volume,
                                              double &largest,
                                              RefVolume *&largest_volume,
                                              double &average,
                                              double &sum );

  static void measure_volume_hydraulic_radius(DLIList <RefVolume*> &ref_volumes,
                                              double &smallest,
                                              RefVolume *&smallest_volume,
                                              double &largest,
                                              RefVolume *&largest_volume,
                                              double &average,
                                              double &sum );
#endif

  static void measure_face_area_and_hydraulic_radius(RefFace *curr_face,
                                                     double &face_area,
                                                     double &face_hydraulic_radius);

  static void measure_volume_volume_and_hydraulic_radius(RefVolume *curr_volume,
                                                         double &volume_volume,
                                                         double &volume_hydraulic_radius);
  
  static void angles_between_volume_surfaces(RefVolume *curr_volume,
                                             double &min_angle,
                                             double &max_angle,
                                             RefFace *&face_min_1,
                                             RefFace *&face_min_2,
                                             RefFace *&face_max_1,
                                             RefFace *&face_max_2);

#ifdef BOYD14
  static void curvature_info_for_face(RefFace *curr_face,
                                      double &max_curvature_change,
                                      double &min_curvature_change,
                                      double &average_curvature_change);
#endif

  static void merged_unmerged_surface_ratio(DLIList <RefVolume*> &ref_volumes,
                                            int &merged, int &unmerged,
                                            double &ratio);

  static void report_intersected_volumes(DLIList <RefVolume*> &volume_list,
                                        DLIList <RefVolume*> &intersection_list);
  
  static void report_intersected_bodies(DLIList <Body*> &ref_bodies,
                                        DLIList <Body*> &intersection_list);

  static RefFace* valid_start( DLIList <RefFace*> &all_faces );

  static void face_list_from_volume_list( DLIList <RefVolume*> &input_vols,
                                          DLIList <RefFace*> &all_faces );
  
  static void find_shells( DLIList <RefVolume*> &input_vols,
                           RefGroup *&owner_groups,
                           int &number_of_shells);

  static void ratio_of_shells_to_volumes(int number_of_shells,
                                         DLIList <RefVolume*> &ref_volumes,
                                         int &number_of_volumes,
                                         double &ratio);
  
  
  static CubitStatus get_boundary_points( RefFace *ref_face,
                                          PointLoopList &boundary_point_loops,
                                          double seg_length_tol);
    ///
    ///Gets the points in ordered loops that are used to facet the boundary
    ///of the surface.
    ///

  static CubitStatus get_curve_facets( RefEdge* curve,
                                       PointList &segments,
                                       double seg_length_tol );
    ///
    ///Get the points that approximate the curve, use the graphics facets.
    ///  

  static CubitStatus convert_to_lines(PointLoopList &boundary_point_loops,
                                      SegLoopList &boundary_line_loops,
                                      RefFace *ref_face );
    ///
    ///Converts the loops of points to line segments.  Can fail if it
    ///can't determine the owner of the line segment. The refface may
    ///be needed to help determine ownership.
    ///

  static void print_surface_measure_summary( DLIList <RefFace*> &ref_faces );
    ///
    ///Prints a summary of all of the surface information.
    ///min curve length, max adjacent curve ratios, min angle, max angle,
    ///min area, max area, min hydraulic radius, min distance between loops,
    ///min distance between a single loop.
    ///

  static void print_volume_measure_summary(DLIList <RefVolume*> &ref_volumes);
    ///
    /// Prints summary of all the volume information, including the
    /// information of the
    /// surfaces in the volumes.
    ///

#ifdef BOYD14
  static void print_body_measure_summary(DLIList <Body*> &bodies);
    ///
    /// Prints summary of all the body information.  volume and surface
    /// info is also printed.
    ///
#endif

  static void find_adjacent_face_ratios(RefVolume *curr_volume,
                                        double &max_face_ratio,
                                        RefFace *&big_face,
                                        RefFace *&small_face);
    ///
    /// Finds the max ratio of adjacent face areas (big face area
    /// divided by  small face area);
    ///

  static void find_small_curves( DLIList <RefVolume*> &ref_vols,
                                 double tol,
                                 DLIList <RefEdge*> &small_curves,
                                 DLIList <double> &small_lengths);
    ///
    /// Finds the curves with lengths less than the given tol.
    ///

  static void find_surfs_with_narrow_regions( DLIList <RefVolume*> &ref_vols,
                                 double tol,
                                 DLIList <RefFace*> &surfs_with_narrow_regions);
    ///
    /// Finds the curves with lengths less than the given tol.
    ///

  static void find_small_faces( DLIList <RefVolume*> &ref_vols,
                                double tol,
                                DLIList <RefFace*> &small_faces);
    ///
    /// Finds the faces with areas less than the given ammount.
    ///

  static void find_small_faces_hydraulic_radius( DLIList <RefVolume*> &ref_vols,
                                                 double tol,
                                                 DLIList <RefFace*> &small_faces,
                                                 DLIList <double> &small_hyd_rad,
                                                 double &percent_planar,
                                                 double &percent_pl_co);
    ///
    /// Finds the faces with hydraulic radii less than tol.
    /// The hydraulic radus is defined by 4*(A/P) where
    /// A is the area of the surface and P is the total perimiter
    /// length around the surface.
    /// Also finds the % surfaces that are planar.   And finds
    /// the % surfaces that are planar and conical.
    /// This is based on the geometry type, not some actual measurement.
    ///
  
  static void find_small_volumes( DLIList <RefVolume*> &ref_vols,
                                  double tol,
                                  DLIList <RefVolume*> &small_volumes);
    ///
    /// Finds the volumes with volumes less than tol.
    ///

  static void find_small_volumes_hydraulic_radius( DLIList <RefVolume*> &ref_vols,
                                                   double tol,
                                                   DLIList <RefVolume*> &small_volumes,
                                                   DLIList <double> &small_hyd_rad);
    ///
    /// Measures the hydraulic radii for the volumes by 6*(V/A) where
    /// V is the volume of the volume and A is the total area of all
    /// the volume's surfaces.  Volumes with small hydraulic radii (compared
    /// to tol) are stored in the small_volumes list.
    ///
  
  static void find_interior_curve_angles( RefVolume* ref_volume,
                                          double upper_bound,
                                          double lower_bound,
                                          DLIList <RefEdge*> &large_edge_angles,
                                          DLIList <RefEdge*> &small_edge_angles,
                                          DLIList <double> &large_angles,
                                          DLIList <double> &small_angles,
                                          int &total_interior,
                                          int &total_fuzzy);
  
    ///
    /// For each surface in the volume, find the interior angles
    /// defined by the curves in the surfaces that are either below or
    /// greater than the the upper and lower bounds.  The angles
    /// and therefor bounds, are tracked in degrees.
    /// Curves that fit into these categories are stored pairwise in the
    /// large or small edge_angles lists.  The corresponding angle
    /// measurements are stored in the respective small and large angles lists.
    ///

  static void find_sharp_tangential_meets( RefVolume *ref_volume,
                                           DLIList <RefEdge*> &small_angle_edge_pairs,
                                           DLIList <RefFace*> &tangential_surface_pairs );
    ///
    /// Find the tangential meetings in the volume.
    /// This specifically looks for surfaces that meet tangentially
    /// that would be a problem.  Usually these are surfaces that
    /// come into a side 180 degrees and on top there is a
    /// sharpe angle.  Usually if there isn't a sharpe curve
    /// angle these meetings are not bad.
    /// Note that this function assumes that you are passing
    /// in sets of curve edges that have small interior angles between them.
    /// It also assumes that the edge pairs are ordered as they would
    /// be found in the surfaces (first edge then next edge).
    /// Basically, call the funciton, find_interior_curve_angles
    /// before calling this function, and pass this function those results.
    ///

  static void find_dihedral_angles( DLIList<RefVolume*> &ref_vols,
                                    double upper_bound,
                                    double lower_bound,
                                    DLIList <RefFace*> &large_face_angles,
                                    DLIList <RefFace*> &small_face_angles,
                                    DLIList <double> &large_angles,
                                    DLIList <double> &small_angles,
                                    int &total_interior,
                                    int &total_fuzzy,
                                    int &total_not_flat);
    ///
    /// Finds the large and small angles or rather the surfaces on the volumes
    /// that make them.
    ///
  

  static void find_dihedral_angles( RefVolume *curr_volume,
                                    double lower_bound,
                                    double upper_bound,
                                    DLIList <RefFace*> &large_face_angles,
                                    DLIList <RefFace*> &small_face_angles,
                                    DLIList <double> &large_angles,
                                    DLIList <double> &small_angles,
                                    int &total_interior,
                                    int &total_fuzzy,
                                    int &total_not_flat);
    ///             
    /// Finds the surfaces that have big and small angles compared to the upper
    /// and lower bounds (these should be in degrees...).
    /// The faces that make these angles are paired (one after the other)
    /// in the lists large and small angles.
    ///

  static void find_close_loops(DLIList <RefVolume*> &ref_vols,
                               DLIList <RefEdge*> &close_edges,
                               DLIList <RefFace*> &close_loop_faces,
                               DLIList <double> &small_lengths,
                               double tol);
    ///
    /// Finds the surfaces with multiple loops that have small, relative
    /// to tol, distances between the loops.
    /// Note that all the lists are syncronized to be connected.
    /// The pairs of close edges with the small_lengths and
    /// the surfaces that they are on.  The close_loop_faces list
    /// may not contain distinct instances of faces, as they
    /// may be repeated.  Note that this will only find the
    /// closest edge pair, within tol, for each loop.
    ///

  static void find_close_loops(RefFace *face,
                               DLIList <RefEdge*> &close_edges,
                               DLIList <double> &small_lengths,
                               double tol);
    ///
    /// Finds the edges of the loops that are within some tol.  The
    /// edges are stored in successive pairs in close_edges while
    /// the actual distances are stored in the small_lengths list.
    /// These lists are syncronized (and are assumed to be empty upon
    /// entrance to the function).
    ///

#ifdef BOYD14
  static void count_entities(DLIList<RefVolume*> &ref_vols,
                             int &num_faces,
                             int &num_curves,
                             int &num_vertices);
    ///
    /// Find the total number of entities in this volume list.
    /// Makes sure not to double count.
    ///
#endif

  static double measure_area(RefFace* curr_face);
    ///
    /// Measures the area of the curr_face.  Uses
    /// the GeomDataObserver class to cache the area on the face.
    /// After calling this function, curr_face will have
    /// a GeomDataObserver observing it.  It will get removed if
    /// the surface is altered or changed.
    ///

  static void find_bad_geometry(RefVolume *volume,
                                DLIList <RefEntity*> &bad_ents);
    ///
    /// This function simply gets the bad entities of the volume.  This
    /// assumes the volume is from some solid modelar where this function
    /// is defined. If not, it will just be empty...
    ///

  static void find_irregular_valence( DLIList <RefVolume*> &ref_volumes,
                                      DLIList <RefVertex*> &irregular_vertices);
  static void find_irregular_valence( RefVolume* ref_volume,
                                      DLIList <RefVertex*> &irregular_vertices);
    ///
    /// Find the irregular vertex valences.
    /// Find things like vertices with valences greater than 4.
    /// Assume for now that the volumes are not merged..
    ///

  static void find_blends( RefVolume *ref_volume,
                           DLIList <RefFace*> &blend_faces,
                           DLIList<DLIList<RefFace*>*> &blend_groups);
    ///   
    ///  Find fillets and rounds.
    ///
    ///  ====IMPORTANT====
    ///  Calling functions BEWARE, YOU MUST DELETE THE LISTS OF
    ///  BLEND SURFACES IN THE BLEND GROUPS AFTER CALLING THIS FUNCTION!!!!
    ///  ====IMPORTANT====
    ///
    ///

  static CubitBoolean is_face_blend(RefFace *ref_face,
                                    RefVolume *ref_volume,
                                    RefEdge *&ref_edge,
                                    RefEdge *&other_edge);
    ///
    /// Determines if a face is a blend surface, returns the
    /// ref_edge on one side of the blend and other_edge on
    /// the opposite side.  For this type of blend, only ref_edge 
    /// must be tangentially meeting with another surface.  
    /// Other_edge must be oriented orthogonally to it and may 
    /// or may not blend with another surface.  This assumes a
    /// rectangular blend surface, without holes.
    ///
 
  static CubitStatus get_centroid( RefFace *ref_face, CubitVector &centroid, double &tot_area );
	///
    /// Takes the area-weighted average of all display facet
    /// centroids and updates the passed in CubitVector to be 
    /// that location.  Also updates the tot_area variable to 
    /// give the total facet area for the passed in ref_face.
    ///

  static CubitStatus center( DLIList<RefFace*> ref_faces );
};
#endif

