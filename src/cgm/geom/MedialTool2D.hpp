//////////////////////////////////////////////////////////////////////
/// Filename      : MedialTool2D.hpp
///
/// Purpose       : 2-D Medial Axis extraction tool. Uses the CMU code
///                 for creating medial. This class acts as the interface
///                 to that.
///
/// Special Notes :
///
/// Creator       : David White 
///
/// Creation Date : 06/21/2002
///
/// Modified by   : Ragunath Sankaranarayanan (09/11/2002)
///
/// Owner         : 
//////////////////////////////////////////////////////////////////////

#ifndef MEDIALTOOL2D_HPP
#define MEDIALTOOL2D_HPP

#include "DLIList.hpp"
#include "CubitGeomConfigure.h"

//---------------------------------------------------------
//Includes for when we really are compiling the medial stuff in.
#ifdef USING_MEDIAL
class MedialVertex;
class MedialSegment;
class Vec3D;

typedef DLIList <MedialVertex*> PointList;
typedef DLIList <MedialSegment*> SegList;
typedef DLIList <PointList*> PointLoopList;
typedef DLIList <SegList*> SegLoopList;
#include "medial/medial_util/MedialDefines.hpp"
class CubitVector;
#endif //USING_MEDIAL


class RefFace;
class RefEdge;

class CUBIT_GEOM_EXPORT MedialTool2D
{
private:
  RefFace *myRefFace;
#ifdef USING_MEDIAL

    /// draw the results to the screen
  void draw_results( MedialVertexPArray &output_points,
                     MedialSegmentPArray &output_segments);
  
  CubitStatus convert_to_segments( PointLoopList &boundary_point_loops,
         MedialSegmentPLoop &boundary_segments);
  
  CubitStatus get_boundary_points( RefFace *ref_face,
                                   PointLoopList &boundary_point_loops );

  CubitStatus get_curve_facets( RefEdge* curve, PointList &segments ) ;

    /// create geometry out of medial results. Generates a DLI list of
    /// RefVertex* and RefEdge*
  CubitStatus create_medial_geometry(MedialVertexPArray &output_points,
       MedialSegmentPArray  &output_segments, DLIList<RefEdge*> &medial_axis );

    /// given a vector in Vec3D format, convert it to CubitVector
  void convert_from_vec3d(Vec3D &vec, CubitVector &new_vec);

    /// given a vector in CubitVector format, convert it to Vec3D
  void convert_to_vec3d(CubitVector& vec, Vec3D &new_vec);

	  /// delete the memory for medial2D class
	void free_medial_memory(MedialVertexPArray &vertices,
													MedialSegmentPArray &segments);

#endif //USING_MEDIAL
  
public:
  MedialTool2D(RefFace* ref_face_ptr);
  ~MedialTool2D();

  CubitStatus create_medial_axis(DLIList <RefEdge*> &medial_axis,
                                 double cell_size,
                                 double angle_tol);
    //- Creates the medial axis of the surface and represents it
    //- by a set of curves.
  
};

#endif //MEDIALTOOL2D_HPP

