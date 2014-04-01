//////////////////////////////////////////////////////////////////////
/// Filename      : MedialTool3D.hpp (interface for MedialTool3D)
///
/// Purpose       : 3-D Medial Axis extraction tool. Given a RefVolume
///                 generates the medial axis for that volume
///
/// Special Notes :
///
/// Creator       : Ragunath Sankaranarayanan (Ragu)
///
/// Creation Date : 08/23/2002
///
/// Owner         : Ragunath Sankaranarayanan (Ragu)
//////////////////////////////////////////////////////////////////////


#ifndef MEDIALTOOL3D_HDR
#define MEDIALTOOL3D_HDR

#include "DLIList.hpp"
#include "CubitGeomConfigure.h"
class RefVertex;
class RefEdge;
class RefFace;
class RefVolume;

//////////////////////////////////////////////////////////////////////
// Includes for when we really are compiling the medial stuff
//////////////////////////////////////////////////////////////////////
#ifdef USING_MEDIAL
//--------- Medial Extractor 3D Includes -----------------------------
#include "medial/medial_util/MedialDefines.hpp"
class Vec3D;
class Triangle3D;
class MedialVertex3D;
class MedialEdge3D;
class MedialFace3D;
//---------- Cubit Includes ------------------------------------------
class CubitVector;
#endif // USING_MEDIAL

class CUBIT_GEOM_EXPORT MedialTool3D  
{
public:

  /// constructor when a ref volume is given
  MedialTool3D(RefVolume *ref_volume_ptr);

  /// default destructor
  virtual ~MedialTool3D();

  /// create the medial axis of the volume given; The medial axis can
  /// be either a vertex, an edge, a face or a combination of the three;
  /// input parameters are cube_size, accuracy, min_angle_strong and
  /// min_angle_diff_components; refer "Medial3D.hpp" for desciption
  /// of the parameters.
  CubitStatus create_medial_axis_3d( DLIList <RefVertex*> &medial_verts,
                                     DLIList <RefEdge*> &medial_edges,
                                     DLIList <RefFace*> &medial_faces,
                                     bool full_search,
                                     double cube_size,
                                     double accuracy,
                                     double min_ang_strong,
                                     double min_angle_diff_components);
  
private:

  /// pointer to the RefVolume whose medial has to be extracted
  RefVolume *myRefVolume;

//////////////////////////////////////////////////////////////////////
// Functions that gets compiled only when compiling the medial stuff
//////////////////////////////////////////////////////////////////////
#ifdef USING_MEDIAL
  /// given a vector in Vec3D format, convert it to CubitVector format
  void convert_vec3d_to_cubitvector(Vec3D &my_vec3d, CubitVector &my_cubitvec);

  /// given a vector in CubitVector format, convert it to Vec3D format
  void convert_cubitvector_to_vec3d(CubitVector &my_cubitvec, Vec3D &my_vec3d);

  /// get all the boundary vertices and boundary facets of the RefVolume
  /// return them in Vec3D and Triangle3D format to be used by ObjectData3D.
  /// translation from CUBIT format to MedialExtractor format is also done
  CubitStatus get_boundary_facets( Vec3DPArray &boundary_verts,
                                   Triangle3DPArray &boundary_facets);

  /// draw the extracted medial to graphics output in the specified colors
  void draw_medial_results( MedialVertex3DPArray &output_vertices,
                            MedialEdge3DPArray &output_edges,
                            MedialFace3DPArray &output_faces,
                            int vert_color, int edge_color, int face_color);

	/// draw the medial face 3d as a shaded polygon
	void draw_medial_face(MedialFace3D* my_face, int tri_color);

#endif // USING_MEDIAL

};
#endif // MEDIALTOOL3D_HDR

