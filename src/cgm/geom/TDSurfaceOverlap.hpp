//-Class: TDSurfaceOverlap.hpp

#ifndef TD_SURFACE_OVERLAP
#define TD_SURFACE_OVERLAP

#include "ToolData.hpp"
#include "CastTo.hpp"
#include "DLIList.hpp"
#include "RTree.hpp"
#include "AbstractTree.hpp"
#include "CubitGeomConfigure.h"

class RefFace;
class Body;
class SurfaceOverlapFacet;

class CUBIT_GEOM_EXPORT TDSurfaceOverlap : public ToolData
{
public:

  TDSurfaceOverlap( RefFace *ref_face_ptr, unsigned short ang_facet_tol, double abs_facet_tol,
    double gap_max );
  ~TDSurfaceOverlap();

  int has_rtree(){if(aTree) return CUBIT_TRUE; else return CUBIT_FALSE;}

  DLIList<SurfaceOverlapFacet*> *get_facet_list();
  AbstractTree<SurfaceOverlapFacet*> *get_facet_rtree();
  //- Get an AbstractTree containing the facets for this surface.  Note this returns
  //- a pointer to the AbstractTree stored in this class, so the calling code should
  //- not free it.

  DLIList<Body*> *get_body_list();

  static int is_surface_overlap(const ToolData* td)
     {return (CAST_TO(td, const TDSurfaceOverlap) != NULL);}

private:

  RefFace *refFacePtr;
  DLIList<SurfaceOverlapFacet*> facetList;
  AbstractTree<SurfaceOverlapFacet*> *aTree;
  DLIList<Body*> bodyList;
  int bodiesRetrieved;
  int marked;
  unsigned short angFacetTol;
  double absFacetTol;
  double gapMax;
};

#endif // TD_SURFACE_OVERLAP

