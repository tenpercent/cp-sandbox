// TDSurfaceOverlap.cpp

#include "TDSurfaceOverlap.hpp"
#include "GMem.hpp"
#include "RefFace.hpp"
#include "SurfaceOverlapFacet.hpp"

TDSurfaceOverlap::TDSurfaceOverlap( RefFace *ref_face_ptr,
                                    unsigned short ang_facet_tol, double abs_facet_tol,
                                    double gap_max )
{
  refFacePtr = ref_face_ptr;
  aTree = NULL;
  marked = 0;
  bodiesRetrieved = 0;
  angFacetTol = ang_facet_tol;
  absFacetTol = abs_facet_tol;
  gapMax = gap_max;
}

TDSurfaceOverlap::~TDSurfaceOverlap()
{
  while( facetList.size() ) delete facetList.pop();
  if( aTree )
    delete aTree;
}

DLIList<SurfaceOverlapFacet*> *
TDSurfaceOverlap::get_facet_list()
{
  if( facetList.size() )
    return &facetList;

  GMem *gmem_ptr = new GMem;
  CubitStatus stat = refFacePtr->get_graphics( *gmem_ptr, angFacetTol, absFacetTol );

  if( !stat )
  {
    delete gmem_ptr;
    return NULL;
  }

  GPoint* plist = gmem_ptr->point_list();
  int* facet_list = gmem_ptr->facet_list();

  int i;
   
  GPoint p[3];
  for (i = 0; i < gmem_ptr->fListCount; )
  {
    int sides = facet_list[i++];
    if (sides != 3)
    {
      PRINT_WARNING("Skipping n-sided polygone in triangle list"
                    " in TDSurfaceOverlap.\n");
      i += sides;
    }
    else
    {
      p[0] = plist[facet_list[i++]];
      p[1] = plist[facet_list[i++]];
      p[2] = plist[facet_list[i++]];
   
      SurfaceOverlapFacet *facet = new SurfaceOverlapFacet( p );
      facetList.append( facet );
    }
  }

  delete gmem_ptr;

  return &facetList;
}

AbstractTree<SurfaceOverlapFacet*> *
TDSurfaceOverlap::get_facet_rtree()
{
  if( !aTree && !get_facet_list() )
    return NULL;

  if( !aTree )
  {
    aTree = new RTree<SurfaceOverlapFacet*>( gapMax );

    int i;
    for( i=facetList.size(); i--; )
    {
      SurfaceOverlapFacet *facet_ptr = facetList.get_and_step();
      aTree->add( facet_ptr );
    }
  }

  return aTree;
}

DLIList<Body*> *
TDSurfaceOverlap::get_body_list()
{
  if( bodiesRetrieved )
    return &bodyList;

  refFacePtr->bodies( bodyList );
  bodiesRetrieved = 1;
  return &bodyList;
}
