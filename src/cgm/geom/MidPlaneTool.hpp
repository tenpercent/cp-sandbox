//-------------------------------------------------------------------------
// Filename      : MidPlaneTool.hpp
//
// Purpose       :
//
// Special Notes : 
//
// Creator       :
//
// Creation Date :
//-------------------------------------------------------------------------

#ifndef MIDPOINTTOOL_HPP
#define MIDPOINTTOOL_HPP


#include "CubitDefines.h"
#include "CubitGeomConfigure.h"
class CubitVector;
class RefFace;
class Body;

template <class X> class DLIList;

class CUBIT_GEOM_EXPORT MidPlaneTool
{

  public:

  MidPlaneTool(){surfIndex=0;}
  ~MidPlaneTool(){}

  CubitStatus create_midplane( RefFace *ref_face1,
                               RefFace *ref_face2,
			       Body *body_ptr,
			       DLIList <RefFace*> &results_list );
  private:

  CubitStatus sort_surfaces(RefFace* big_face, RefFace* predecessor, DLIList <RefFace*> &mid_surfaces);

  int surfIndex;

  CubitStatus get_other_surfs(RefFace* base_face, RefFace* predecessor, DLIList <RefFace*> &other_surfs);

};

#endif

