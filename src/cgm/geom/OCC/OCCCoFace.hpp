//-------------------------------------------------------------------------
// Filename      : OCCCoFace.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Jane Hu 
//
// Creation Date : 04/23/08
//
//-------------------------------------------------------------------------

#ifndef OCCCOFACE_HPP
#define OCCCOFACE_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitEntity.hpp"
#include <TopoDS_Face.hxx>
#include "DLIList.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class GeometryQueryEngine;
class TopologyEntity;
class TopologyBridge;
class LoopSM;

class OCCBody;
class OCCLump;
class OCCShell;
class OCCSurface;
class OCCLoop;
class OCCCurve;
class OCCPoint;
// ********** END FORWARD DECLARATIONS     **********

class OCCCoFace
{
public:
  
  OCCCoFace(OCCSurface *surf_ptr, OCCShell *shell_ptr, CubitSense sense);
    //- A constructor
  
  virtual ~OCCCoFace() ;
    //- The destructor
    
  inline CubitSense sense(){return faceSense;}
    //- returns the sense of the underlying coface wrt the underlying face

  inline void set_sense(CubitSense sense) {faceSense = sense;}

  inline OCCSurface *surface()
    {return mySurface;}
    //- get the surface associated with this coedge
    
  inline void set_surface(OCCSurface* surface) {mySurface = surface;}

  inline OCCShell* shell() const { return myShell; }

  inline void set_shell(OCCShell * shell) {myShell = shell;}

  GeometryQueryEngine*
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.

  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );
protected: 
  
private:
  OCCSurface *mySurface;
  OCCShell *myShell;
  CubitSense faceSense;
};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

