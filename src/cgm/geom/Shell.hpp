//-------------------------------------------------------------------------
// Filename      : Shell.hpp
//
// Purpose       : Represents a shell of a model.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef SHELL_HPP
#define SHELL_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "GroupingEntity.hpp"
#include "CubitDefines.h"
#include "CubitBox.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class ShellSM;
class RefFace;
class CoFace;
class RefVolume;
template <class X> class DLIList;
// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT Shell : public GroupingEntity
{
public :

  Shell() ;
    //- The default constructor.

  Shell(ShellSM* OSMEPtr) ;
    //- The constructor with a pointer to an other solid model entity.

  DagType dag_type() const { return DagType::shell_type(); }

  CoFace* get_co_face_ptr( RefFace* ref_face_ptr );
    //R CoFace*
    //R- A pointer to a child CoFace.
    //I ref_face_ptr
    //I- A pointer to the RefFace who's corresponding 
    //I- CoFace is to be found.
    //- If there is more than one CoFace in the shell for
    //- the specified RefFace, the first one encountered
    //- will be returned.
			
  RefVolume* get_ref_volume_ptr( );
    //R RefVolume*
    //R- A pointer to the parent RefVolume.

  CubitBox bounding_box();
    //R CubitBox
    //R- The bounding box of the Shell
    //- Calculate the union of the bounding boxes of the
    //- RefFaces of this Shell, and return it is the 
    //- bounding box of the Shell.
      
  ShellSM* get_shell_sm_ptr() const;
  
  CubitBoolean is_sheet();
    //- Return true if shell is a sheet.
      
protected: 

private:
  
    Shell( const Shell& );
    void operator=( const Shell&);
} ;

// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

