//-------------------------------------------------------------------------
// Filename      : CoFace.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef COFACE_H
#define COFACE_H

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "SenseEntity.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class RefFace;
class RefVolume;
class Shell;

// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT CoFace : public SenseEntity
{
   public :

      CoFace() ;
      //- The default constructor

      virtual ~CoFace() ;
      //- The destructor

      CoFace(RefFace* facePtr, CubitSense sense) ;
      //I facePtr
      //I- The pointer to a face.
      //I sense
      //I- The sense of this CoFace. 
      //- The constructor with a pointer to a face and the sense of this
      //- CoFace wrt its parent Face.

      DagType dag_type() const { return DagType::co_face_type(); }

      RefFace* get_ref_face_ptr() ;
      //R RefFace*
      //R- A pointer to the RefFace which the current sense
      //R- entity is associated with.
      //- This function returns a pointer to the RefFace which
      //- the current CoFace is associated with.

      RefVolume* get_ref_volume();
      //R RefVolume*
      //R- A pointer to the RefVolume that is associated with this CoFace.
      //- There is only one refVolume with this CoFace.  So if more than one
      //- are found this function will assert, if there are none found,
      //- NULL is returned.
			
			Shell* get_shell_ptr() ;
			//R Shell*
			//R- A pointer to the parent Shell of this CoFace, or NULL 
			//R- if there is no parent.

      virtual void switch_child_notify(ModelEntity const* newChild, 
                                       ModelEntity const* oldChild);
      //R void
      //I newChild
      //I- A pointer to the new child
      //I oldChild
      //I- A pointer to the old child
      //- This function is called after a child of a ModelEntity is
      //- switched. The sense of a CoFace may change if one of its 
      //- RefFaces changes. This function takes care of that. If the sense
      //- of the RefFaces that were switched is same, nothing is done. If
      //- the RefFaces are of opposite sense, the sense of this object is
      //- switched, i.e. if it was FORWARD, it is made REVERSED, and vice
      //- versa. 

      
   protected: 

   private:
    CoFace( const CoFace& );
    void operator=( const CoFace&);
} ;


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

