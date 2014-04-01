//-------------------------------------------------------------------------
// Filename      : CoEdge.hpp
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

#ifndef COEDGE_HPP
#define COEDGE_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "SenseEntity.hpp"
#include "ToolDataUser.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class RefEdge;
class RefFace;
class Loop;
class CoEdgeSM;
// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT  CoEdge : public SenseEntity,
               public ToolDataUser
{
   public :

      CoEdge() ;
      //- The default constructor

      virtual ~CoEdge() ;
      //- The destructor 

      CoEdge(RefEdge* edgePtr, CubitSense sense) ;
      //I edgePtr
      //I- The pointer to a edge.
      //I sense
      //I- The sense of this CoEdge. 
      //- The constructor with a pointer to a edge and the sense of this
      //- CoEdge.
      
      CoEdge(CoEdgeSM* OSMEPtr);

      DagType dag_type() const { return DagType::co_edge_type(); }

      RefEdge* get_ref_edge_ptr() ;
      //R RefEdge*
      //R- A pointer to the RefEdge which the current sense
      //R- entity is associated with.
      //- This function returns a pointer to the RefEdge which
      //- the current CoEdge is associated with.

      RefFace* get_ref_face();
      //R RefFace*
      //R- A pointer to the RefFace on which this co-edge is a part of.

      Loop* get_loop_ptr();
			//R Loop*
			//R- A pointer to the parent Loop of this CoEdge.
			
      virtual void switch_child_notify(ModelEntity const* newChild, 
                                       ModelEntity const* oldChild) ;
      //R void
      //I newChild
      //I- A pointer to the new child
      //I oldChild
      //I- A pointer to the old child
      //- This function is called after a child of a ModelEntity is
      //- switched. The sense of a CoEdge may change if one of its 
      //- RefEdges changes. This function takes care of that. If the sense
      //- of the RefEdges that were switched is same, nothing is done. If
      //- the RefEdges are of opposite sense, the sense of this object is
      //- switched, i.e. if it was FORWARD, it is made REVERSED, and vice
      //- versa. 

      CoEdgeSM* get_co_edge_sm_ptr() const;
      
      CubitStatus set_co_edge_sm_ptr(CoEdgeSM*);

      CubitBoolean about_spatially_equal( CoEdge* other_coedge,
                                          CubitSense relative_sense,
                                          double tolerance_factor = 1.0,
                                          CubitBoolean notify_refEntity = 
                                            CUBIT_FALSE );
                                            
      RefVertex* start_vertex();
      
      RefVertex* end_vertex();

   protected: 

   private:
    CoEdge( const CoEdge& );
    void operator=( const CoEdge&);
} ;


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

