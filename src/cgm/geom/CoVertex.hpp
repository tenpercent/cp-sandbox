//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : CoVertex.hpp
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

#ifndef COVERTEX_HPP
#define COVERTEX_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN MOTIF INCLUDES         **********
// ********** END MOTIF INCLUDES           **********

// ********** BEGIN OPEN INVENTOR INCLUDES **********
// ********** END OPEN INVENTOR INCLUDES   **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "CubitDefines.h"
#include "SenseEntity.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class RefVertex ;
class RefEdge ;
// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT CoVertex : public SenseEntity
{
   public :

      CoVertex() ;
      //- The default constructor

      CoVertex(RefVertex* vertexPtr) ;
      //I vertexPtr
      //I- The pointer to a vertex.
      //- The constructor.

      virtual ~CoVertex() ;
      //- The destructor

      DagType dag_type() const { return DagType::co_vertex_type(); }
  
#ifdef BOYD14
      RefEdge* get_ref_edge() ;
      //R RefEdge*
      //R- A pointer to the parent RefEdge
      //- This function returns a pointer to the parent RefEdge
      //- of this CoVertex.  It will assert that the CoVertex
      //- has exactly one parent RefEdge.
#endif

      RefVertex* get_ref_vertex_ptr()  ;
      //R RefVertex*
      //R- A pointer to the RefVertex which the current sense
      //R- entity is associated with.
      //- This function returns a pointer to the RefVertex which
      //- the current CoVertex is associated with.

   protected: 

   private:
    CoVertex( const CoVertex& );
    void operator=( const CoVertex&);
} ;


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

