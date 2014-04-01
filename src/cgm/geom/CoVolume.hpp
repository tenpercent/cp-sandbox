//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : CoVolume.hpp
//
// Purpose       : A spurious class that represents the sense of a volume.
//                 We put this to keep a consistent hierarchy of classes.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef COVOLUME_HPP
#define COVOLUME_HPP

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

class RefVolume ;
class Body;

// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT CoVolume : public SenseEntity
{
   public :

      CoVolume() ;
      //- The default constructor

      CoVolume(RefVolume* volumePtr) ;
      //I volumePtr
      //I- The pointer to a volume.
      //- The constructor.

      virtual ~CoVolume() ;
      //- The destructor

      DagType dag_type() const { return DagType::co_volume_type(); }

      RefVolume* get_ref_volume_ptr() ;
      //R RefVolume*
      //R- A pointer to the RefVolume which the current sense
      //R- entity is associated with.
      //- This function returns a pointer to the RefVolume which
      //- the current CoVolume is associated with.
      
#ifdef BOYD14
      Body* get_body_ptr();
      //R Body*
      //R- A pointer to the Body which the current sense 
      //R- entity is associated with.
#endif

   protected: 

   private:
    CoVolume( const CoVolume& );
    void operator=( const CoVolume&);
} ;


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

