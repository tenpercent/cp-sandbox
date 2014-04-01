//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : BodySM.hpp
//
// Purpose       : To declare the Body solid model class.
//
// Special Notes :
//
// Creator       : Stephen J. Verzi
//
// Creation Date : 02/26/97
//
// Owner         : Stephen J. Verzi
//-------------------------------------------------------------------------

#ifndef BODYSM_HPP
#define BODYSM_HPP

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

#include "TopologyBridge.hpp"
class CubitTransformMatrix;
class CubitVector;
class CubitBox;

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT BodySM : public TopologyBridge
{
   public :

      BodySM() ;
      //- The default constructor

      virtual ~BodySM() ;
      //- The destructor
  
      virtual const type_info& topology_entity_type_info() const;
   

      virtual CubitStatus get_transforms( CubitTransformMatrix &tfm)= 0 ;
      //R CubitStatus
      //R- CUBIT_SUCCESS/FAILURE
      //- return the BODY transformation matrix


      virtual CubitStatus mass_properties( CubitVector& centroid,
                                           double& volume ) = 0;
      
      virtual CubitPointContainment point_containment( const CubitVector& pos ) = 0;

      CubitBox bounding_box();
      
   protected: 

   private:
} ;


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

