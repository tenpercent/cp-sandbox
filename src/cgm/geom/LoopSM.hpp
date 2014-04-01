//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : LoopSM.hpp
//
// Purpose       : To declare the Loop solid model class.
//
// Special Notes :
//
// Creator       : Stephen J. Verzi
//
// Creation Date : 02/26/97
//
// Owner         : Stephen J. Verzi
//-------------------------------------------------------------------------

#ifndef LOOPSM_HPP
#define LOOPSM_HPP

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

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT LoopSM : public TopologyBridge
{
public :
  
  LoopSM() ;
    //- The default constructor
  
  virtual ~LoopSM() ;
    //- The destructor
  
  virtual const type_info& topology_entity_type_info() const;
  
  virtual CubitStatus get_angle_metric(double& /*angle_metric*/)
    { return CUBIT_FAILURE; }
    //- Sets the value of angle_metric, ONLY if the LoopSM
    //- has a convenient way of calculating/storing it.
    //- Returns CUBIT_SUCCESS if it set the value of angle_metric,
    //- CUBIT_FAILURE otherwise.
  
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

