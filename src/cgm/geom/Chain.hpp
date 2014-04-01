//-------------------------------------------------------------------------
// Filename      : Chain.hpp
//
// Purpose       : Represents a chain of vertices of an edge of a model.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef CHAIN_HPP
#define CHAIN_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitUtil.hpp"

#include "GroupingEntity.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class RefVertex;
// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT Chain : public GroupingEntity
{
public :
  
  Chain() ;
    //- The default constructor
  
  virtual ~Chain() ;
    //- The destructor

  DagType dag_type() const { return DagType::chain_type(); }
  
  RefVertex* start_vertex();
  RefVertex* end_vertex();
    //R RefVertex*
    //R- Returned RefVertex pointer
    //- These functions return the RefVertex associated with 
    //- the first/last CoVertex'es of this Chain.
    
  CoVertex* start_co_vertex();
  CoVertex* end_co_vertex();
  
protected: 
  
private:
    Chain( const Chain& );
    void operator=( const Chain&);
};

// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

