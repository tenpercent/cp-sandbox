//-------------------------------------------------------------------------
// Filename      : CAMergeStatus.hpp
//
// Purpose       : Save autoMergeStatus from RefEntities.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/06/02
//-------------------------------------------------------------------------

#ifndef CA_MERGE_STATUS_HPP
#define CA_MERGE_STATUS_HPP

#include "CubitAttrib.hpp"

// TODO - decide where to define enum AutoMergeStatus
#include "RefEntity.hpp"

class CUBIT_GEOM_EXPORT CAMergeStatus : public CubitAttrib
{
  private:
  
  AutoMergeStatus status;


  public:
  
  CAMergeStatus( RefEntity* attrib_owner );
  
  CAMergeStatus( RefEntity*, CubitSimpleAttrib* );

  virtual ~CAMergeStatus();
  
  virtual CubitSimpleAttrib* cubit_simple_attrib();
  
  virtual CubitStatus actuate();
  
  virtual CubitStatus update();
  
  virtual CubitStatus reset();
  
  virtual int int_attrib_type();
  
  virtual const type_info& entity_type_info() const;
  
};

CubitAttrib* CAMergeStatus_creator(RefEntity* entity, CubitSimpleAttrib *p_csa);

#endif

