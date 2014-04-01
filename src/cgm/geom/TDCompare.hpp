//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename      : TDCompare.hpp
//
// Purpose       : This class is meant to be the container of temporary
//                 data during compare.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 11/25/96
//
// Owner         : Raikanta Sahu
//-------------------------------------------------------------------------

#ifndef TD_COMPARE_HPP
#define TD_COMPARE_HPP

// ********** BEGIN STANDARD INCLUDES         **********
// ********** END STANDARD INCLUDES           **********

// ********** BEGIN MOTIF INCLUDES            **********
// ********** END MOTIF INCLUDES              **********

// ********** BEGIN OPEN INVENTOR INCLUDES    **********
// ********** END OPEN INVENTOR INCLUDES      **********

// ********** BEGIN ACIS INCLUDES             **********
// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********
#include "ToolData.hpp"
// ********** END CUBIT INCLUDES              **********

// ********** BEGIN FORWARD DECLARATIONS      **********
class RefEntity ;
// ********** END FORWARD DECLARATIONS        **********

// ********** BEGIN MACRO DEFINITIONS         **********
// ********** END MACRO DEFINITIONS           **********

// ********** BEGIN ENUM DEFINITIONS          **********
// ********** END ENUM DEFINITIONS            **********

class TDCompare : public ToolData
{

// ********** BEGIN FRIEND CLASS DECLARATIONS **********
// ********** END FRIEND CLASS DECLARATIONS   **********

   public:

      TDCompare() ;
      //- Default constructor. 

      virtual ~TDCompare() ;
      //- Destructor

      static int is_compare(const ToolData* td) 
	  {return (dynamic_cast<TDCompare*>(const_cast<ToolData*>(td)) != NULL);}
  
      void set_compare_partner(RefEntity* partner) ;
      RefEntity* get_compare_partner() const ; 

      
   protected: 

   private:

      RefEntity* comparePartner_ ;
} ;

// ********** BEGIN INLINE FUNCTIONS          **********

//-------------------------------------------------------------------------
// Purpose       : Default constructor
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 11/25/96
//-------------------------------------------------------------------------
inline
TDCompare::TDCompare() : comparePartner_(NULL) 
{
}

//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 11/25/96
//-------------------------------------------------------------------------
inline
TDCompare::~TDCompare()
{
}

//-------------------------------------------------------------------------
// Purpose       : Set the RefEntity that the one that this TDCompare belongs
//                 to compares with.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 11/25/96
//-------------------------------------------------------------------------
inline
void TDCompare::set_compare_partner(RefEntity* partner)
{
   comparePartner_ = partner ;
}

//-------------------------------------------------------------------------
// Purpose       : Get the "compared" partner to the RefEntity that owns
//                 this TDCompare object.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 11/25/96
//-------------------------------------------------------------------------
inline
RefEntity* TDCompare::get_compare_partner() const
{
   return comparePartner_ ;
}

// ********** END INLINE FUNCTIONS            **********

// ********** BEGIN FRIEND FUNCTIONS          **********
// ********** END FRIEND FUNCTIONS            **********

// ********** BEGIN EXTERN FUNCTIONS          **********
// ********** END EXTERN FUNCTIONS            **********

// ********** BEGIN HELPER CLASS DECLARATIONS **********
// ********** END HELPER CLASS DECLARATIONS   **********

#endif //TD_COMPARE_HPP

