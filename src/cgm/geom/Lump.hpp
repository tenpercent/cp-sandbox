//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : Lump.hpp
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

#ifndef LUMP_HPP
#define LUMP_HPP


#include "CubitDefines.h"
#include "GeometryEntity.hpp"
#include "GeometryDefines.h"
#include "CubitGeomConfigure.h"


class CUBIT_GEOM_EXPORT Lump : public GeometryEntity
{
   public :

      Lump() ;
      //- The default constructor

      virtual ~Lump() ;
      //- The destructor
  
      virtual const type_info& topology_entity_type_info() const;

  virtual GeometryType geometry_type()
  {return UNDEFINED_LUMP_TYPE;};
      //R GeometryType (enum)
      //R- The enumerated type of the geometric representation

   virtual CubitStatus mass_properties( CubitVector &centroid, double &volume ) = 0;

   protected: 

   private:
} ;


#endif

