//- Class: RefEntityNameMap
//-
//- Description: This class maintains an map between a RefEntity
//- name and the corresponding RefEntity. Note that there may be multiple
//- names for a RefEntity, but a single name may only refer to one RefEntity.
//-
//- Owner: Greg Sjaardema
//- Checked by: 
//- Version: $Id: 

#ifndef REFENTITYNAMEMAP_HPP
#define REFENTITYNAMEMAP_HPP

#include "CubitString.hpp"
#include "CubitGeomConfigure.h"

class RefEntity;

class CUBIT_GEOM_EXPORT RefEntityNameMap {
public:
  RefEntityNameMap(const CubitString &initial_key,
			   const RefEntity *intial_value);
  //- Constructor

  RefEntity    *value() const;
  void          value(const RefEntity *new_value);
  CubitString   key()   const;
  CubitString  *ptr_to_key();

private:
  CubitString keyField;
  RefEntity  *valueField;
};

inline CubitString  RefEntityNameMap::key()   const {return keyField;}
inline CubitString *RefEntityNameMap::ptr_to_key()  
{return &keyField;}
inline RefEntity   *RefEntityNameMap::value() const {return valueField;}
inline void         RefEntityNameMap::value(const RefEntity *new_value)
{
  valueField = (RefEntity *)new_value;
}

inline RefEntityNameMap::RefEntityNameMap(const CubitString &initial_key,
					  const RefEntity *intial_value)
  {
    keyField   = initial_key;
    valueField = (RefEntity *)intial_value;
  }
#endif

