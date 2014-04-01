
#include "CubitCoordinateSystem.hpp"

#include <assert.h>

#include "CubitMessage.hpp"

//! default constructor, initializes to rectangular type with an identity transform
CubitCoordinateSystem::CubitCoordinateSystem(int id)
 : mId(id), mType(Rectangular), mMatrix(CubitTransformMatrix()), mReference(0)
{
}

CubitCoordinateSystem::~CubitCoordinateSystem()
{
  if(mReference)
    mReference->mUses.remove(this);
  
  // if we have uses, make them reference this reference
  while(mUses.size())
  {
    CubitCoordinateSystem* sys = mUses.get();
    sys->set_transform(sys->mMatrix * this->mMatrix);
    sys->mReference = this->mReference;
    if(mReference)
      mReference->mUses.append(sys);
    mUses.remove();
  }
}

CubitCoordinateSystem::Type CubitCoordinateSystem::get_type() const
{
  return mType;
}

void CubitCoordinateSystem::set_type(CubitCoordinateSystem::Type type)
{
  mType = type;
}


//! return the transform of this coordinate system
const CubitTransformMatrix& CubitCoordinateSystem::get_transform() const
{
  return mMatrix;
}

void CubitCoordinateSystem::set_transform(const CubitTransformMatrix& mat)
{
  mMatrix = mat;
}


//! return the transform of this coordinate system 
//! with respect to the global coordinate system
CubitTransformMatrix CubitCoordinateSystem::get_concatenated_transform() const
{
  CubitTransformMatrix mat;
  mat = mMatrix;
  CubitCoordinateSystem* ref;
  for(ref = mReference; ref != 0; ref = ref->mReference)
    mat = mat * ref->mMatrix;
  return mat;
}


