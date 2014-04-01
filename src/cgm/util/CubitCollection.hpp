//- Class:        CubitCollection
//- Description:  An generic base class to create collections
//-               of pointers to any type of object.
//- Owner:        Bill Bohnhoff
//- Checked by:
//- Version: $Id: 

#ifndef CUBIT_COLLECTION_HPP
#define CUBIT_COLLECTION_HPP

#include "CubitContainer.hpp"
#include "CubitUtilConfigure.h"


class CUBIT_UTIL_EXPORT CubitCollection: public CubitContainer
{
  public:

    //- Heading: Constructors and Destructor

    CubitCollection();
    CubitCollection(int num_items);
    virtual ~CubitCollection();
   
  protected:

};

#endif // CUBIT_COLLECTION_HPP

