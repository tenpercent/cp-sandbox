//- Class:        CubitContainer
//- Description:  An abstract base class for generic containers
//-               of pointers to any type of object.
//- Owner:        Bill Bohnhoff
//- Checked by:
//- Version: $Id: 

#ifndef CUBIT_CONTAINER_HPP
#define CUBIT_CONTAINER_HPP

#include "CubitUtilConfigure.h"


class CUBIT_UTIL_EXPORT CubitContainer
{
  public:

    //- Heading: Constructors and Destructor

    CubitContainer();
    CubitContainer(int num_items);
    virtual ~CubitContainer();
   
#ifdef BOYD15
    int number_items() const;
    //- returns the number of items in the container
#endif

  protected:

    int numberItems;
};

#endif // CUBIT_CONTAINER_HPP

