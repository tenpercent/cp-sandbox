//- Class:        CubitContainer
//- Description:  An abstract class template to create containers
//-               of pointers to any object.
//- Owner:        Bill Bohnhoff
//- Checked by:
//- Version: $Id: 

#include "CubitDefines.h"
#include "CubitContainer.hpp"
#include "CubitMessage.hpp"


CubitContainer::CubitContainer()
{
  numberItems = 0;
}

CubitContainer::CubitContainer(int i)
{
  if(i<1)
  {
    PRINT_ERROR("Improper number of items "
                                "specified for CubitContainer\n");
    exit(1);
  }

  numberItems = i;
}

CubitContainer::~CubitContainer() { }


