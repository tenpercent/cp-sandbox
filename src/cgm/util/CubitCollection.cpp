//- Class:        CubitCollection
//- Description:  An generic base class to create collections
//-               of pointers to any object.
//- Owner:        Bill Bohnhoff
//- Checked by:
//- Version: $Id: 

#include "CubitCollection.hpp"


CubitCollection::CubitCollection()
{
  numberItems = 0;
}

CubitCollection::CubitCollection(int i): CubitContainer(i) { }

CubitCollection::~CubitCollection() { }


