//- Class:       IndexedDouble
//- Description: IndexedDouble associates an index with a double value for
//-              use with sorted or other lists.
//- Owner:       Scott Mitchell
//- Version: $Id: 

#ifndef INDEXED_DOUBLE
#define INDEXED_DOUBLE

class IndexedDouble {
public: 
  int myIndex;
  int index() {return myIndex;}
  
  double myDouble;
  double val() {return myDouble;}

  IndexedDouble( int index_set, double double_set )
  { myIndex = index_set; myDouble = double_set; }
  //- constructor
  
};

#endif

