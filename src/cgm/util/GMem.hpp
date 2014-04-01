#ifndef GMEM_HPP
#define GMEM_HPP

// Include for CubitBoolean
#include "CubitDefines.h"
#include "CubitUtilConfigure.h"

//- A point in 3D space.
struct GPoint
{
   float x;
   float y;
   float z;
};

typedef struct
{
   GPoint at;
   GPoint from;
   GPoint up;
   float width;
   float height;
   CubitBoolean isCentered;
} CubitViewParam;

class CUBIT_UTIL_EXPORT GMem
{
private:
     // These are made private to make them essentially read-only,
     // although you still have direct access to the array elements.
   GPoint *pointList; // x, y, z vector
   int ptsSize;   // size of points vector
   int fListSize; // size of facetList vector
   int *facetList; // number of points in face, indices of points in face
   CubitBoolean pointsConsolidated;
   
   void consolidate_few_points( double tolerance );
   void consolidate_many_points( double tolerance );
   
public:
   int pointListCount;  // number valid points stored in points vector
   int fListCount;      // number valid integers stored in facetList.
    // If you always use triangles, this will be 4 times the number of polygons.
   
     // Constructor/Destructor take care of dynamic memory.
   GMem();
   ~GMem();
   
     // These are high-level functions that don't require you
     // to know exactly what the internals are
   void allocate_tri(int num_tri);
#ifdef BOYD15
   void allocate_more_tri(int num_tri);
#endif
   void allocate_polylines(int num_lines);
   void allocate_more_polylines(int num_lines);
   void clean_out();
   
     // The rest give you a little more direct control over the internals
   
     // Whether points at the same coords have been merged
   void points_consolidated(CubitBoolean yes_no)
      { pointsConsolidated = yes_no; }
   CubitBoolean points_consolidated()
      { return pointsConsolidated; }
   void consolidate_points( double tolerance );
   
     // Access to the arrays without letting the client
     // change the pointer directly.
   GPoint* point_list()
      { return pointList; }
   int* facet_list()
      { return facetList; }
   
     // Swap one array for another.  Requiring a function call
     // instead of direct access makes it more likely that
     // dynamic memory and count values will be handled correctly.
     //    Important!!!!!!!!!!!
     // The array passed in must have been allocated by 'new',
     // and must not be deleted outside of GMem!!!!!
   void replace_point_list(GPoint new_point_list[], int num_valid_points,
                           int array_size);
   void replace_facet_list(int new_facet_list[], int num_valid_entries,
                           int array_size);
   
     // Read only access to the array sizes
   int point_list_size()
      { return ptsSize; }
   int facet_list_size()
      { return fListSize; }

     // Copy constructor and operator=
   GMem(const GMem& from);
   GMem& operator=(const GMem& from);
};

#endif

