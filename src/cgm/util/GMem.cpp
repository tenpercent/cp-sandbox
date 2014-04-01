#include "GMem.hpp"
#include <string.h> // To define NULL
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "OctTree.hpp"

GMem::GMem()
{
   ptsSize = 0;
   fListSize = 0;
   pointListCount = 0;
   fListCount = 0;
   pointList = NULL;
   facetList = NULL;
   pointsConsolidated = CUBIT_FALSE;
}

GMem::~GMem()
{
   if (pointList)
      delete [] pointList;
   if (facetList)
      delete [] facetList;
}

void GMem::allocate_tri(int num_tri)
{
   if (num_tri*3 > ptsSize) 
   {
      if (pointList)
         delete [] pointList;
      ptsSize = 3 * num_tri;
      pointList = new GPoint[ptsSize];
   } 
   if (num_tri*4 > fListSize) 
   {
      if (facetList)
         delete [] facetList;
      fListSize = 4 * num_tri;
      facetList = new int[fListSize]; 
   }
}

// Allocates space for a polyline consisting of
// {num_lines} segments.  This is num_lines + 1 pointList.
void GMem::allocate_polylines(int num_lines)
{
   num_lines++;
   if (num_lines > ptsSize)
   {
      if (pointList)
         delete [] pointList;
      ptsSize = num_lines;
      pointList = new GPoint[ptsSize];
   }
}

// Allocates additional space for a bigger polyline.
// The total number of points will now be
// pointListCount + num_additional_segs.
void GMem::allocate_more_polylines(int num_additional_segs)
{
     // Make sure we need to do something
   if (num_additional_segs <= 0)
      return;
     // If there are no points, we need space for one additional point
   if (pointListCount == 0)
      num_additional_segs++;
     // Make sure we need to grow
   if (ptsSize >= pointListCount + num_additional_segs)
      return;

     // Make a big enough array of points
   GPoint *new_points = new GPoint[pointListCount + num_additional_segs];
     // Copy the old data to the new array
   if (pointList)
   {
      memcpy (new_points, pointList, pointListCount*sizeof(GPoint));
      delete []pointList;
   }
     // Store the new array in 'pointList'
   pointList = new_points;
   ptsSize = pointListCount + num_additional_segs;
}

void GMem::clean_out()
{
   pointListCount = 0;
   fListCount = 0;
}


void GMem::replace_point_list(GPoint new_point_list[], int num_valid_points,
                              int array_size)
{
     // First, delete the old memory
   if (pointList)
      delete [] pointList;
     // Now replace it with the new
   pointList = new_point_list;
     // Set the array size
   ptsSize = array_size;
     // Set the number of valid entries
   pointListCount = num_valid_points;
}

void GMem::replace_facet_list(int new_facet_list[], int num_valid_entries,
                              int array_size)
{
     // First, delete the old memory
   if (facetList)
      delete [] facetList;
     // Now replace it with the new
   facetList = new_facet_list;
     // Set the array size
   fListSize = array_size;
     // Set the number of valid entries
   fListCount = num_valid_entries;
}

// copy constructor
GMem::GMem(const GMem& from)
{
     // Set the array pointers to NULL
   pointList = NULL;
   facetList = NULL;
     // Set one equal to the other
   *this = from;
}

// Equality operator
GMem& GMem::operator=(const GMem& from)
{
   if (this != &from)
   {
        // Make a copy of the point array
      GPoint* temp1 = new GPoint[from.ptsSize];
      memcpy (temp1, from.pointList, from.ptsSize*sizeof(GPoint));
        // Put it in the receiving GMem
      replace_point_list(temp1, from.pointListCount, from.ptsSize);
        // Make a copy of the facet array
      int* temp2 = new int[from.fListSize];
      memcpy (temp2, from.facetList, from.fListSize*sizeof(int));
        // Put it in the receiving GMem
      replace_facet_list(temp2, from.fListCount, from.fListSize);
        // Set whether it's consolidated
      pointsConsolidated = from.pointsConsolidated;
   }
   return *this;
}


void GMem::consolidate_points( double tolerance )
{
  if ( pointListCount < 1000 )
    consolidate_few_points( tolerance );
  else 
    consolidate_many_points( tolerance );
}

void GMem::consolidate_few_points( double tolerance )
{
  const double tolsqr = tolerance * tolerance;
  int* index_map = new int[pointListCount];
    
    // Consolidate the point list.  index_map is used
    // to maintain a map between the old index of a point
    // (the index into index_map) and the new index of a
    // point (the value in index_map).
  int write = 0, read, comp;
  for( read = 0; read < pointListCount; read++ )
  {
    const GPoint& pti = pointList[read];
    for( comp = 0; comp < write; comp++ )
    {
      const GPoint& ptj = pointList[comp];
      double x = pti.x - ptj.x;
      double y = pti.y - ptj.y;
      double z = pti.z - ptj.z;
      if( (x*x+y*y+z*z) <= tolsqr )
        break;
    }
    
    index_map[read] = comp;
    if( comp == write )
    {
      pointList[comp] = pointList[read];
      write++;
    }
  }
  pointListCount = write;  

    // Update the facet list using values from index_map.
  int *itor = facetList;
  const int* end = facetList + fListCount;
  while( itor < end )
    for( int count = *(itor++); count--; itor++ )
      *itor = index_map[*itor];

  delete [] index_map;
  pointsConsolidated = CUBIT_TRUE;
}

class GPointOctTreeEval {
  public: static inline CubitVector coordinates(GPoint* p)
    { return CubitVector( p->x, p->y, p->z ); }
};

void GMem::consolidate_many_points( double tolerance )
{
  const double tolsqr = tolerance * tolerance;
  
    // build OctTree
  DLIList<GPoint*> point_list(pointListCount);
  GPoint* p_itor = pointList;
  GPoint* p_end  = p_itor + pointListCount;
  while( p_itor < p_end )
    point_list.append( p_itor++ );
  OctTree<GPoint, GPointOctTreeEval> tree( point_list, tolerance );
  point_list.clean_out();
  
    // Consolidate the point list.  index_map is used
    // to maintain a map between the old index of a point
    // (the index into index_map) and the new index of a
    // point (the value in index_map).
  int* index_map = new int[pointListCount];
  GPoint* new_array = new GPoint[pointListCount];
  int read, write = 0;
  for ( read = 0; read < pointListCount; read++)
  {
    GPoint* pt = pointList + read;
    CubitVector v(pt->x, pt->y, pt->z);
    index_map[read] = write;
    point_list.clean_out();
    tree.nodes_near( v, point_list );
    while ( point_list.size() )
    {
      GPoint* p = point_list.pop();
      int index = p - pointList;
      assert( (index >= 0) && (index < pointListCount) );
      CubitVector v2(p->x, p->y, p->z);
      if ( (index < read) && ((v - v2).length_squared() < tolsqr) )
      {
        index_map[read] = index_map[index];
        break;
      }
    }
    
    if ( index_map[read] == write )
    {
      new_array[write++] = pointList[read];
    }
  }
  pointListCount = write;  
  delete [] pointList;
  pointList = new_array;
  
    // Update the facet list using values from index_map.
  int *itor = facetList;
  const int* end = facetList + fListCount;
  while( itor < end )
    for( int count = *(itor++); count--; itor++ )
      *itor = index_map[*itor];

  delete [] index_map;
  pointsConsolidated = CUBIT_TRUE;
}
