//------------------------------------------------------------------
//- Class:   AbstractTree
//- Author:  Kevin Albrecht 
//- Created: 2 October 2003
//-
//- Description:
//-   Abstract class to act as superclass of spacial data structures
//-   such as KDDTree and RTree.
//------------------------------------------------------------------

#ifndef ABSTRACTTREE_HPP
#define ABSTRACTTREE_HPP

#include "CubitVector.hpp"
#include "CubitBox.hpp"
#include "DLIList.hpp"

template <class Z> class AbstractTree 
{
  public:
    //// Other members
    typedef double (*DistSqFunc)(CubitVector &a, Z& b);

    //// Pure virtual methods; required by all subclasses 
    virtual CubitStatus add (Z data) = 0;
    virtual void set_tol (double tol) = 0;
    virtual double get_tol () = 0;
    virtual CubitBoolean remove (Z data) = 0;
    virtual CubitStatus find (const CubitBox &range_box, DLIList <Z> &range_members) = 0;
    virtual CubitStatus k_nearest_neighbor (CubitVector &q, int k, double &closest_dist,
                                            DLIList<Z> &nearest_neighbors,
                                            DistSqFunc dist_sq_point_data) = 0;

    //// Optional methods for subclasses   
    virtual CubitStatus balance () { return CUBIT_SUCCESS; };

    //// Destructor
    virtual ~AbstractTree () {};
};

#endif
