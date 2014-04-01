//-------------------------------------------------------------------------
// Filename      : OCCPoint.hpp
//
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Steven J. Owen
//
// Creation Date : 08/02/96
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

#ifndef POINT_OCCPOINT_HPP
#define POINT_OCCPOINT_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********
// ********** BEGIN CUBIT INCLUDES         **********
#include "Point.hpp"
#include <stdio.h>
#include "gp_Pnt.hxx"
#include "TopoDS_Vertex.hxx"
#include "BRepBuilderAPI_MakeVertex.hxx"
#include "OCCAttribSet.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class CubitSimpleAttrib;
class BRepBuilderAPI_Transform;
class BRepAlgoAPI_BooleanOperation;
// ********** END FORWARD DECLARATIONS     **********

class OCCPoint : public Point
{
private:

  TopoDS_Vertex *myTopoDSVertex;
  CubitBoolean myMarked ; 
  OCCAttribSet attribSet;

public :
  
  OCCPoint(const CubitVector &location );
    //I- CubitVector &location
    //I- location of point (creates a CubiPoint).
    //I- DLIList<Curve*> curves
    //I- curves attaced to point

  OCCPoint(TopoDS_Vertex* thePoint ):myTopoDSVertex(thePoint), myMarked(CUBIT_FALSE){};
     //I- gp_Pnt *thePoint
     //I- pointer to the TopoDS_Vertex associated with OCCPoint
     //I- DLIList<Curve*> curves
     //I- curves attaced to point

  OCCPoint(gp_Pnt& thePoint );
    //I- gp_Pnt *thePoint
    //I- pointer to the TopoDS_Vertex associated with OCCPoint
    //I- DLIList<Curve*> curves
    //I- curves attaced to point

  virtual ~OCCPoint();
    //- The destructor

  void set_myMarked(CubitBoolean marked) {myMarked = marked;}

  TopoDS_Vertex *get_TopoDS_Vertex(){return myTopoDSVertex; }
  void set_TopoDS_Vertex(TopoDS_Vertex vertex);

  virtual void append_simple_attribute_virt(CubitSimpleAttrib*);
    //R void
    //I
    //I-
    //I- that is to be appended to this OSME object.
    //- The purpose of this function is to append a
    //- attribute to the OSME. The  is attached to each of the
    //- underlying solid model entities this one points to.

  virtual void remove_simple_attribute_virt(CubitSimpleAttrib*);
    //R void
    //I CubitSimpleAttrib*
    //I- A reference to a CubitSimpleAttrib object which is the object
    //I- that is to be removed to this OSME object.
    //- The purpose of this function is to remove a simple
    //- attribute from the OSME. The attribute is attached to each of the
    //- underlying solid model entities this one points to.

  virtual void remove_all_simple_attribute_virt();
    //R void
    //I-
    //- The purpose of this function is to remove all simple
    //- attributes from the OSME.

  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>&);
  virtual CubitStatus get_simple_attribute(const CubitString& name,
                                           DLIList<CubitSimpleAttrib*>&);
    //R CubitSimpleAttrib*
    //R- the returned cubit simple attribute.
    //- The purpose of this function is to get the attributes
    //- of the geometry entity. The name is attached to the underlying solid
    //- model entity(ies) this one points to.
    //- MJP Note:
    //- This is the code that implements the requirement that names
    //- of VGI Entities propagate across solid model boolean
    //- operations.  The success of this relies, of course, on the underlying
    //- solid modeler being able to propagate attributes across
    //- such operations on its entities. If it cannot, then "names"
    //- of VGI entities will not propagate.

  virtual CubitVector coordinates() const;
    //R CubitVector
    //R- Contains the coordinate values {x y z} of this Point
    //- Returns the spatial coordinates of this Point.
  
  CubitBoolean is_equal( OCCPoint & other, double Tol);

  double distance( OCCPoint & other);
  double SquareDistance ( OCCPoint & other);

  virtual CubitBox bounding_box() const ;
    // see comments in GeometryEntity.hpp
  
  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  virtual void get_parents_virt( DLIList<TopologyBridge*>& parents );
  virtual void get_children_virt( DLIList<TopologyBridge*>& children );

  void update_OCC_entity( BRepBuilderAPI_Transform *aBRepTrsf,
                          BRepAlgoAPI_BooleanOperation *op = NULL);
};


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif
