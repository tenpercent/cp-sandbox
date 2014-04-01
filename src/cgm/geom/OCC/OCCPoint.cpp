//-------------------------------------------------------------------------
// Filename      : OCCPoint.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Steven J. Owen
//
// Creation Date : 07/15/00
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "OCCPoint.hpp"
#include "OCCCurve.hpp"
#include "OCCQueryEngine.hpp"
#include "CastTo.hpp"
#include "CubitSimpleAttrib.hpp"
#include "BRep_Tool.hxx"
#include "TopExp.hxx"
#include "TopoDS_Edge.hxx"
#include "TopoDS.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "TopTools_IndexedDataMapOfShapeListOfShape.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "BRepAlgoAPI_BooleanOperation.hxx"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the location
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/08/07
//-------------------------------------------------------------------------
OCCPoint::OCCPoint( const CubitVector &location )
{
  gp_Pnt pt = gp_Pnt( location.x(), location.y(), location.z());
  TopoDS_Vertex theVertex = BRepBuilderAPI_MakeVertex(pt);
  myTopoDSVertex = new TopoDS_Vertex(theVertex);
}

 //-------------------------------------------------------------------------
// Purpose       : The constructor with a gp_Pnt point
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/19/07
//-------------------------------------------------------------------------
OCCPoint::OCCPoint( gp_Pnt& thePoint )
{
  TopoDS_Vertex theVertex = BRepBuilderAPI_MakeVertex(thePoint);
  myTopoDSVertex = new TopoDS_Vertex(theVertex);
}

//-------------------------------------------------------------------------
// Purpose       : The destructor. 
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/08/07
//-------------------------------------------------------------------------
OCCPoint::~OCCPoint() 
{
  if (myTopoDSVertex)
    delete (TopoDS_Vertex*)myTopoDSVertex;
}

void OCCPoint::set_TopoDS_Vertex(TopoDS_Vertex vertex)
{
  if(vertex.IsEqual(*myTopoDSVertex))
    return;
  TopoDS_Vertex* the_vertex = new TopoDS_Vertex(vertex);
  if(myTopoDSVertex)
    delete (TopoDS_Vertex*)myTopoDSVertex;
  myTopoDSVertex = the_vertex;
}
//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
void OCCPoint::append_simple_attribute_virt(CubitSimpleAttrib *csa)
  { OCCAttribSet::append_attribute(csa, *myTopoDSVertex); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
void OCCPoint::remove_simple_attribute_virt(CubitSimpleAttrib *csa)
  { OCCAttribSet::remove_attribute(csa, *myTopoDSVertex); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
void OCCPoint::remove_all_simple_attribute_virt()
  { OCCAttribSet::remove_attribute(NULL, *myTopoDSVertex); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
CubitStatus OCCPoint::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                               csa_list)
  { return OCCAttribSet::get_attributes(*myTopoDSVertex,csa_list);
  }

CubitStatus OCCPoint::get_simple_attribute(const CubitString& name,
                                     DLIList<CubitSimpleAttrib*>& csa_list )
  { return OCCAttribSet::get_attributes( name, *myTopoDSVertex, csa_list );
  }

//-------------------------------------------------------------------------
// Purpose       : Returns the coordinates of this Point. 
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/08/07
//-------------------------------------------------------------------------
CubitVector OCCPoint::coordinates() const
{
  gp_Pnt pt = BRep_Tool::Pnt(*myTopoDSVertex);
  CubitVector p(pt.X(), pt.Y(), pt.Z());
  return p;
}

CubitBoolean OCCPoint::is_equal(OCCPoint & other, double Tol)
{
  gp_Pnt pt = BRep_Tool::Pnt(*myTopoDSVertex);
  TopoDS_Vertex otherVertex = *(other.get_TopoDS_Vertex());
  const gp_Pnt otherPnt = BRep_Tool::Pnt(otherVertex);
  return pt.IsEqual(otherPnt, Tol);
}   
  
double OCCPoint::distance(OCCPoint & other)
{
  gp_Pnt pt = BRep_Tool::Pnt(*myTopoDSVertex);
  return pt.Distance(BRep_Tool::Pnt(*(other.get_TopoDS_Vertex())));
}

double OCCPoint::SquareDistance (OCCPoint & other)
{
  gp_Pnt pt = BRep_Tool::Pnt(*myTopoDSVertex);
  return pt.SquareDistance(BRep_Tool::Pnt(*(other.get_TopoDS_Vertex())));
}
//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: AcisGeometryEngine
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
GeometryQueryEngine* OCCPoint::get_geometry_query_engine() const
{
  return OCCQueryEngine::instance();
}                 

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
CubitBox OCCPoint::bounding_box() const 
{
  CubitVector temp_vector = this->coordinates();
  CubitBox temp_box(temp_vector);
  return temp_box;
}


void OCCPoint::get_parents_virt( DLIList<TopologyBridge*>& parents ) 
{
  OCCQueryEngine* oqe = (OCCQueryEngine*) get_geometry_query_engine();
  OCCCurve * curve = NULL;
  DLIList <OCCCurve* > *curves = oqe->CurveList;
  TopTools_IndexedDataMapOfShapeListOfShape M;
  for(int i = 0; i <  curves->size(); i++)
  {
     curve = curves->get_and_step();
     TopExp::MapShapesAndAncestors(*(curve->get_TopoDS_Edge()),
                                   TopAbs_VERTEX, TopAbs_EDGE, M);
     if (!M.Contains(*(get_TopoDS_Vertex())))
       continue;

     const TopTools_ListOfShape& ListOfShapes =
                                M.FindFromKey(*(get_TopoDS_Vertex()));
     if (!ListOfShapes.IsEmpty())
     {
         TopTools_ListIteratorOfListOfShape it(ListOfShapes) ;
         for (;it.More(); it.Next())
         {
           TopoDS_Edge Edge = TopoDS::Edge(it.Value());
           int k = oqe->OCCMap->Find(Edge);
           parents.append_unique((OCCPoint*)(oqe->OccToCGM->find(k))->second);
         }
     }
  }
}
void OCCPoint::get_children_virt( DLIList<TopologyBridge*>& ) 
  {  }

//----------------------------------------------------------------
// Function: to update the core vertex
//           for any movement of the body/surface/curve/vertex.
// Author: Jane Hu
//----------------------------------------------------------------
void OCCPoint::update_OCC_entity( BRepBuilderAPI_Transform *aBRepTrsf,
                                  BRepAlgoAPI_BooleanOperation *op)
{
  if (this->myMarked == CUBIT_TRUE)
    return;

  assert(aBRepTrsf != NULL || op != NULL);

  TopoDS_Shape shape;
  if(aBRepTrsf)
    shape = aBRepTrsf->ModifiedShape(*get_TopoDS_Vertex());
  else
  {
    TopTools_ListOfShape shapes;
    shapes.Assign(op->Modified(*get_TopoDS_Vertex()));
    if(shapes.Extent() == 0)
      shapes.Assign(op->Generated(*get_TopoDS_Vertex()));
    if(shapes.Extent() == 1)
      shape = shapes.First();
    else if(shapes.Extent() > 1)
    {
      //update all attributes first.
      TopTools_ListIteratorOfListOfShape it;
      it.Initialize(shapes);
      for(; it.More(); it.Next())
      {
        shape = it.Value();
        OCCQueryEngine::instance()->copy_attributes(*get_TopoDS_Vertex(),shape);
      }
      shape = shapes.First();
    }
    else if(op->IsDeleted(*get_TopoDS_Vertex()))
      ;
    else
      return ;
  }
  TopoDS_Vertex vertex;
  if(!shape.IsNull())
    vertex = TopoDS::Vertex(shape);

  OCCQueryEngine::instance()->update_OCC_map(*myTopoDSVertex, vertex);

  set_myMarked(CUBIT_TRUE);
  set_TopoDS_Vertex(vertex);
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
