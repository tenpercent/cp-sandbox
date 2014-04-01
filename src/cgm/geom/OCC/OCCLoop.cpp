//-------------------------------------------------------------------------
// Filename      : OCCLoop.cpp
//
// Purpose       : Loops for Facet-based geometry module
//
// Special Notes :
//
// Creator       : Steven J. Owen
//
// Creation Date : 12/06/00
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

#include "CastTo.hpp"
#include "OCCLoop.hpp"
#include "OCCQueryEngine.hpp"
#include "CoEdgeSM.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCShell.hpp"
#include "OCCSurface.hpp"
#include "OCCCoEdge.hpp"
#include "OCCCurve.hpp"
#include "OCCPoint.hpp"

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include "TopExp_Explorer.hxx"
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
// Purpose       : The constructor with a pointer to the TopoDS_Wire.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCLoop::OCCLoop( TopoDS_Wire *theWire )
{
  myTopoDSWire = theWire;
}

//-------------------------------------------------------------------------
// Purpose       : The default destructor.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCLoop::~OCCLoop()
{
  disconnect_all_curves();
  if (myTopoDSWire)
    delete (TopoDS_Wire*)myTopoDSWire;
}

void OCCLoop::set_TopoDS_Wire(TopoDS_Wire loop)
{
   if(loop.IsEqual(*myTopoDSWire))
     return;
   TopoDS_Wire* the_wire = new TopoDS_Wire(loop);
   if(myTopoDSWire)
     delete (TopoDS_Wire*)myTopoDSWire;
   myTopoDSWire = the_wire;
}
//-------------------------------------------------------------------------
// Purpose       : Tear down topology
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 11/29/07
//-------------------------------------------------------------------------
void OCCLoop::disconnect_all_curves()
{
  myCoEdgeList.clean_out();
}

OCCCoEdge* OCCLoop::remove_coedge(OCCCoEdge *coedge)
{ 
  return myCoEdgeList.remove(coedge);
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the
//                 underlying solid model entity this one points to.
//
//
// Special Notes :
//
//-------------------------------------------------------------------------
void OCCLoop::append_simple_attribute_virt(CubitSimpleAttrib* /*csattrib_ptr*/)
{
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple
//                 attribute attached to this geometry entity. The name is
//                 removed from the underlying BODY this points to.
//
// Special Notes :
//
//-------------------------------------------------------------------------
void OCCLoop::remove_simple_attribute_virt(CubitSimpleAttrib* /*csattrib_ptr*/)
{
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes :
//
//-------------------------------------------------------------------------
void OCCLoop::remove_all_simple_attribute_virt()
{
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the
//                 attributes attached to this geometry entity. The name is
//                 attached to the underlying BODY this points to.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitStatus OCCLoop::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                                 /*cubit_simple_attrib_list*/)
{
  return CUBIT_FAILURE;
}
CubitStatus OCCLoop::get_simple_attribute(const CubitString&,
                                          DLIList<CubitSimpleAttrib*>&)
  { return CUBIT_FAILURE; }

void OCCLoop::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  OCCQueryEngine* oqe = (OCCQueryEngine*) get_geometry_query_engine();
  OCCSurface * surf = NULL;
  DLIList <OCCSurface* > *surfs = oqe->SurfaceList;
  TopTools_IndexedDataMapOfShapeListOfShape M;
  for(int i = 0; i <  surfs->size(); i++)
  {
     surf = surfs->get_and_step();
     TopExp::MapShapesAndAncestors(*(surf->get_TopoDS_Face()),
                                   TopAbs_WIRE, TopAbs_FACE, M);
     if (!M.Contains(*(get_TopoDS_Wire())))
	continue;

     const TopTools_ListOfShape& ListOfShapes =
                                M.FindFromKey(*(get_TopoDS_Wire()));
     if (!ListOfShapes.IsEmpty())
     {
         TopTools_ListIteratorOfListOfShape it(ListOfShapes) ;
         for (;it.More(); it.Next())
         {
           TopoDS_Face Face = TopoDS::Face(it.Value());
           int k = oqe->OCCMap->Find(Face);
           parents.append_unique((OCCSurface*)(oqe->OccToCGM->find(k))->second);
         }
     }
  }
}

void OCCLoop::get_children_virt( DLIList<TopologyBridge*>& children )
{
  CAST_LIST_TO_PARENT(myCoEdgeList, children);
}

//-------------------------------------------------------------------------
// Purpose       : compute bounding box of loop
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox OCCLoop::bounding_box() const
{
   CubitBox box;
   for (int i = myCoEdgeList.size(); i > 0; i--)
   {
      DLIList<OCCCoEdge*> coedges = myCoEdgeList;
      OCCCoEdge* coedge = coedges.get_and_step();
      box |= coedge->curve()->bounding_box();
   }
   return box;
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: OCCQueryEngine
//
// Special Notes :
//
//-------------------------------------------------------------------------
GeometryQueryEngine* OCCLoop::get_geometry_query_engine() const
{
   return OCCQueryEngine::instance();
}                


//----------------------------------------------------------------
// Function: to update the core Loop
//           for any movement of the body/surface.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCLoop::update_OCC_entity( BRepBuilderAPI_Transform *aBRepTrsf,
                                        BRepAlgoAPI_BooleanOperation *op)
{
  assert(aBRepTrsf != NULL || op != NULL);

  TopoDS_Shape shape;
  CubitBoolean need_update = CUBIT_TRUE;
  if (aBRepTrsf)
    shape = aBRepTrsf->ModifiedShape(*get_TopoDS_Wire());
  else
  {
    TopTools_ListOfShape shapes; 
    shapes.Assign(op->Modified(*get_TopoDS_Wire()));
    if(shapes.Extent() == 0)
         shapes.Assign(op->Generated(*get_TopoDS_Wire()));
    if(shapes.Extent())
      shape = shapes.First();
    else if (op->IsDeleted(*get_TopoDS_Wire()))
      ;
    else
      need_update = CUBIT_FALSE;
  }

  //set the curves
  for (int i = 1; i <= myCoEdgeList.size(); i++)
  {
     OCCCurve *curve = CAST_TO(myCoEdgeList.get_and_step()->curve(), OCCCurve);
     curve->update_OCC_entity(aBRepTrsf, op);
  }
  TopoDS_Wire loop;
  if (need_update)
  {
    loop = TopoDS::Wire(shape);
    OCCQueryEngine::instance()->update_OCC_map(*myTopoDSWire, loop);
  }
  return CUBIT_SUCCESS;
}

