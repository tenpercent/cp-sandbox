//-------------------------------------------------------------------------
// Filename      : OCCBody.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 7/18/00
//
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CastTo.hpp"
#include "BodySM.hpp"
#include "Body.hpp"
#include "OCCBody.hpp"
#include "CubitSimpleAttrib.hpp"
#include "OCCQueryEngine.hpp"
#include "DLIList.hpp"
#include "Surface.hpp"
#include "OCCSurface.hpp"
#include "CubitTransformMatrix.hpp"
#include "OCCPoint.hpp"
#include "OCCCurve.hpp"
#include "OCCCoEdge.hpp"
#include "OCCLoop.hpp"
#include "OCCShell.hpp"
#include "OCCLump.hpp"
#include "OCCModifyEngine.hpp"

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include "BRepBuilderAPI_Transform.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "gp_Ax1.hxx"
#include "gp_Ax2.hxx"
#include "Bnd_Box.hxx"
#include "BRepBndLib.hxx"
#include "TopExp_Explorer.hxx"
#include "BRep_Builder.hxx"
#include "TopoDS.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "Standard_Boolean.hxx"
#include "LocOpe_SplitShape.hxx"
//-------------------------------------------------------------------------
// Purpose       : A constructor with a list of lumps that are attached.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCBody::OCCBody(TopoDS_CompSolid *theShape, CubitBoolean isSheetBody,
                 OCCSurface* surface, OCCShell* shell)
{
  myTopoDSShape = theShape;
  IsSheetBody = isSheetBody; 
  mySheetSurface = surface;
  myShell = shell;
  update_bounding_box();
}

void OCCBody::set_sheet_surface(OCCSurface* surface)
{
  IsSheetBody = CUBIT_TRUE;
  mySheetSurface = surface;
  if(surface == NULL)
    IsSheetBody = CUBIT_FALSE;
}

void OCCBody::lumps(DLIList<Lump*>& my_lumps)
{
  myLumps.clean_out();
  myLumps += my_lumps;
}

void OCCBody::set_TopoDS_Shape( TopoDS_CompSolid theshape)
{
  if(theshape.IsEqual(*myTopoDSShape))
    return;

  TopoDS_CompSolid * the_comp = new TopoDS_CompSolid(theshape);
  if(myTopoDSShape)
    delete (TopoDS_CompSolid*)myTopoDSShape;
  myTopoDSShape = the_comp;
}

OCCBody::OCCBody(DLIList<Lump*>& my_lumps)
{
  myLumps += my_lumps;
  TopoDS_CompSolid* new_top = make_CompSolid(my_lumps);
  myTopoDSShape = new_top;
  IsSheetBody = CUBIT_FALSE;
  myShell = NULL;
  update_bounding_box();
}

TopoDS_CompSolid* OCCBody::make_CompSolid(DLIList<Lump*>& my_lumps)
{
  BRep_Builder B;
  TopoDS_CompSolid Co;
  B.MakeCompSolid(Co);
  for(int i = 0; i < myLumps.size(); i ++)
  {
     TopoDS_Solid * solid = CAST_TO(myLumps.get_and_step(), OCCLump)->get_TopoDS_Solid();
     B.Add(Co, *solid);
  }
  TopoDS_CompSolid* new_top = new TopoDS_CompSolid(Co);
  return new_top; 
}

OCCBody::~OCCBody() 
{
  if (myTopoDSShape)
    delete (TopoDS_CompSolid*)myTopoDSShape;
}

GeometryQueryEngine* OCCBody::get_geometry_query_engine() const
{
  return OCCQueryEngine::instance();
}

void OCCBody::append_simple_attribute_virt(CubitSimpleAttrib *csa)
{ 
  if(IsSheetBody)
  {
    TopoDS_Face * face = mySheetSurface->get_TopoDS_Face();
    OCCAttribSet::append_attribute(csa, *face);
  }
  else if(myShell)
  {
    TopoDS_Shell* shell = myShell->get_TopoDS_Shell();
    OCCAttribSet::append_attribute(csa, *shell);
  }
  else
    OCCAttribSet::append_attribute(csa, *myTopoDSShape); 
}
  
void OCCBody::remove_simple_attribute_virt(CubitSimpleAttrib *csa)
{ 
  if(IsSheetBody)
  {
    TopoDS_Face * face = mySheetSurface->get_TopoDS_Face();
    OCCAttribSet::remove_attribute(csa, *face);
  }
  else if(myShell)
  {
    TopoDS_Shell* shell = myShell->get_TopoDS_Shell();
    OCCAttribSet::remove_attribute(csa, *shell);
  }
  else
    OCCAttribSet::remove_attribute(csa, *myTopoDSShape); 
}


void OCCBody::remove_all_simple_attribute_virt()
  { remove_simple_attribute_virt(NULL); }
  
CubitStatus OCCBody::get_simple_attribute(DLIList<CubitSimpleAttrib*>& csa_list)
{ 
  if(IsSheetBody)
  {
    TopoDS_Face * face = mySheetSurface->get_TopoDS_Face();
    return OCCAttribSet::get_attributes(*face, csa_list);
  }
  else if(myShell)
  {
    TopoDS_Shell* shell = myShell->get_TopoDS_Shell();
    return OCCAttribSet::get_attributes(*shell, csa_list);
  }
  return OCCAttribSet::get_attributes(*myTopoDSShape,csa_list); 
}

CubitStatus OCCBody::get_simple_attribute( const CubitString& name,
                                          DLIList<CubitSimpleAttrib*>& csa_list )
{ 
  if(IsSheetBody)
  {
    TopoDS_Face * face = mySheetSurface->get_TopoDS_Face();
    return OCCAttribSet::get_attributes(name, *face, csa_list);
  }
  else if(myShell)
  {
    TopoDS_Shell* shell = myShell->get_TopoDS_Shell();
    return OCCAttribSet::get_attributes(name, *shell, csa_list);
  }
  return OCCAttribSet::get_attributes( name, *myTopoDSShape, csa_list ); 
}

CubitStatus OCCBody::get_transforms( CubitTransformMatrix &tfm )
{
  return CUBIT_SUCCESS;
}


//----------------------------------------------------------------
// Function: copy
// Description: create a new copy of the body.
// Author: sjowen
//----------------------------------------------------------------
BodySM* OCCBody::copy()
{
  return (BodySM*)NULL;
}
    
//----------------------------------------------------------------
// Function: move
// Description: translate the body and its child entities
//
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::move(double dx, double dy, double dz)
{
  gp_Vec aVec(dx, dy, dz);
  gp_Trsf aTrsf;
  aTrsf.SetTranslation(aVec);

  
  BRepBuilderAPI_Transform aBRepTrsf(aTrsf);
  if(IsSheetBody)
  {
    TopoDS_Face * face = mySheetSurface->get_TopoDS_Face();
    aBRepTrsf.Perform(*face);
    OCCQueryEngine::instance()->update_entity_shape(mySheetSurface, &aBRepTrsf);
  }
  else if(myShell)
  {
    TopoDS_Shell* shell = myShell->get_TopoDS_Shell();
    aBRepTrsf.Perform(*shell);
    myShell->update_OCC_entity(&aBRepTrsf);
  }
  else
  {
    aBRepTrsf.Perform(*myTopoDSShape);
    update_OCC_entity(&aBRepTrsf);
  }

  // calculate for bounding box
  update_bounding_box();

  return CUBIT_SUCCESS;
}


//----------------------------------------------------------------
// Function: rotate
// Description: rotate the body and its child entities
//
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::rotate( double x, double y, double z, 
                             double angle )//in radians
{
  gp_Pnt aOrigin(0,0,0);
  gp_Dir aDir(x, y, z);
  gp_Ax1 anAxis(aOrigin, aDir);

  //a is angular value of rotation in radians
  gp_Trsf aTrsf;
  aTrsf.SetRotation(anAxis, angle);

  BRepBuilderAPI_Transform aBRepTrsf(aTrsf);
  if(IsSheetBody)
  {
    TopoDS_Face * face = mySheetSurface->get_TopoDS_Face();
    aBRepTrsf.Perform(*face);
    OCCQueryEngine::instance()->update_entity_shape(mySheetSurface, &aBRepTrsf);
  }
  else if(myShell)
  {
    TopoDS_Shell* shell = myShell->get_TopoDS_Shell();
    aBRepTrsf.Perform(*shell);
    myShell->update_OCC_entity(&aBRepTrsf);
  }
  else
  {
    aBRepTrsf.Perform(*myTopoDSShape);
    update_OCC_entity(&aBRepTrsf);
  }

  // calculate for bounding box
  update_bounding_box();

  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Function: scale
// Description: scale the body and its child entities
//              use a constant scale factor
//
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::scale(double scale_factor )
{
  gp_Trsf aTrsf;
  aTrsf.SetScaleFactor(scale_factor);

  BRepBuilderAPI_Transform aBRepTrsf(aTrsf);
  if(IsSheetBody)
  {
    TopoDS_Face * face = mySheetSurface->get_TopoDS_Face();
    aBRepTrsf.Perform(*face);
    OCCQueryEngine::instance()->update_entity_shape(mySheetSurface, &aBRepTrsf);
  }
  else if(myShell)
  {
    TopoDS_Shell* shell = myShell->get_TopoDS_Shell();
    aBRepTrsf.Perform(*shell);
    myShell->update_OCC_entity(&aBRepTrsf);
  }
  else
  {
    aBRepTrsf.Perform(*myTopoDSShape);
    update_OCC_entity(&aBRepTrsf);
  }

  // calculate for bounding box
  update_bounding_box();

  return CUBIT_SUCCESS;  
}

//----------------------------------------------------------------
// Function: scale
// Description: scale the body and its child entities
//
// Author: 
//----------------------------------------------------------------
CubitStatus OCCBody::scale(double scale_factor_x,
                             double scale_factor_y,
                             double scale_factor_z )
{
  PRINT_ERROR("non-uniform scaling is not supported on OCC bodies");
  return CUBIT_FAILURE;
}

//----------------------------------------------------------------
// Function: reflect
// Description: reflect the body about an plane given the normal
//              vector.
//
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::reflect( double reflect_axis_x,
                              double reflect_axis_y,
                              double reflect_axis_z )
{
  gp_Pnt aOrigin(0,0,0);
  gp_Dir aDir(reflect_axis_x, reflect_axis_y,reflect_axis_z); 
  gp_Ax2 anAx2(aOrigin, aDir);

  gp_Trsf aTrsf;
  aTrsf.SetMirror(anAx2);

  BRepBuilderAPI_Transform aBRepTrsf(aTrsf);
  if(IsSheetBody)
  {
    TopoDS_Face * face = mySheetSurface->get_TopoDS_Face();
    aBRepTrsf.Perform(*face);
    OCCQueryEngine::instance()->update_entity_shape(mySheetSurface, &aBRepTrsf);
  }
  else if(myShell)
  {
    TopoDS_Shell* shell = myShell->get_TopoDS_Shell();
    aBRepTrsf.Perform(*shell);
    myShell->update_OCC_entity(&aBRepTrsf);
  }
  else
  {
    aBRepTrsf.Perform(*myTopoDSShape);
    update_OCC_entity(&aBRepTrsf);
  }

  // update underlining OCC entities
  update_bounding_box();
  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Function: private function to update the core compsolid and      
//           for any movement of the body.
// Note:     input shape must have the same number of CompSolids
//           as the body's lumps number.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::update_OCC_entity( BRepBuilderAPI_Transform *aBRepTrsf,
                                       BRepAlgoAPI_BooleanOperation *op) 
{
  if(IsSheetBody || myShell)
    return CUBIT_FAILURE;

  assert(aBRepTrsf != NULL || op != NULL);

  TopoDS_CompSolid compsolid;
  if(aBRepTrsf)
  {
    TopoDS_Shape shape = aBRepTrsf->ModifiedShape(*myTopoDSShape);
    TopoDS_CompSolid compsolid = TopoDS::CompSolid(shape);
  
    if(OCCQueryEngine::instance()->OCCMap->IsBound(*myTopoDSShape) )
       OCCQueryEngine::instance()->update_OCC_map(*myTopoDSShape, shape);
    else if (!shape.IsEqual(*myTopoDSShape))
       set_TopoDS_Shape(compsolid);
  }

  //Boolean operation works only on one lump body
  //set the lumps
  DLIList<Lump *> lumps;
  lumps = this->lumps();
  for (int i = 1; i <= lumps.size(); i++)
  {
     OCCLump *lump = CAST_TO(lumps.get_and_step(), OCCLump);
     lump->update_OCC_entity(aBRepTrsf, op);
  }

  if (aBRepTrsf && !compsolid.IsNull())
    set_TopoDS_Shape(compsolid);

  update_bounding_box(); 

  //unset marks.
  DLIList<OCCCurve*> curves;
  DLIList<OCCPoint*> points;
  get_all_curves(curves);
  get_all_points(points);

  for(int i = 0; i < curves.size(); i++)
    curves.get_and_step()->set_myMarked(CUBIT_FALSE);

  for(int i = 0; i < points.size(); i++)
    points.get_and_step()->set_myMarked(CUBIT_FALSE);
  return CUBIT_SUCCESS;
}

//----------------------------------------------------------------
// Function: TopoDS_Shape level function to update the core Body
//           for any Boolean operation of the body.
// Author: Jane Hu
//----------------------------------------------------------------
CubitStatus OCCBody::update_OCC_entity(TopoDS_Shape& old_shape,
                                       TopoDS_Shape& new_shape,
                                       BRepBuilderAPI_MakeShape *op,
                                       LocOpe_SplitShape* sp)
{
  //set the Shells
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(old_shape, TopAbs_SOLID, M);
  TopTools_ListOfShape shapes;
  TopoDS_Shape shape;

  for(int ii=1; ii<=M.Extent(); ii++)
  {
    TopoDS_Solid solid = TopoDS::Solid(M(ii));

    TopTools_ListOfShape shapes;
    if(op)
    {
      shapes.Assign(op->Modified(solid));
      if(shapes.Extent() == 0)
         shapes.Assign(op->Generated(solid));
    }
    else if(sp)
      shapes.Assign(sp->DescendantShapes(solid));

    if (shapes.Extent() == 1)
      shape = shapes.First();

    else if(shapes.Extent() > 1)
    {
      //update all attributes first.
      TopTools_ListIteratorOfListOfShape it;
      it.Initialize(shapes);
      for(; it.More(); it.Next())
      {
        shape = it.Value();
        OCCQueryEngine::instance()->copy_attributes(old_shape, shape);
      } 
      shape = shapes.First();
    }
    else if(op->IsDeleted(solid))
    {
       TopTools_IndexedMapOfShape M_new;
       TopExp::MapShapes(new_shape, TopAbs_SOLID, M_new);
       if (M_new.Extent()== 1)
         shape = M_new(1);
       else
         shape.Nullify();
    }
    else
    {
       shape = solid;
       continue;
    }

    if(shapes.Extent() > 0 || (op && op->IsDeleted(solid)))
      OCCLump::update_OCC_entity(solid, shape, op, sp);
  }
  if(!old_shape.IsSame(new_shape))
    OCCQueryEngine::instance()->update_OCC_map(old_shape, new_shape);
  return CUBIT_SUCCESS;
}
//----------------------------------------------------------------
// Function: update_bounding_box
// Description: calculate for bounding box of this OCCBody
//
// Author: janehu
//----------------------------------------------------------------
void OCCBody::update_bounding_box() 
{
     Bnd_Box box;
     TopoDS_Shape shape;
     if(IsSheetBody)
       shape = *(mySheetSurface->get_TopoDS_Face());
     else if(myShell)
       shape = *(myShell->get_TopoDS_Shell());
     else
       shape=*myTopoDSShape;

     //calculate the bounding box
     BRepBndLib::Add(shape, box);
     double min[3], max[3];

     //get values
     box.Get(min[0], min[1], min[2], max[0], max[1], max[2]);

    //update boundingbox.
    boundingbox.reset(min, max);
}

//----------------------------------------------------------------
// Function: get_bounding_box
// Description: get the  bounding box of this OCCBody
//
// Author: janehu
//----------------------------------------------------------------
CubitBox OCCBody::get_bounding_box()
{
  return boundingbox ;
}

void OCCBody::get_parents_virt( DLIList<TopologyBridge*>& ) 
  {}
  
void OCCBody::get_children_virt( DLIList<TopologyBridge*>& lumps )
{
  if (IsSheetBody)
  { 
    lumps.append(mySheetSurface->my_lump());
    return;
  }

  if (myShell)
  {
    lumps.append(myShell->my_lump());
    return;
  }

  TopoDS_Shape shape = *myTopoDSShape;
  if(OCCQueryEngine::instance()->OCCMap->IsBound(shape)==CUBIT_FALSE )
  {
    DLIList<Lump*> occ_lumps;
    occ_lumps = this->lumps();
    lumps.append( occ_lumps.get()); 
    return;
  }

  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSShape, TopAbs_SOLID, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
     TopologyBridge *lump = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
     if (lump)
       lumps.append_unique(lump);
  }
}


//-------------------------------------------------------------------------
// Purpose       : Find centroid
//
// Special Notes : 
//
// Author       : Jane Hu  
//
// Creation Date : 11/30/07
//-------------------------------------------------------------------------
CubitStatus OCCBody::mass_properties( CubitVector& centroid, 
                                      double& volume )
{
  if(IsSheetBody || myShell)
    return CUBIT_FAILURE;
  GProp_GProps myProps;
  TopoDS_Shape shape = *myTopoDSShape;
  if(OCCQueryEngine::instance()->OCCMap->IsBound(shape)==CUBIT_FALSE )
  {
    DLIList<Lump *> lumps;
    lumps = this->lumps();
    assert (lumps.size() == 1);
    OCCLump *lump = CAST_TO(lumps.get(), OCCLump);
    shape = *(lump->get_TopoDS_Solid());
  }
  BRepGProp::VolumeProperties(shape, myProps);
  volume = myProps.Mass();
  gp_Pnt pt = myProps.CentreOfMass(); 
  centroid.set(pt.X(), pt.Y(), pt.Z());
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Used to be OCCQueryEngine::is_point_in_body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitPointContainment OCCBody::point_containment( const CubitVector &point )
{
  if(IsSheetBody || myShell)
    return CUBIT_PNT_UNKNOWN;

  CubitPointContainment pc_value;
  OCCLump *lump;

  int i;
  for(i=myLumps.size(); i--;)
  {
    lump = dynamic_cast<OCCLump*>(myLumps.get_and_step()); 
    pc_value = lump->point_containment( point );
    if( pc_value == CUBIT_PNT_INSIDE )
      return CUBIT_PNT_INSIDE;
    else if( pc_value == CUBIT_PNT_BOUNDARY )
      return CUBIT_PNT_BOUNDARY;
  }

  return CUBIT_PNT_OUTSIDE;
}

//-------------------------------------------------------------------------
// Purpose       : return all surfaces in this body. 
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 01/10/08
//-------------------------------------------------------------------------

void OCCBody::get_all_surfaces(DLIList<OCCSurface*> &surfaces)
{
  if(IsSheetBody)
  {
    surfaces.append(mySheetSurface);
    return;
  }

  TopoDS_Shape shape;
  if (myShell)
    shape = *(myShell->get_TopoDS_Shell());
  else
  {
    shape = *myTopoDSShape;

    if(OCCQueryEngine::instance()->OCCMap->IsBound(shape) == CUBIT_FALSE )
    {
      DLIList<Lump *> lumps;
      lumps = this->lumps();
      assert (lumps.size() == 1);
      OCCLump *lump = CAST_TO(lumps.get(), OCCLump);
      shape = *(lump->get_TopoDS_Solid());
    }
  }
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(shape, TopAbs_FACE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
       TopologyBridge *surface = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
       OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
       if (occ_surface)
            surfaces.append_unique(occ_surface);
  }
}

//-------------------------------------------------------------------------
// Purpose       : return all curves in this body.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 01/10/08
//-------------------------------------------------------------------------

void OCCBody::get_all_curves(DLIList<OCCCurve*> &curves)
{
  TopoDS_Shape shape;
  if(IsSheetBody)
    shape = *(mySheetSurface->get_TopoDS_Face());
  else if(myShell)
    shape = *(myShell->get_TopoDS_Shell());
  else 
  {
    shape = *myTopoDSShape;

    if(OCCQueryEngine::instance()->OCCMap->IsBound(shape) == CUBIT_FALSE )
    {  
      DLIList<Lump *> lumps;
      lumps = this->lumps();
      assert (lumps.size() == 1);
      OCCLump *lump = CAST_TO(lumps.get(), OCCLump);
      shape = *(lump->get_TopoDS_Solid());
    }
  }
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(shape, TopAbs_EDGE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
       TopologyBridge *curve = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
       OCCCurve* occ_curve = CAST_TO(curve, OCCCurve);
       if (occ_curve)
         curves.append_unique(occ_curve);
  }
}

//-------------------------------------------------------------------------
// Purpose       : return all points in this body.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 01/10/08
//-------------------------------------------------------------------------

void OCCBody::get_all_points(DLIList<OCCPoint*> &points)
{
  TopoDS_Shape shape;
  if(IsSheetBody)
    shape = *(mySheetSurface->get_TopoDS_Face());
  else if(myShell)
    shape = *(myShell->get_TopoDS_Shell());
  else 
  {
    shape = *myTopoDSShape;

    if(OCCQueryEngine::instance()->OCCMap->IsBound(shape) == CUBIT_FALSE )
    {
      DLIList<Lump *> lumps;
      lumps = this->lumps();
      assert (lumps.size() == 1);
      OCCLump *lump = CAST_TO(lumps.get(), OCCLump);
      shape = *(lump->get_TopoDS_Solid());
    }
  }

  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(shape, TopAbs_VERTEX, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
       TopologyBridge *vertex = OCCQueryEngine::instance()->occ_to_cgm(M(ii));
       OCCPoint* occ_point = CAST_TO(vertex, OCCPoint);
       if (occ_point)
         points.append_unique(occ_point);
  }
}

