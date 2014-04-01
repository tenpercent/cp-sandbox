//-------------------------------------------------------------------------
// Filename      : OCCLump.hpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 7/18/2000
//-------------------------------------------------------------------------

#ifndef FACET_LUMP_HPP
#define FACET_LUMP_HPP

// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "Lump.hpp"
#include "OCCAttribSet.hpp"
#include <stdio.h>
#include "TopoDS_Solid.hxx"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class TopologyEntity;
class BodySM;
class OCCAttrib;

class OCCBody;
class OCCShell;
class OCCSurface;
class OCCLoop;
class OCCCoEdge;
class OCCCurve;
class OCCPoint;
class BRepBuilderAPI_Transform;
class BRepBuilderAPI_MakeShape;
class BRepAlgoAPI_BooleanOperation;
class LocOpe_SplitShape;
// ********** END FORWARD DECLARATIONS     **********

class OCCLump : public Lump
{
public:
  
  OCCLump(TopoDS_Solid *theSolid, OCCSurface* surface = NULL, 
          OCCShell* shell = NULL);
  
  OCCSurface* my_sheet_surface() {return mySheetSurface;}
  OCCShell * my_shell() {return myShell;}
  void set_surface(OCCSurface* surface) {mySheetSurface = surface;}
  void set_shell(OCCShell* shell) {myShell = shell;}

  virtual ~OCCLump();
    //- The destructor

  void add_body(BodySM* new_body)
    {myBodyPtr = new_body;}
    
  TopoDS_Solid *get_TopoDS_Solid(){ return myTopoDSSolid; }
  void set_TopoDS_Solid(TopoDS_Solid solid); 

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
  
  virtual CubitBox bounding_box() const ;
  
  virtual GeometryQueryEngine* 
  get_geometry_query_engine() const;
    //R GeometryQueryEngine*
    //R- A pointer to the geometric modeling engine associated with
    //R- the object.
    //- This function returns a pointer to the geometric modeling engine
    //- associated with the object.
  
  virtual CubitStatus merge( GeometryEntity* /*GEPtr*/)
    {
      PRINT_ERROR("BUG: In OCCLump::merge\n"
                  "     This function should not be called at all\n"
                  "  This is a Bug -- please report it!\n");
      return CUBIT_FAILURE;
    }
  
  virtual TopologyEntity* unmerge(DLIList<RefVolume*>)
    {
      PRINT_ERROR( "BUG: In OCCLump::unmerge\n"
                   "     This function should not be called\n"
                   "  This is a Bug -- please report it!\n" );
      return (TopologyEntity*)NULL;
    }
  
  virtual double measure();
    //R double
    //R- The numeric value of the measure (its units depend on the dimension
    //R- of the RefEntity being "measured")
    //- A generic geometric extent function.
    //- Returns volume for Lump, area for Surface, length for Curve and 
    //- 1.0 for Point

  virtual void get_parents_virt( DLIList<TopologyBridge*>& parents );
  virtual void get_children_virt( DLIList<TopologyBridge*>& children );

  inline BodySM* get_body() const { return myBodyPtr; }
    
  inline void remove_body() {myBodyPtr = 0;}
 
  CubitStatus mass_properties( CubitVector& centroid, double& volume );

  CubitPointContainment point_containment( const CubitVector &point );

  CubitStatus update_OCC_entity( BRepBuilderAPI_Transform *aBRepTrsf,
                                BRepAlgoAPI_BooleanOperation *op = NULL);
  static CubitStatus update_OCC_entity(TopoDS_Solid& old_shape,
                                       TopoDS_Shape& new_shape,
                                       BRepBuilderAPI_MakeShape *op,
                                       LocOpe_SplitShape* sp = NULL);
protected: 
  
private:
  BodySM *myBodyPtr;

  TopoDS_Solid *myTopoDSSolid;

  OCCAttribSet attribSet;
    //List of OCCAttrib*'s instead of CubitSimpleAttribs 

  OCCSurface *mySheetSurface;
  OCCShell * myShell;
} ;


// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********

// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

#endif

