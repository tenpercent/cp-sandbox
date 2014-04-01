//-------------------------------------------------------------------------
// Filename      : OCCAttribSet.hpp
//
// Purpose       : Common attrib functions for MBG
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/01/03
//-------------------------------------------------------------------------

#ifndef FACET_BRIDGE_HPP
#define FACET_BRIDGE_HPP

#include <DLIList.hpp>
#include "TDF_Label.hxx"

class CubitSimpleAttrib;
class CubitString;
class TopoDS_Shape;

class OCCAttribSet 
{

  public:
  
    OCCAttribSet() {}; 
    
    ~OCCAttribSet() { ; }
    static void FindShape(TopoDS_Shape& shape,
                          TDF_Label& aLabel,
                          CubitBoolean& found);
    
    static CubitBoolean find_attribute(TDF_Label child,
                                       CubitSimpleAttrib* csa);

    static void append_attribute( CubitSimpleAttrib*, TopoDS_Shape& shape );
    
    //remove this simple attrib from all shapes. useful when it's a shared
    //feature like materials
    static void remove_attribute( CubitSimpleAttrib* );
    
    //remove this simple attrib from the shape attribs.
    static void remove_attribute(CubitSimpleAttrib*, TopoDS_Shape& shape );
 
    //remove this shape's label from the lable tree. 
    static void remove_attribute( TopoDS_Shape& shape);

    static void get_attributes(TDF_Label &lab,
                               DLIList<CubitSimpleAttrib*>& list);

    static CubitStatus get_attributes( TopoDS_Shape& shape,
                                       DLIList<CubitSimpleAttrib*>& ) ;
    
    static CubitStatus get_attributes( const CubitString& name,
                                       TopoDS_Shape& shape,
                                DLIList<CubitSimpleAttrib*>& ) ;
    
    static int attribute_count() ;
    
  private:
};

#endif
