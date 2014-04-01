// File:        OCCShapeAttributeSet.hxx
// Created:     Thur Jul 10  2008
// Author:      Jane Hu

#ifndef _OCCShapeAttributeSet_HeaderFile
#define _OCCShapeAttributeSet_HeaderFile

class TopoDS_Shape;
class TDF_Label;

#ifndef _TopTools_LocationSet_HeaderFile
#include <TopTools_LocationSet.hxx>
#endif
#ifndef _BRep_Builder_HeaderFile
#include <BRep_Builder.hxx>
#endif
#ifndef _TopTools_IndexedMapOfShape_HeaderFile
#include <TopTools_IndexedMapOfShape.hxx>
#endif
#ifndef _GeomTools_SurfaceSet_HeaderFile
#include <GeomTools_SurfaceSet.hxx>
#endif
#ifndef _GeomTools_CurveSet_HeaderFile
#include <GeomTools_CurveSet.hxx>
#endif
#ifndef _GeomTools_Curve2dSet_HeaderFile
#include <GeomTools_Curve2dSet.hxx>
#endif
#ifndef _TColStd_IndexedMapOfTransient_HeaderFile
#include <TColStd_IndexedMapOfTransient.hxx>
#endif
#ifndef _TopAbs_ShapeEnum_HeaderFile
#include <TopAbs_ShapeEnum.hxx>
#endif
#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif

//! Contains a Shape and all  its subshapes, locations <br>
//!          and geometries, and attributes. <br>
//! <br>
//  It's re-write version of BRepTools_ShapeSet
class OCCShapeAttributeSet {

public:

    void* operator new(size_t,void* anAddress) 
      {
        return anAddress;
      }
    void* operator new(size_t size) 
      { 
        return Standard::Allocate(size); 
      }
    void  operator delete(void *anAddress) 
      { 
        if (anAddress) Standard::Free((Standard_Address&)anAddress); 
      }
 // Methods PUBLIC
 // 

//! Builds an empty ShapeAttributeSet. <br>
OCCShapeAttributeSet();

OCCShapeAttributeSet(const BRep_Builder& B);

~OCCShapeAttributeSet(){Clear();} ;

//! Stores the goemetry of <S>. <br>
void AddGeometry(const TopoDS_Shape& S) ;

//! Writes the attributs of <S>  on the stream <OS> in a <br>
//!          format that can be read back by Read. <br>
void WriteAttribute(const TopoDS_Shape& S,
                    Standard_OStream& OS,
                    TDF_Label& l_attr)const;

void  ReadAttribute(TopoDS_Shape& S,
                    Standard_IStream&   IS,
                    TDF_Label& l_attr)const;

//! Stores <S> and its sub-shape. Returns the index of <S>. <br>
//!          The method AddGeometry is called on each sub-shape. <br>
Standard_Integer Add(const TopoDS_Shape& S) ;

void Write(Standard_OStream& OS)const;

void Read(Standard_IStream& IS, bool print);

void Read(TopoDS_Shape& S,Standard_IStream& IS, const int nbshapes,
          TDF_Label* l_attr = NULL) const;

//! Writes the geometry of  me  on the stream <OS> in a <br>
//!          format that can be read back by Read. <br>
void WriteGeometry(Standard_OStream& OS) const;

void ReadGeometry(Standard_IStream& IS);

//! Writes the geometry of <S>  on the stream <OS> in a <br>
//!          format that can be read back by Read. <br>
void WriteGeometry(const TopoDS_Shape& S,Standard_OStream& OS) const;

void ReadGeometry(const TopAbs_ShapeEnum T,
                  Standard_IStream&      IS,
                  TopoDS_Shape&          S);

//! Writes   on  <OS>   the shape   <S>.    Writes the <br>
//!          orientation, the index of the TShape and the index <br>
//!          of the Location. <br>
void Write(const TopoDS_Shape& S,
           Standard_OStream& OS,
           TDF_Label* l_attr = NULL) const;

//! Writes the 3d polygons <br>
//!          on the stream <OS> in a format that can <br>
//!          be read back by Read. <br>
Standard_EXPORT   void WritePolygon3D(Standard_OStream& OS,const Standard_Boolean Compact = Standard_True) const;

void ReadPolygon3D(Standard_IStream& IS);

//! Writes the polygons on triangulation <br>
//!          on the stream <OS> in a format that can <br>
//!          be read back by Read. <br>
Standard_EXPORT   void WritePolygonOnTriangulation(Standard_OStream& OS,const Standard_Boolean Compact = Standard_True) const;

void ReadPolygonOnTriangulation(Standard_IStream& IS);

//! Writes the triangulation <br>
//!          on the stream <OS> in a format that can <br>
//!          be read back by Read. <br>
Standard_EXPORT   void WriteTriangulation(Standard_OStream& OS,const Standard_Boolean Compact = Standard_True) const;

void ReadTriangulation(Standard_IStream& IS);

//! Clears the content of the set.
void  Clear();

void Check(const TopAbs_ShapeEnum T,TopoDS_Shape& S);

int NbShapes() const;

protected:
 // Methods PROTECTED
 // 


 // Fields PROTECTED
 //


private: 

 // Methods PRIVATE
 // 


 // Fields PRIVATE
 //
BRep_Builder myBuilder;
TopTools_IndexedMapOfShape myShapes;
TopTools_LocationSet myLocations;
Standard_Integer myFormatNb;
GeomTools_SurfaceSet mySurfaces;
GeomTools_CurveSet myCurves;
GeomTools_Curve2dSet myCurves2d;
TColStd_IndexedMapOfTransient myPolygons2D;
TColStd_IndexedMapOfTransient myPolygons3D;
TColStd_IndexedMapOfTransient myTriangulations;
TColStd_IndexedMapOfTransient myNodes;
};





// other Inline functions and methods (like "C++: function call" methods)
//


#endif
