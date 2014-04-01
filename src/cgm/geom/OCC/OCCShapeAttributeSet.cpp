// File:        OCCShapeAttributeSet.cxx
// Created:     Thur Jul 10  2008
// Author:      Jane Hu

#include <Standard_Stream.hxx>
#include <BRepTools.hxx>
#include "OCCShapeAttributeSet.hpp"
#include "CubitSimpleAttrib.hpp"
#include "OCCAttribSet.hpp"
//#include <Poly.hxx>
#include <TopoDS.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_HArray1OfInteger.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <BRep_GCurve.hxx>
#include <Handle_BRep_CurveOnClosedSurface.hxx>
#include <Handle_BRep_CurveOnSurface.hxx>
#include <TopAbs_ShapeEnum.hxx>
//#include <BRep_TFace.hxx>
//#include <BRep_TEdge.hxx>
//#include <BRep_TVertex.hxx>
#include <Handle_BRep_GCurve.hxx>
#include <BRep_Tool.hxx>
#include <TDF_ChildIterator.hxx>
#include <Handle_TDataStd_Shape.hxx>
#include <TopTools_LocationSet.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS_Iterator.hxx>
//#include <BRep_CurveRepresentation.hxx>
#include <Poly_Polygon3D.hxx>
//#include <BRep_Polygon3D.hxx>
//#include <BRep_PolygonOnSurface.hxx>
//#include <BRep_PolygonOnClosedSurface.hxx>
//#include <BRep_PolygonOnTriangulation.hxx>
//#include <BRep_PolygonOnClosedTriangulation.hxx>
#include <BRep_CurveOnSurface.hxx>
#include <BRep_CurveOnClosedSurface.hxx>
//#include <BRep_ListOfCurveRepresentation.hxx>
//#include <BRep_ListIteratorOfListOfCurveRepresentation.hxx>
#include <BRep_PointOnCurve.hxx>
#include <BRep_PointOnCurveOnSurface.hxx>
#include <BRep_PointOnSurface.hxx>
//#include <BRep_ListIteratorOfListOfPointRepresentation.hxx>
#include <TDF_Label.hxx>
#include <TDataStd_Shape.hxx>
#include <Handle_TDataStd_Name.hxx>
#include <TDataStd_Name.hxx>
#include <Handle_TDataStd_ExtStringArray.hxx>
#include <TDataStd_ExtStringArray.hxx>
#include <TCollection_ExtendedString.hxx>
#include <Handle_TDataStd_IntegerArray.hxx>
#include <TDataStd_IntegerArray.hxx>
#include <Handle_TDataStd_RealArray.hxx>
#include <TDataStd_RealArray.hxx>
#include <Handle_BRep_TVertex.hxx>
#include <BRep_ListIteratorOfListOfPointRepresentation.hxx>
#include <BRep_TVertex.hxx>
#include <BRep_PointRepresentation.hxx>
#include <Handle_BRep_TEdge.hxx>
#include <BRep_ListIteratorOfListOfCurveRepresentation.hxx>
#include <Handle_BRep_TFace.hxx>
#include <BRep_TFace.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_CurveRepresentation.hxx>

#include <TopoDS_Vertex.hxx>

//#include <TColgp_HArray1OfPnt.hxx>
//#include <TColgp_HArray1OfPnt2d.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <Poly_Triangulation.hxx>
#include <Poly_PolygonOnTriangulation.hxx>


#ifdef MacOS
#define strcasecmp(p,q) strcmp(p,q)
#elseif WNT
#define strcasecmp strcmp
#elseif AIX
#include <string.h>
#endif

const char* dVersion  = "CASCADE Topology V1, (c) Matra-Datavision";
const char* dVersion2 = "CASCADE Topology V2, (c) Matra-Datavision";

//=======================================================================
//function : PrintShapeEnum
//purpose  :
//=======================================================================

static void PrintShapeEnum(const TopAbs_ShapeEnum T,
                           Standard_OStream& S,
                           Standard_Boolean C)
{
  switch(T) {

  case TopAbs_VERTEX :
    if (C) S << "Ve"; else S << "VERTEX   ";
    break;

  case TopAbs_EDGE :
    if (C) S << "Ed"; else S << "EDGE     ";
    break;

  case TopAbs_WIRE :
    if (C) S << "Wi"; else S << "WIRE     ";
    break;

  case TopAbs_FACE :
    if (C) S << "Fa"; else S << "FACE     ";
    break;

  case TopAbs_SHELL :
    if (C) S << "Sh"; else S << "SHELL    ";
    break;

  case TopAbs_SOLID :
    if (C) S << "So"; else S << "SOLID    ";
    break;

  case TopAbs_COMPSOLID :
    if (C) S << "CS"; else S << "COMPSOLID";
    break;

  case TopAbs_COMPOUND :
    if (C) S << "Co"; else S << "COMPOUND ";
    break;

  case TopAbs_SHAPE :
    if (C) S << "Sp"; else S << "SHAPE";
    break;
  }
}

//=======================================================================
//function : PrintRegularity
//purpose  :
//=======================================================================

static void PrintRegularity(const GeomAbs_Shape C,
                            Standard_OStream& OS)
{
  switch (C) {

  case GeomAbs_C0 :
    OS << "C0";
    break;

  case GeomAbs_G1 :
    OS << "G1";
    break;

  case GeomAbs_C1 :
    OS << "C1";
    break;

  case GeomAbs_G2 :
    OS << "G2";
    break;

  case GeomAbs_C2 :
    OS << "C2";
    break;

  case GeomAbs_C3 :
    OS << "C3";
    break;

  case GeomAbs_CN :
    OS << "CN";
    break;

  }
}

//=======================================================================
//function : PrintOrientation
//purpose  :
//=======================================================================

static void PrintOrientation(const TopAbs_Orientation O,
                             Standard_OStream& S,
                             Standard_Boolean C)
{
  switch(O) {

  case TopAbs_FORWARD :
    if (C) S << "+"; else S << "FORWARD";
    break;

  case TopAbs_REVERSED :
    if (C) S << "-"; else S << "REVERSED";
    break;

  case TopAbs_INTERNAL :
    if (C) S << "i"; else S << "INTERNAL";
    break;

  case TopAbs_EXTERNAL :
    if (C) S << "e"; else S << "EXTERNAL";
    break;
  }
}

//=======================================================================
//function : ReadShapeEnum
//purpose  :
//=======================================================================

static TopAbs_ShapeEnum ReadShapeEnum(Standard_IStream& IS)
{
  std::string buffer;
  IS >> buffer;

  switch (buffer[0]) {

  case 'V' :
    return TopAbs_VERTEX;

  case 'E' :
    return TopAbs_EDGE;

  case 'W' :
    return TopAbs_WIRE;

  case 'F' :
    return TopAbs_FACE;

  case 'S' :
    if (buffer[1] == 'h')
      return TopAbs_SHELL;
    else
      return TopAbs_SOLID;

  case 'C' :
    if (buffer[1] == 'S')
      return TopAbs_COMPSOLID;
    else
      return TopAbs_COMPOUND;

  }
  return TopAbs_COMPOUND;
}

//=======================================================================
//function : ReadRegularity
//purpose  :
//=======================================================================

static GeomAbs_Shape ReadRegularity(Standard_IStream& IS)
{
  std::string buffer;
  IS >> buffer;
  switch (buffer[0]) {

  case 'C' :
    switch (buffer[1]) {

    case '0' :
      return GeomAbs_C0;

    case '1' :
      return GeomAbs_C1;

    case '2' :
      return GeomAbs_C2;

    case '3' :
      return GeomAbs_C3;

    case 'N' :
      return GeomAbs_CN;
    }
    break;

  case 'G' :
    switch (buffer[1]) {

    case '1' :
      return GeomAbs_G1;

    case '2' :
      return GeomAbs_G2;

    }
    break;
  }
  return GeomAbs_C0;
}

//=======================================================================
//function : OCCShapeAttributeSet
//purpose  :
//=======================================================================

OCCShapeAttributeSet::OCCShapeAttributeSet()
  :myFormatNb(1)
{
}

//=======================================================================
//function : OCCShapeAttributeSet
//purpose  :
//=======================================================================

OCCShapeAttributeSet::OCCShapeAttributeSet (const BRep_Builder& B)
  :myBuilder(B)
{
}
//=======================================================================
//function : Add
//purpose  :
//=======================================================================

Standard_Integer  OCCShapeAttributeSet::Add(const TopoDS_Shape& S)
{
  if (S.IsNull()) return 0;
  myLocations.Add(S.Location());
  TopoDS_Shape S2 = S;
  S2.Location(TopLoc_Location());
  Standard_Integer index = myShapes.FindIndex(S2);
  if (index == 0) {
    AddGeometry(S2);

    for (TopoDS_Iterator its(S2,Standard_False,Standard_False);
         its.More(); its.Next())
      Add(its.Value());
    index = myShapes.Add(S2);
  }
 
  return index;
}

//=======================================================================
//function : AddGeometry
//purpose  :
//=======================================================================

void OCCShapeAttributeSet::AddGeometry(const TopoDS_Shape& S)
{
  // Add the geometry

  if (S.ShapeType() == TopAbs_VERTEX) {

    Handle(BRep_TVertex) TV = Handle(BRep_TVertex)::DownCast(S.TShape());
    BRep_ListIteratorOfListOfPointRepresentation itrp(TV->Points());

    while (itrp.More()) {
      const Handle(BRep_PointRepresentation)& PR = itrp.Value();

      if (PR->IsPointOnCurve()) {
        myCurves.Add(PR->Curve());
      }

      else if (PR->IsPointOnCurveOnSurface()) {
        myCurves2d.Add(PR->PCurve());
        mySurfaces.Add(PR->Surface());
      }

      else if (PR->IsPointOnSurface()) {
        mySurfaces.Add(PR->Surface());
      }

      myLocations.Add(PR->Location());
      itrp.Next();
    }

  }
  else if (S.ShapeType() == TopAbs_EDGE) {

    // Add the curve geometry
    Handle(BRep_TEdge) TE = Handle(BRep_TEdge)::DownCast(S.TShape());
    BRep_ListIteratorOfListOfCurveRepresentation itrc(TE->Curves());

    while (itrc.More()) {
      const Handle(BRep_CurveRepresentation)& CR = itrc.Value();
      if (CR->IsCurve3D()) {
        if (!CR->Curve3D().IsNull()) {
          myCurves.Add(CR->Curve3D());
          myLocations.Add(CR->Location());
        }
      }
      else if (CR->IsCurveOnSurface()) {
        mySurfaces.Add(CR->Surface());
        myCurves2d.Add(CR->PCurve());
        myLocations.Add(CR->Location());
        if (CR->IsCurveOnClosedSurface())
          myCurves2d.Add(CR->PCurve2());
      }
      else if (CR->IsRegularity()) {
        mySurfaces.Add(CR->Surface());
        myLocations.Add(CR->Location());
        mySurfaces.Add(CR->Surface2());
        myLocations.Add(CR->Location2());
      }
      itrc.Next();
    }
  }

  else if (S.ShapeType() == TopAbs_FACE) {

    // Add the surface geometry
    Handle(BRep_TFace) TF = Handle(BRep_TFace)::DownCast(S.TShape());
    if (!TF->Surface().IsNull())  mySurfaces.Add(TF->Surface());

    myLocations.Add(TF->Location());
  }
}

//=======================================================================
//function : ReadAttribute
//purpose  :
//=======================================================================

void  OCCShapeAttributeSet::ReadAttribute(TopoDS_Shape& S,
                                          Standard_IStream&   IS,
                                          TDF_Label& l_attr)const
{
  std::string buffer, type, stringdata;
  DLIList<CubitString*> strings;
  DLIList<double> doubles;
  DLIList<int> ints;
  do {
    IS >> buffer; 
    std::string::size_type i = buffer.find_first_of("*");
    type = buffer.substr( 0, i );

    CubitString* string_prt = new CubitString(type.c_str());
    strings.clean_out();
    strings.append(string_prt);

    char s = ' ';
    do {
      IS >> buffer;
      i = buffer.find_first_of("*");
      if(i > 0 && i < buffer.size())
      {
        stringdata = buffer.substr( 0, i );
        CubitString* string_prt2 = new CubitString(stringdata.c_str());
        strings.append(string_prt2);
      }

      //check if next data is still string
      IS.get(); //' '
      IS.get(s); //either '*' or 'number' or a char
      IS.unget();
    } while(s!= '*' && !(s >= '0' && s <= '9'));
    
    int tmp_int;
    double  tmp_dbl; 
    doubles.clean_out();
    ints.clean_out();
    IS.get(s); //either '*' or 'number'
    while (s != '\n')
    {
      while(s != '*') // integer attributes
      {
        IS.unget(); 
        IS >> tmp_int;
        ints.append( tmp_int );
        IS.get(); //' '
        IS.get(s); //either '*' or 'number'
      }

      IS.get(); //' '
      IS.get(s); //either '*' or 'number'
      while(s != '*') // double attributes
      {
        IS.unget();
        IS >> tmp_dbl;
        doubles.append( tmp_dbl );
        IS.get(); //' '
        IS.get(s); //either '*' or 'number' 
      }

      IS.get(s); //'\n' 
    }
    
    CubitSimpleAttrib *tmp_attrib = new CubitSimpleAttrib(&strings, &doubles, &ints);

    for(int i = 0; i < strings.size(); i++)
      delete strings.get_and_step();

    OCCAttribSet::append_attribute(tmp_attrib, S);
    delete tmp_attrib;

    IS >> buffer;
  }while(buffer[0] != '*');
}

//=======================================================================
//function : WriteAttribute
//purpose  :
//=======================================================================

void  OCCShapeAttributeSet::WriteAttribute(const TopoDS_Shape& S,
                                           Standard_OStream&   OS,
                                           TDF_Label& l_attr)const
{
  if(l_attr.IsNull())
    return;

  Standard_Boolean found = Standard_False;
  TDF_Label myLabel;
  for (TDF_ChildIterator it1(l_attr,Standard_False); it1.More(); it1.Next())
  {
    //find the same shape attribute first
    myLabel = it1.Value();

    Handle_TDataStd_Shape attr_shape;
    TopoDS_Shape exsiting_shape;
    if(TDataStd_Shape::Find(myLabel, attr_shape))
      exsiting_shape = attr_shape->Get(myLabel);

    if(!exsiting_shape.IsNull())
    {
      if(exsiting_shape.IsPartner(S))
      {
        found = Standard_True; 
        break;
      }
    }
  }
  if(!found)
  {
    OS << "\n*";
    return;
  }

  if(!myLabel.HasChild())
  {
    OS << "\n*";
    return;
  }

  for (TDF_ChildIterator it2(myLabel,Standard_False); it2.More(); it2.Next())
  {
    TDF_Label child = it2.Value();
    //Write out all attributes
    Handle_TDataStd_Name attr_name;
    TCollection_ExtendedString name_string;
    if(child.FindAttribute(TDataStd_Name::GetID(), attr_name))
    {
      OS << "\n";
      OS << "CGM_ATTRIB ";
      name_string = attr_name->Get(); 
      name_string.Print(OS);
      OS << "* " ;
    }
    else
      continue;

    Handle_TDataStd_ExtStringArray attr_strings;
    if(child.FindAttribute(TDataStd_ExtStringArray::GetID(), attr_strings))
    {
      Standard_Integer i = attr_strings->Lower();
      TCollection_ExtendedString string;
      int size = attr_strings->Upper();
      for(; i <= size; i++)
      {
         string = attr_strings->Value(i);     
         string.Print(OS);
         if(i < size ) 
           OS << " ";
      }
    }
    OS << "* " ;

    Handle_TDataStd_IntegerArray attr_ints;
    
    if(child.FindAttribute(TDataStd_IntegerArray::GetID(), attr_ints))
    {
      for(Standard_Integer i = attr_ints->Lower(); i <= attr_ints->Upper(); i++)
        OS << attr_ints->Value(i) << " ";
    }
    
    OS << '*' << ' ' ;

    Handle_TDataStd_RealArray attr_doubles;
    if(child.FindAttribute(TDataStd_RealArray::GetID(), attr_doubles))
    {
      Standard_Integer i = attr_doubles->Lower();
      for(;i <= attr_doubles->Upper(); i++)
        OS << attr_doubles->Value(i) << " ";
    }
    OS << '*' ; 
  }
  OS << "\n*";
}

//=======================================================================
//function : Write
//purpose  :
//=======================================================================

void  OCCShapeAttributeSet::Write(Standard_OStream& OS)const
{
  //on sauvegarde l'ancien LC_NUMERIC
  
  char *oldnum,*plocal ;
  plocal =   setlocale(LC_NUMERIC, NULL) ;
  oldnum = new char[strlen(plocal)+1] ;
  strcpy(oldnum,plocal);

  // on positionne LC_NUMERIC a "C" (point decimal)
  setlocale(LC_NUMERIC, "C") ;

  int  prec = OS.precision(15);

  // write the copyright
  if (myFormatNb == 2)
    OS << "\n" << dVersion2 << endl;
  else
    OS << "\n" << dVersion << endl;

  //-----------------------------------------
  // write the locations
  //-----------------------------------------
  myLocations.Write(OS);

  //-----------------------------------------
  // write the geometry
  //-----------------------------------------

  WriteGeometry(OS);

  //-----------------------------------------
  // write the shapes
  //-----------------------------------------
  Standard_Integer i, nbShapes = myShapes.Extent();

  OS << "\nTShapes " << nbShapes << "\n";


  // subshapes are written first
  for (i = 1; i <= nbShapes; i++) {

    const TopoDS_Shape& S = myShapes(i);

    // Type
    PrintShapeEnum(S.ShapeType(),OS,Standard_True);
    OS << "\n";

    // Geometry
    WriteGeometry(S,OS);

    // Flags
    OS << "\n";
    OS << (S.Free()       ? 1 : 0);
    OS << (S.Modified()   ? 1 : 0);
    OS << (S.Checked()    ? 1 : 0);
    OS << (S.Orientable() ? 1 : 0);
    OS << (S.Closed()     ? 1 : 0);
    OS << (S.Infinite()   ? 1 : 0);
    OS << (S.Convex()     ? 1 : 0);
    OS << "\n";

    // sub-shapes

    Standard_Integer l = 0;
    TopoDS_Iterator its(S,Standard_False,Standard_False);
    while (its.More()) {
      Write(its.Value(),OS);
      l++;
      if (l == 10) {
        OS << "\n";
        l = 0;
      }
      its.Next();
    }
    Write(TopoDS_Shape(),OS); // Null shape to end the list
    OS << "\n";
  }

  OS << endl;
  OS.precision(prec);

  // on remet le LC_NUMERIC a la precedente valeur
  setlocale(LC_NUMERIC, oldnum) ;
  delete[] oldnum;
}

//=======================================================================
//function : Read
//purpose  :
//=======================================================================

void  OCCShapeAttributeSet::Read(Standard_IStream& IS,
                                 bool print_results) 
{
 // on sauvegarde l'ancien LC_NUMERIC
  char *oldnum,*plocal ;
  plocal =   setlocale(LC_NUMERIC, NULL) ;
  oldnum = new char[strlen(plocal)+1] ;
  strcpy(oldnum,plocal);

  Clear();

  // Check the version
  char vers[101];
  do {
    IS.getline(vers,100,'\n');
    // BUC60769 PTV 18.10.2000: remove possible '\r' at the end of the line
    //Standard_Integer lv = strlen(vers);
    //char *pm;
    //if(pm = strchr(vers,'\r'))
    //  *pm ='\0';

    for (Standard_Integer lv = (strlen(vers)- 1); lv > 1 && (vers[lv] == '\r' || vers[lv] == '\n') ;lv--)
      vers[lv] = '\0';

  } while ( ! IS.fail() && strcmp(vers,dVersion) && strcmp(vers,dVersion2) );
  if (IS.fail()) {
    if (print_results)
      cout << "File was not written with this version of the topology"<<endl;
    setlocale(LC_NUMERIC, oldnum) ;
    delete[] oldnum;
    return;
  }
  if (strcmp(vers,dVersion2) == 0) myFormatNb = 2;
  else myFormatNb = 1;

  //-----------------------------------------
  // read the locations
  //-----------------------------------------

  myLocations.Read(IS);

  //-----------------------------------------
  // read the geometry
  //-----------------------------------------

  ReadGeometry(IS);

  //-----------------------------------------
  // read the shapes
  //-----------------------------------------

  std::string buffer;
  IS >> buffer;
  if (buffer != "TShapes") {
    if (print_results)
      cout << "Not a TShape table"<<endl;
    setlocale(LC_NUMERIC, oldnum) ;
    delete[] oldnum;
    return;
  }

  Standard_Integer i, nbShapes;
  IS >> nbShapes;

  for (i = 1; i <= nbShapes; i++) {

    TopoDS_Shape S;

    //Read type and create empty shape.
    TopAbs_ShapeEnum T = ReadShapeEnum(IS);
    ReadGeometry(T,IS,S);

    // Set the flags
    IS >> buffer;

    // sub-shapes
    TopoDS_Shape SS;
    do {
      Read(SS,IS,nbShapes);
      if (!SS.IsNull())
        myBuilder.Add(S,SS);
    } while(!SS.IsNull());

    S.Free      (buffer[0] == '1');
    S.Modified  (buffer[1] == '1');

    if (myFormatNb == 2)
      S.Checked   (buffer[2] == '1');
    else
      S.Checked   (Standard_False);     // force check at reading..

    S.Orientable(buffer[3] == '1');
    S.Closed    (buffer[4] == '1');
    S.Infinite  (buffer[5] == '1');
    S.Convex    (buffer[6] == '1');

    // check

    if (myFormatNb == 1)
      Check(T,S);

    myShapes.Add(S);
  }

  setlocale(LC_NUMERIC, oldnum) ;
  delete[] oldnum;
}

//=======================================================================
//function : WriteGeometry
//purpose  :
//=======================================================================

void  OCCShapeAttributeSet::WriteGeometry(Standard_OStream& OS) const
{
  myCurves2d.Write(OS);
  myCurves.Write(OS);
  WritePolygon3D(OS);
  WritePolygonOnTriangulation(OS);
  mySurfaces.Write(OS);
  WriteTriangulation(OS);
}

//=======================================================================
//function : ReadGeometry
//purpose  :
//=======================================================================

void  OCCShapeAttributeSet::ReadGeometry(Standard_IStream& IS)
{
  myCurves2d.Read(IS);
  myCurves.Read(IS);
  ReadPolygon3D(IS);
  ReadPolygonOnTriangulation(IS);
  mySurfaces.Read(IS);
  ReadTriangulation(IS);
}

//=======================================================================
//function : WritePolygon3D
//purpose  :
//=======================================================================

void OCCShapeAttributeSet::WritePolygon3D(Standard_OStream&      OS,
                                        const Standard_Boolean Compact)const
{
  Standard_Integer i, j, nbpol = myPolygons3D.Extent();
  if (Compact)
    OS << "Polygon3D " << nbpol << endl;
  else {
    OS << " -------\n";
    OS <<"Dump of " << nbpol << " Polygon3Ds\n";
    OS << " -------\n";
  }

  Handle(Poly_Polygon3D) P;
  for (i = 1; i <= nbpol; i++) {
    P = Handle(Poly_Polygon3D)::DownCast(myPolygons3D(i));
    if (Compact) {
      OS << P->NbNodes() << " ";
      OS << ((P->HasParameters()) ? "1" : "0") << "\n";
    }
    else {
      OS << "  "<< i << " : Polygon3D with " << P->NbNodes() << " Nodes\n";
      OS << ((P->HasParameters()) ? "with" : "without") << " parameters\n";
    }


    // write the deflection
    if (!Compact) OS << "Deflection : ";
    OS << P->Deflection() << "\n";

    // write the nodes
    if (!Compact) OS << "\nNodes :\n";

    Standard_Integer i1, nbNodes = P->NbNodes();
    const TColgp_Array1OfPnt& Nodes = P->Nodes();
    for (j = 1; j <= nbNodes; j++) {
      if (!Compact) OS << setw(10) << j << " : ";
      if (!Compact) OS << setw(17);
      OS << Nodes(j).X() << " ";
      if (!Compact) OS << setw(17);
      OS << Nodes(j).Y() << " ";
      if (!Compact) OS << setw(17);
      OS << Nodes(j).Z();
      if (!Compact) OS << "\n";
      else OS << " ";
    }
    OS <<"\n";

    if (P->HasParameters()) {
      if (!Compact) OS << "\nParameters :\n";
      const TColStd_Array1OfReal& Param = P->Parameters();
      for ( i1 = 1; i1 <= nbNodes; i1++ ) {
        OS << Param(i1) << " ";
      }
      OS <<"\n";
    }
  }
}

//=======================================================================
//function : WritePolygonOnTriangulation
//purpose  :
//=======================================================================

void OCCShapeAttributeSet::WritePolygonOnTriangulation(
                                          Standard_OStream&      OS,
                                          const Standard_Boolean Compact)const
{
  Standard_Integer i, j, nbpOntri = myNodes.Extent();
  if (Compact)
    OS << "PolygonOnTriangulations " << nbpOntri << endl;
  else {
    OS << " -------\n";
    OS <<"Dump of " << nbpOntri << " PolygonOnTriangulations\n";
    OS << " -------\n";
  }

  Handle(Poly_PolygonOnTriangulation) Poly;
  Handle(TColStd_HArray1OfReal) Param;

  for (i=1; i<=nbpOntri; i++) {
    Poly = Handle(Poly_PolygonOnTriangulation)::DownCast(myNodes(i));
    const TColStd_Array1OfInteger& Nodes = Poly->Nodes();
    if (!Compact) {
      OS << "  "<< i << " : PolygonOnTriangulation with " << Nodes.Length() << " Nodes\n";
    }
    else OS << Nodes.Length()<<" ";
    if (!Compact) OS <<"  ";
    for (j=1; j <= Nodes.Length(); j++) OS << Nodes.Value(j) << " ";
    OS << "\n";

    // writing parameters:
    Param = Poly->Parameters();
    if (Compact) OS <<"p ";

    // write the deflection
    if (!Compact) OS << "  Deflection : ";
    OS <<Poly->Deflection() << " ";
    if (!Compact) OS << "\n";

    if (!Param.IsNull()) {
      if (!Compact) {
        OS << "  "<< "Parameters :";
      }
      else OS << "1 " ;
      if (!Compact) OS <<"  ";
      for (j=1; j <= Param->Length(); j++) OS << Param->Value(j) << " ";
      OS << "\n";
    }
    else OS <<"0 \n";
  }

}

//=======================================================================
//function : WriteTriangulation
//purpose  :
//=======================================================================

void OCCShapeAttributeSet::WriteTriangulation(Standard_OStream&      OS,
                                        const Standard_Boolean Compact)const
{
  Standard_Integer i, j, nbNodes, nbtri = myTriangulations.Extent();
  Standard_Integer nbTriangles = 0, n1, n2, n3;
  if (Compact)
    OS << "Triangulations " << nbtri << endl;
  else {
    OS << " -------\n";
    OS <<"Dump of " << nbtri << " Triangulations\n";
    OS << " -------\n";
  }

  Handle(Poly_Triangulation) T;
  for (i = 1; i <= nbtri; i++) {
    T = Handle(Poly_Triangulation)::DownCast(myTriangulations(i));
    if (Compact) {
      OS << T->NbNodes() << " " << T->NbTriangles() << " ";
      OS << ((T->HasUVNodes()) ? "1" : "0") << " ";
    }
    else {
      OS << "  "<< i << " : Triangulation with " << T->NbNodes() << " Nodes and "
         << T->NbTriangles() <<" Triangles\n";
      OS << "      "<<((T->HasUVNodes()) ? "with" : "without") << " UV nodes\n";
    }

    // write the deflection

    if (!Compact) OS << "  Deflection : ";
    OS <<T->Deflection() << "\n";

    // write the 3d nodes

    if (!Compact) OS << "\n3D Nodes :\n";

    nbNodes = T->NbNodes();
    const TColgp_Array1OfPnt& Nodes = T->Nodes();
    for (j = 1; j <= nbNodes; j++) {
      if (!Compact) OS << setw(10) << j << " : ";
      if (!Compact) OS << setw(17);
      OS << Nodes(j).X() << " ";
      if (!Compact) OS << setw(17);
      OS << Nodes(j).Y() << " ";
      if (!Compact) OS << setw(17);
      OS << Nodes(j).Z();
      if (!Compact) OS << "\n";
      else OS << " ";
    }

    if (T->HasUVNodes()) {
      if (!Compact) OS << "\nUV Nodes :\n";
      const TColgp_Array1OfPnt2d& UVNodes = T->UVNodes();
      for (j = 1; j <= nbNodes; j++) {
        if (!Compact) OS << setw(10) << j << " : ";
        if (!Compact) OS << setw(17);
        OS << UVNodes(j).X() << " ";
        if (!Compact) OS << setw(17);
        OS << UVNodes(j).Y();
        if (!Compact) OS << "\n";
        else OS << " ";
      }
    }

    if (!Compact) OS << "\nTriangles :\n";
    nbTriangles = T->NbTriangles();
    const Poly_Array1OfTriangle& Triangles = T->Triangles();
    for (j = 1; j <= nbTriangles; j++) {
      if (!Compact) OS << setw(10) << j << " : ";
      Triangles(j).Get(n1, n2, n3);
      if (!Compact) OS << setw(10);
      OS << n1 << " ";
      if (!Compact) OS << setw(10);
      OS << n2 << " ";
      if (!Compact) OS << setw(10);
      OS << n3;
      if (!Compact) OS << "\n";
      else OS << " ";
    }
    OS << "\n";
  }

}
//=======================================================================
//function : WriteGeometry
//purpose  :
//=======================================================================

void  OCCShapeAttributeSet::WriteGeometry(const TopoDS_Shape& S,
                                        Standard_OStream&   OS)const
{
  // Write the geometry

  if (S.ShapeType() == TopAbs_VERTEX) {

    // Write the point geometry
    TopoDS_Vertex V = TopoDS::Vertex(S);
    OS << BRep_Tool::Tolerance(V) << "\n";
    gp_Pnt p = BRep_Tool::Pnt(V);
    OS<<p.X()<<" "<<p.Y()<<" "<<p.Z()<<"\n";

    Handle(BRep_TVertex) TV = Handle(BRep_TVertex)::DownCast(S.TShape());
    BRep_ListIteratorOfListOfPointRepresentation itrp(TV->Points());

    while (itrp.More()) {
      const Handle(BRep_PointRepresentation)& PR = itrp.Value();

      OS << PR->Parameter();
      if (PR->IsPointOnCurve()) {
        OS << " 1 " << myCurves.Index(PR->Curve());
      }

      else if (PR->IsPointOnCurveOnSurface()) {
        OS << " 2 " <<  myCurves2d.Index(PR->PCurve());
        OS << " " << mySurfaces.Index(PR->Surface());
      }

      else if (PR->IsPointOnSurface()) {
        OS << " 3 " << PR->Parameter2() << " ";
        OS << mySurfaces.Index(PR->Surface());
      }

      OS << " " << myLocations.Index(PR->Location());
      OS << "\n";

      itrp.Next();
    }

    OS << "0 0\n"; // end representations

  }

  else if (S.ShapeType() == TopAbs_EDGE) {

    // Write the curve geometry

    Handle(BRep_TEdge) TE = Handle(BRep_TEdge)::DownCast(S.TShape());

    OS << " " << TE->Tolerance() << " ";
    OS << ((TE->SameParameter()) ? 1 : 0) << " ";
    OS << ((TE->SameRange())     ? 1 : 0) << " ";
    OS << ((TE->Degenerated())   ? 1 : 0) << "\n";

    Standard_Real first, last;
    BRep_ListIteratorOfListOfCurveRepresentation itrc = TE->Curves();
    while (itrc.More()) {
      const Handle(BRep_CurveRepresentation)& CR = itrc.Value();
      if (CR->IsCurve3D()) {
        if (!CR->Curve3D().IsNull()) {
          Handle(BRep_GCurve) GC = Handle(BRep_GCurve)::DownCast(itrc.Value());
          GC->Range(first, last);
          OS << "1 ";                               // -1- Curve 3D
          OS << " "<<myCurves.Index(CR->Curve3D());
          OS << " "<<myLocations.Index(CR->Location());
          OS << " "<<first<<" "<<last;
          OS << "\n";
        }
      }
      else if (CR->IsCurveOnSurface()) {
        Handle(BRep_GCurve) GC = Handle(BRep_GCurve)::DownCast(itrc.Value());
        GC->Range(first, last);
        if (!CR->IsCurveOnClosedSurface())
          OS << "2 ";                             // -2- Curve on surf
        else
          OS << "3 ";                             // -3- Curve on closed surf
        OS <<" "<<myCurves2d.Index(CR->PCurve());
        if (CR->IsCurveOnClosedSurface()) {
          OS <<" " << myCurves2d.Index(CR->PCurve2());
          PrintRegularity(CR->Continuity(),OS);
        }
        OS << " " << mySurfaces.Index(CR->Surface());
        OS << " " << myLocations.Index(CR->Location());
        OS << " "<<first<<" "<<last;
        OS << "\n";

        // Write UV Points // for XML Persistence higher performance
        if (myFormatNb == 2)
        {
          gp_Pnt2d Pf,Pl;
          if (CR->IsCurveOnClosedSurface()) {
            Handle(BRep_CurveOnClosedSurface) COCS =
              Handle(BRep_CurveOnClosedSurface)::DownCast(CR);
            COCS->UVPoints2(Pf,Pl);
          }
          else {
            Handle(BRep_CurveOnSurface) COS =
              Handle(BRep_CurveOnSurface)::DownCast(CR);
            COS->UVPoints(Pf,Pl);
          }
          OS << Pf.X() << " " << Pf.Y() << " " << Pl.X() << " " << Pl.Y() << "\n";
        }
      }
      else if (CR->IsRegularity()) {
        OS << "4 ";                              // -4- Regularity
        PrintRegularity(CR->Continuity(),OS);
        OS << " "<<mySurfaces.Index(CR->Surface());
        OS << " "<<myLocations.Index(CR->Location());
        OS << " "<<mySurfaces.Index(CR->Surface2());
        OS << " "<<myLocations.Index(CR->Location2());
        OS << "\n";
      }
      itrc.Next();
    }
    OS << "0\n"; // end of the list of representations
  }
   else if (S.ShapeType() == TopAbs_FACE) {

    Handle(BRep_TFace) TF = Handle(BRep_TFace)::DownCast(S.TShape());
    const TopoDS_Face& F = TopoDS::Face(S);

    if (!(TF->Surface()).IsNull()) {
      OS << ((BRep_Tool::NaturalRestriction(F)) ? 1 : 0);
      OS << " ";
      // Write the surface geometry
      OS << " " <<TF->Tolerance();
      OS << " " <<mySurfaces.Index(TF->Surface());
      OS << " " <<myLocations.Index(TF->Location());
      OS << "\n";
    }
  }

}


//=======================================================================
//function : Write
//purpose  :
//=======================================================================

void  OCCShapeAttributeSet::Write(const TopoDS_Shape& S,
                               Standard_OStream& OS,
                               TDF_Label* l_attr)const
{
  if (S.IsNull()) OS << "*";
  else {
    PrintOrientation(S.Orientation(),OS,Standard_True);
    OS << myShapes.Extent() - myShapes.FindIndex(S.Located(TopLoc_Location())) + 1;
    OS << " " << myLocations.Index(S.Location()) << " ";
  }
  //Write Attributes
  Standard_Integer i, nbShapes = myShapes.Extent();
  if(l_attr != NULL)
  {
    for ( i = 1; i <= nbShapes; i++)
    {
      const TopoDS_Shape& Sh = myShapes(i);
      WriteAttribute(Sh, OS, *l_attr);
    }
  }
}


//=======================================================================
//function : ReadGeometry
//purpose  :
//=======================================================================

void  OCCShapeAttributeSet::ReadGeometry(const TopAbs_ShapeEnum T,
                                         Standard_IStream&      IS,
                                         TopoDS_Shape&          S)
{
  // Read the geometry

  Standard_Integer val,c,pc,pc2,s,s2,l,l2,t, pt, pt2;
  Standard_Real tol,X,Y,Z,first,last,p1,p2;
  Standard_Real PfX,PfY,PlX,PlY;
  gp_Pnt2d aPf, aPl;
  Standard_Boolean closed;
#ifndef DEB
  GeomAbs_Shape reg = GeomAbs_C0;
#else
  GeomAbs_Shape reg;
#endif
  switch (T) {


    //---------
    // vertex
    //---------

  case TopAbs_VERTEX :
    {
      TopoDS_Vertex& V = TopoDS::Vertex(S);

      // Read the point geometry
      IS >> tol;
      IS >> X >> Y >> Z;
      myBuilder.MakeVertex(V,gp_Pnt(X,Y,Z),tol);
      Handle(BRep_TVertex) TV = Handle(BRep_TVertex)::DownCast(V.TShape());

      BRep_ListOfPointRepresentation& lpr = TV->ChangePoints();
      TopLoc_Location L;

      do {
        IS >> p1 >> val;

        Handle(BRep_PointRepresentation) PR;
        switch (val) {

        case 1 :
          {
            IS >> c;

            if (myCurves.Curve(c).IsNull())
              break;

            Handle(BRep_PointOnCurve) POC =
              new BRep_PointOnCurve(p1,
                                    myCurves.Curve(c),
                                    L);
            PR = POC;
          }
          break;

        case 2 :
          {
            IS >> pc >> s;

            if (myCurves2d.Curve2d(pc).IsNull() ||
                mySurfaces.Surface(s).IsNull())
              break;

            Handle(BRep_PointOnCurveOnSurface) POC =
              new BRep_PointOnCurveOnSurface(p1,
                                             myCurves2d.Curve2d(pc),
                                             mySurfaces.Surface(s),
                                             L);
            PR = POC;
          }
          break;

        case 3 :
          {
            IS >> p2 >> s;

            if (mySurfaces.Surface(s).IsNull())
              break;

            Handle(BRep_PointOnSurface) POC =
              new BRep_PointOnSurface(p1,p2,
                                      mySurfaces.Surface(s),
                                      L);
            PR = POC;
          }
          break;
        }

        if (val > 0) {
          IS >> l;
          if (!PR.IsNull()) {
            PR->Location(myLocations.Location(l));
            lpr.Append(PR);
          }
        }
      } while (val > 0);
    }
    break;


    //---------
    // edge
    //---------


    case TopAbs_EDGE :

      // Create an edge
      {
        TopoDS_Edge& E = TopoDS::Edge(S);

        myBuilder.MakeEdge(E);

        // Read the curve geometry
        IS >> tol;
        IS >> val;
        myBuilder.SameParameter(E,(val == 1));
        IS >> val;
        myBuilder.SameRange(E,(val == 1));
        IS >> val;
        myBuilder.Degenerated(E,(val == 1));

        do {
          IS >> val;
          switch (val) {

          case 1 :                               // -1- Curve 3D
            IS >> c >> l;
            if (!myCurves.Curve(c).IsNull()) {
              myBuilder.UpdateEdge(E,myCurves.Curve(c),
                                   myLocations.Location(l),tol);
            }
            IS >> first >> last;
            if (!myCurves.Curve(c).IsNull()) {
              Standard_Boolean Only3d = Standard_True;
              myBuilder.Range(E,first,last,Only3d);
            }
            break;


          case 2 :                               // -2- Curve on surf
          case 3 :                               // -3- Curve on closed surf
            closed = (val == 3);
            IS >> pc;
            if (closed) {
              IS >> pc2;
              reg = ReadRegularity(IS);
            }

            // surface, location
            IS >> s >> l;

            // range
            IS >> first >> last;

            // read UV Points // for XML Persistence higher performance
            if (myFormatNb == 2)
            {
              IS >> PfX >> PfY >> PlX >> PlY;
              aPf = gp_Pnt2d(PfX,PfY);
              aPl = gp_Pnt2d(PlX,PlY);
            }

            if (myCurves2d.Curve2d(pc).IsNull() ||
                (closed && myCurves2d.Curve2d(pc2).IsNull()) ||
                mySurfaces.Surface(s).IsNull())
              break;

            if (closed) {
              if (myFormatNb == 2)
                myBuilder.UpdateEdge(E,myCurves2d.Curve2d(pc),
                                     myCurves2d.Curve2d(pc2),
                                     mySurfaces.Surface(s),
                                     myLocations.Location(l),tol,
                                     aPf, aPl);
              else
                myBuilder.UpdateEdge(E,myCurves2d.Curve2d(pc),
                                     myCurves2d.Curve2d(pc2),
                                     mySurfaces.Surface(s),
                                     myLocations.Location(l),tol);

              myBuilder.Continuity(E,
                                   mySurfaces.Surface(s),
                                   mySurfaces.Surface(s),
                                   myLocations.Location(l),
                                   myLocations.Location(l),
                                   reg);
            }
            else
            {
              if (myFormatNb == 2)
                myBuilder.UpdateEdge(E,myCurves2d.Curve2d(pc),
                                     mySurfaces.Surface(s),
                                     myLocations.Location(l),tol,
                                     aPf, aPl);
              else
                myBuilder.UpdateEdge(E,myCurves2d.Curve2d(pc),
                                     mySurfaces.Surface(s),
                                     myLocations.Location(l),tol);
            }
            myBuilder.Range(E,
                            mySurfaces.Surface(s),
                            myLocations.Location(l),
                            first,last);
            break;

          case 4 :                               // -4- Regularity
            reg = ReadRegularity(IS);
            IS >> s >> l >> s2 >> l2;
            if (mySurfaces.Surface(s).IsNull() ||
                mySurfaces.Surface(s2).IsNull())
              break;
            myBuilder.Continuity(E,
                                 mySurfaces.Surface(s),
                                 mySurfaces.Surface(s2),
                                 myLocations.Location(l),
                                 myLocations.Location(l2),
                                 reg);
            break;

          case 5 :   // -5- Polygon3D
            IS >> c >> l;
            myBuilder.UpdateEdge(E,Handle(Poly_Polygon3D)::DownCast(myPolygons3D(c)), myLocations.Location(l));
            break;

          case 6 :
          case 7 :
            closed = (val == 7);
            IS >> pt;
            if (closed) {
              IS >> pt2;
            }
            IS >> t >> l;
            if (closed) {
              myBuilder.UpdateEdge
                (E, Handle(Poly_PolygonOnTriangulation)::DownCast(myNodes(pt)),
                 Handle(Poly_PolygonOnTriangulation)::DownCast(myNodes(pt2)),
                 Handle(Poly_Triangulation)::DownCast(myTriangulations(t)),
                 myLocations.Location(l));
            }
            else {
              myBuilder.UpdateEdge
                (E,Handle(Poly_PolygonOnTriangulation)::DownCast(myNodes(pt)),
                 Handle(Poly_Triangulation)::DownCast(myTriangulations(t)),
                 myLocations.Location(l));
            }
            // range

            break;

          }
        } while (val > 0);
      }
    break;


    //---------
    // wire
    //---------

  case TopAbs_WIRE :
    myBuilder.MakeWire(TopoDS::Wire(S));
    break;


    //---------
    // face
    //---------
  case TopAbs_FACE :
    {
    // create a face :
    TopoDS_Face& F = TopoDS::Face(S);
    myBuilder.MakeFace(F);

    IS >> val; // natural restriction
    if (val == 0 || val == 1) {
      IS >> tol >> s >> l;
      if (!mySurfaces.Surface(s).IsNull()) {
        myBuilder.UpdateFace(TopoDS::Face(S),
                             mySurfaces.Surface(s),
                             myLocations.Location(l),tol);
        myBuilder.NaturalRestriction(TopoDS::Face(S),(val == 1));
      }
    }

    // BUC60769
    std::string line;
    std::getline( IS, line );
    std::getline( IS, line );
    std::istringstream str( line );

    if (str.get() == '2') {
      // cas triangulation
      str >> s;
      myBuilder.UpdateFace(TopoDS::Face(S),
                           Handle(Poly_Triangulation)::DownCast(myTriangulations(s)));
    }
//    else IS.seekg(pos);
    }
    break;


    //---------
    // shell
    //---------

  case TopAbs_SHELL :
    myBuilder.MakeShell(TopoDS::Shell(S));
    break;


    //---------
    // solid
    //---------

  case TopAbs_SOLID :
    myBuilder.MakeSolid(TopoDS::Solid(S));
    break;


    //---------
    // compsolid
    //---------

  case TopAbs_COMPSOLID :
    myBuilder.MakeCompSolid(TopoDS::CompSolid(S));
    break;


    //---------
    // compound
    //---------

  case TopAbs_COMPOUND :
    myBuilder.MakeCompound(TopoDS::Compound(S));
    break;

  default:
    break;
  }

}

//=======================================================================
//function : ReadPolygonOnTriangulation
//purpose  :
//=======================================================================

void OCCShapeAttributeSet::ReadPolygonOnTriangulation(Standard_IStream& IS)
{
    std::string buffer;
    IS >> buffer;
    if (buffer.find("PolygonOnTriangulations") == std::string::npos)
      return;

    Standard_Integer i, j, val, nbpol = 0, nbnodes =0;
    Standard_Integer hasparameters;
    Standard_Real par;
    Handle(TColStd_HArray1OfReal) Param;
    Handle(Poly_PolygonOnTriangulation) Poly;
    IS >> nbpol;
    for (i=1; i<=nbpol; i++) {
      IS >> nbnodes;
      TColStd_Array1OfInteger Nodes(1, nbnodes);
      for (j = 1; j <= nbnodes; j++) {
        IS >> val;
        Nodes(j) = val;
      }
      IS >> buffer;
      Standard_Real def;
      IS >> def;
      IS >> hasparameters;
      if (hasparameters) {
        TColStd_Array1OfReal Param1(1, nbnodes);
        for (j = 1; j <= nbnodes; j++) {
          IS >> par;
          Param1(j) = par;
        }
        Poly = new Poly_PolygonOnTriangulation(Nodes, Param1);
      }
      else Poly = new Poly_PolygonOnTriangulation(Nodes);
      Poly->Deflection(def);
      myNodes.Add(Poly);
    }
}

//=======================================================================
//function : ReadPolygon3D
//purpose  :
//=======================================================================

void OCCShapeAttributeSet::ReadPolygon3D(Standard_IStream& IS)
{
    std::string buffer;
    //  Standard_Integer i, j, p, val, nbpol, nbnodes, hasparameters;
    Standard_Integer i, j, p, nbpol=0, nbnodes =0, hasparameters = Standard_False;  Standard_Real d, x, y, z;

    IS >> buffer;
    if (buffer.find("Polygon3D") == std::string::npos)
      return;
  
    Handle(Poly_Polygon3D) P;
    IS >> nbpol;
    for (i=1; i<=nbpol; i++) {
      IS >> nbnodes;
      IS >> hasparameters;
      TColgp_Array1OfPnt Nodes(1, nbnodes);
      IS >> d;
      for (j = 1; j <= nbnodes; j++) {
        IS >> x >> y >> z;
        Nodes(j).SetCoord(x,y,z);
      }
      if (hasparameters) {
        TColStd_Array1OfReal Param(1,nbnodes);
        for (p = 1; p <= nbnodes; p++) {
          IS >> Param(p);
        }
        P = new Poly_Polygon3D(Nodes, Param);
      }
      else P = new Poly_Polygon3D(Nodes);
      P->Deflection(d);
      myPolygons3D.Add(P);
    }
}

//=======================================================================
//function : ReadTriangulation
//purpose  :
//=======================================================================

void OCCShapeAttributeSet::ReadTriangulation(Standard_IStream& IS)
{
  std::string buffer;
  //  Standard_Integer i, j, val, nbtri;
  Standard_Integer i, j, nbtri =0;
  Standard_Real d, x, y, z;
  Standard_Integer nbNodes =0, nbTriangles=0;
  Standard_Boolean hasUV= Standard_False;

  Handle(Poly_Triangulation) T;

  IS >> buffer;
  if (buffer.find("Triangulations") != std::string::npos) {
    IS >> nbtri;
    for (i=1; i<=nbtri; i++) {
      IS >> nbNodes >> nbTriangles >> hasUV;
      IS >> d;

      TColgp_Array1OfPnt Nodes(1, nbNodes);
      TColgp_Array1OfPnt2d UVNodes(1, nbNodes);

      for (j = 1; j <= nbNodes; j++) {
        IS >> x >> y >> z;
        Nodes(j).SetCoord(x,y,z);
      }

      if (hasUV) {
        for (j = 1; j <= nbNodes; j++) {
          IS >> x >> y;
          UVNodes(j).SetCoord(x,y);
        }
      }


      // read the triangles
      Standard_Integer n1,n2,n3;
      Poly_Array1OfTriangle Triangles(1, nbTriangles);
      for (j = 1; j <= nbTriangles; j++) {
        IS >> n1 >> n2 >> n3;
        Triangles(j).Set(n1,n2,n3);
      }

      if (hasUV) T =  new Poly_Triangulation(Nodes,UVNodes,Triangles);
      else T = new Poly_Triangulation(Nodes,Triangles);

      T->Deflection(d);

      myTriangulations.Add(T);
    }
  }
}

//=======================================================================
//function : Clear
//purpose  :
//=======================================================================

void  OCCShapeAttributeSet::Clear()
{
  mySurfaces.Clear();
  myCurves.Clear();
  myCurves2d.Clear();
  myPolygons3D.Clear();
  myPolygons2D.Clear();
  myNodes.Clear();
  myTriangulations.Clear();
  myShapes.Clear();
  myLocations.Clear();
}

//=======================================================================
//function : Read
//purpose  :
//=======================================================================

void  OCCShapeAttributeSet::Read(TopoDS_Shape& S,
                                 Standard_IStream& IS,
                                 const int nbshapes,
                                 TDF_Label* label )const
{
  std::string buffer, buffer_attr;
  IS >> buffer;
  if (buffer[0] == '*')
    S = TopoDS_Shape();
  else {
    char type;
    int num;
    std::istringstream buffstr(buffer);
    buffstr >> type >> num;
    S = myShapes(nbshapes - num + 1);
    switch (type) {

    case '+' :
      S.Orientation(TopAbs_FORWARD);
      break;

    case '-' :
      S.Orientation(TopAbs_REVERSED);
      break;

    case 'i' :
      S.Orientation(TopAbs_INTERNAL);
      break;

    case 'e' :
      S.Orientation(TopAbs_EXTERNAL);
      break;
    }

    Standard_Integer l;
    IS >> l;
    S.Location(myLocations.Location(l));
  }
  if(label != NULL)
  {
    Standard_Integer i, nbShapes = myShapes.Extent();
    for ( i = 1; i <= nbShapes; i++)
    {
      TopoDS_Shape Sh = myShapes(i);
      IS >> buffer_attr;
      if(buffer_attr[0] != '*' && buffer_attr[0] != 'C')
        break;
      if(buffer_attr[0] == '*') //empty attributes for this shape
        continue;
      ReadAttribute(Sh, IS,*label);
    }
  }
}

//=======================================================================
//function : Check
//purpose  :
//=======================================================================

void OCCShapeAttributeSet::Check(const TopAbs_ShapeEnum T,
                                 TopoDS_Shape&          S)
{
  if (T == TopAbs_FACE) {
    const TopoDS_Face& F = TopoDS::Face(S);
    BRepTools::Update(F);
  }
}

//=======================================================================
//function : NbShapes
//purpose  :
//=======================================================================

int  OCCShapeAttributeSet::NbShapes() const
{
  return myShapes.Extent();
}

