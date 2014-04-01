//-------------------------------------------------------------------------
// Filename      : OCCAttribSet.cpp
//
// Purpose       : Common attribute functionality for MBG
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/01/03
//-------------------------------------------------------------------------

#include "OCCAttribSet.hpp"
//#include "OCCAttrib.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitSimpleAttrib.hpp"
#include "CubitFileIOWrapper.hpp"
#include "Handle_TDataStd_Shape.hxx"
#include "TCollection_ExtendedString.hxx"
#include "Handle_TDataStd_Name.hxx"
#include "Handle_TDataStd_ExtStringArray.hxx"
#include "Handle_TDataStd_RealArray.hxx"
#include "Handle_TDataStd_IntegerArray.hxx"
#include "TDataStd_Name.hxx"
#include "TDataStd_ExtStringArray.hxx"
#include "TDataStd_RealArray.hxx"
#include "TDataStd_IntegerArray.hxx"
#include "TDataStd_Shape.hxx"
#include "TopoDS_Shape.hxx"
#include "TDF_ChildIterator.hxx"
#include <vector>

void OCCAttribSet::FindShape(TopoDS_Shape& shape,
                             TDF_Label& aLabel,
                             CubitBoolean& found) 
{
  for (TDF_ChildIterator it1(OCCQueryEngine::instance()->mainLabel, CUBIT_FALSE);
    it1.More(); it1.Next())
  {
    //find the same shape attribute first
    aLabel = it1.Value();

    Handle_TDataStd_Shape attr_shape;
    TopoDS_Shape old_shape;
    if(aLabel.FindAttribute(TDataStd_Shape::GetID(), attr_shape))
      old_shape = attr_shape->Get(aLabel);
    if(old_shape.IsPartner(shape))
    {
      found = CUBIT_TRUE;
      break;
    }
  }
}

void OCCAttribSet::append_attribute( CubitSimpleAttrib* csa, TopoDS_Shape& shape )
{
  //Add attributes on child of myLabel
  //1. add shape attribute, first check to make sure there's no shape attribute
  CubitBoolean found = CUBIT_FALSE;
  TDF_Label aLabel;

  FindShape(shape, aLabel, found);

  if(!found)
  { 
    aLabel = OCCQueryEngine::instance()->mainLabel.NewChild();
    Handle_TDataStd_Shape attr_shape = TDataStd_Shape::Set(aLabel, shape);
    //myLabel.AddAttribute(attr_shape);
    //test if the attribute has been attached
    assert(aLabel.IsAttribute(TDataStd_Shape::GetID()));
  }

  //2. add type attribute , below attributes are added on child lable of myLabel
  //First check if the attributes are already exist, if so, create new data.
  TDF_Label lab;
  Handle_TDataStd_Name attr_name;
  TCollection_ExtendedString name_string;
  CubitString* type = csa->string_data_list()->get_and_step();
  TCollection_ExtendedString cstring( (Standard_CString)type->c_str() );
  found = CUBIT_FALSE;
  for (TDF_ChildIterator it(aLabel,CUBIT_FALSE); it.More(); it.Next()) 
  {
    lab = it.Value();
    found = find_attribute(lab, csa);
    if(!found)
    {
      //check if the type attribute is the same, disable it
      TCollection_ExtendedString old_string;
      if(lab.FindAttribute(TDataStd_Name::GetID(), attr_name))
      {
        old_string = attr_name->Get();

        if(old_string == cstring)
          lab.ForgetAllAttributes();
      }
    }
    if(found)
      break;
  }

  //if not exist,  create new child label
  if (!found)
  {
    lab = aLabel.NewChild();
    attr_name = TDataStd_Name::Set(lab, cstring);
    assert(lab.IsAttribute(TDataStd_Name::GetID()));
  }
     
  else//nothing to do
    return;

  //3. add string attribute
  DLIList<CubitString*>* strings = csa->string_data_list();
  Handle_TDataStd_ExtStringArray attr_string;
  int size = strings->size()-1;
  if(strings && size > 0)
  {
    //set the length of String Array .
    attr_string = TDataStd_ExtStringArray::Set(lab, 0, size-1);
  }
  strings->reset();
  strings->step(); //type string filters out
  for(int i = 0; strings && i < size; i++)
  {
    char const* string1 = strings->get_and_step()->c_str();
    TCollection_ExtendedString 
       cstring(string1, Standard_True );
    attr_string->SetValue(i, cstring) ;
  }
  
  if (strings && size > 0)
    assert(lab.IsAttribute(TDataStd_ExtStringArray::GetID()));

  //4. add double attribute
  DLIList<double*>* doubles = csa->double_data_list();
  Handle_TDataStd_RealArray attr_double;
  size = doubles->size();
  if(doubles && size > 0)
  {
    //set the length of double array .
    attr_double = TDataStd_RealArray::Set(lab,0, size-1);
  } 
  doubles->reset();
  for(int i = 0; doubles && i < size; i++)
    attr_double->SetValue(i, *(doubles->get_and_step()));

  if(doubles && size > 0)
    assert(lab.IsAttribute(TDataStd_RealArray::GetID()));
    
  //5. add int attribute
  DLIList<int*>* ints = csa->int_data_list();
  Handle_TDataStd_IntegerArray attr_int;
  size = ints->size();
  if(ints && size > 0)
  {
    //set the length of int array be 11, can be extended.
    attr_int = TDataStd_IntegerArray::Set(lab, 0, size-1);
  }
  ints->reset();
  for(int i = 0; ints && i < size; i++)
    attr_int->SetValue(i, *(ints->get_and_step()));

  if(ints && size > 0)
    assert(lab.IsAttribute(TDataStd_IntegerArray::GetID()));
}

CubitBoolean OCCAttribSet::find_attribute(TDF_Label child,
                                          CubitSimpleAttrib* csa)
{
  DLIList<int*>* ints = csa->int_data_list();
  CubitString type = csa->character_type();
  char const* string1 = type.c_str();
  Standard_Boolean isMultiByte = Standard_True;
  TCollection_ExtendedString cstring( string1, isMultiByte );
  DLIList<double*>* doubles = csa->double_data_list();
  DLIList<CubitString*>* strings = csa->string_data_list();

  //find the same type attribute first
  Handle_TDataStd_Name attr_name;
  TCollection_ExtendedString old_string;
  if(child.FindAttribute(TDataStd_Name::GetID(), attr_name))
  {
    old_string = attr_name->Get();

    if(old_string != cstring)
      return CUBIT_FALSE;
  }

  else
    return CUBIT_FALSE;

  //continue to compare the rest attributes.
  CubitBoolean is_same = CUBIT_TRUE;
  if(ints->size() > 0 )
  {
    Handle_TDataStd_IntegerArray attr_ints;
    if(child.FindAttribute(TDataStd_IntegerArray::GetID(), attr_ints) &&
       attr_ints->Length() == ints->size())
    {
      for(int i = 0; i < ints->size(); i++)
        if(attr_ints->Value(i) != *ints->get_and_step())
        {
          is_same = CUBIT_FALSE;
          break;
        }
    }
  }
  if(!is_same)
    return is_same;

  if(doubles->size() > 0 )
  {
    Handle_TDataStd_RealArray attr_doubles;
    if(child.FindAttribute(TDataStd_RealArray::GetID(), attr_doubles) &&
       attr_doubles->Length() == doubles->size())
    {
      for(int i = 0; i < doubles->size(); i++)
        if(attr_doubles->Value(i) != *doubles->get_and_step())
        {
          is_same = CUBIT_FALSE;
          break;
        }
    }
  }

  if(!is_same)
    return is_same;

  if(strings->size() > 0 )
  {
    Handle_TDataStd_ExtStringArray attr_strings;
    if(child.FindAttribute(TDataStd_ExtStringArray::GetID(), attr_strings) &&
       attr_strings->Length() == strings->size()-1)
    {
      strings->step();//filter out type string
      for(int i = 0; i < strings->size()-1; i++)
      {
        CubitString astring = *strings->get_and_step();
        char const* string1 = astring.c_str();
        TCollection_ExtendedString string( string1 , Standard_True);
        if(attr_strings->Value(i) != string)
        {
          is_same = CUBIT_FALSE;
          break;
        }
      }
    }
  }
  return is_same  ;
}

void OCCAttribSet::remove_attribute( CubitSimpleAttrib* csa)
{
  CubitString type = csa->character_type();
  if(type.length() == 0)
      return;

  TDF_Label myLabel;
  for (TDF_ChildIterator it(OCCQueryEngine::instance()->mainLabel, CUBIT_FALSE); it.More(); it.Next())
  {
    //find the same type attribute first
    myLabel = it.Value();

    //forget csa attribute from the document
    for (TDF_ChildIterator it1(myLabel, CUBIT_FALSE); it1.More(); it1.Next())
    {
      //find the same type attribute first
      TDF_Label child = it1.Value();

      CubitBoolean found = find_attribute(child, csa);

      if(found)
        child.ForgetAllAttributes( );
    }
  }
}

void OCCAttribSet::remove_attribute( TopoDS_Shape& shape)
{
  CubitBoolean found = CUBIT_FALSE;
  TDF_Label myLabel;

  FindShape(shape, myLabel, found);

  if(!found)
    return;

  myLabel.ForgetAllAttributes();
}

void OCCAttribSet::remove_attribute( CubitSimpleAttrib* csa, 
                                     TopoDS_Shape& shape)
{
  CubitBoolean found = CUBIT_FALSE;
  TDF_Label myLabel;

  FindShape(shape, myLabel, found);

  if(!found)
    return;

  for (TDF_ChildIterator it1(myLabel, CUBIT_FALSE); it1.More(); it1.Next())
  {
    //find the same type attribute first
    TDF_Label child = it1.Value();

    if(!csa)
      child.ForgetAllAttributes( );

    else
    {
      CubitBoolean found = find_attribute(child, csa);

      if(found)
        child.ForgetAllAttributes( );
    }
  }
}

void OCCAttribSet::get_attributes(TDF_Label &lab,
                                  DLIList<CubitSimpleAttrib*>& list)
{
  DLIList<CubitString*> strings;
  CubitString* string;
  DLIList<double> doubles;
  DLIList<int> ints;

  Handle_TDataStd_Name attr_name;
  TCollection_ExtendedString name_string;
  Handle_TDataStd_ExtStringArray attr_string;
  Handle_TDataStd_RealArray attr_doubles;
  Handle_TDataStd_IntegerArray attr_ints;

  if(lab.FindAttribute(TDataStd_Name::GetID(), attr_name))
  {
    name_string = attr_name->Get();
    int length = name_string.Length();
    std::vector<char> temp_type(length+1);
    for (int i = 1 ;  i <= length; i++)
    {
      Standard_ExtCharacter c = name_string.Value(i);
      temp_type[i-1] = ToCharacter(c);
    }
    temp_type[length] = '\0';
    string = new CubitString(&temp_type[0]);
    strings.append(string);
  }
  if(lab.FindAttribute(TDataStd_ExtStringArray::GetID(), attr_string))
  {
    int length = attr_string->Length();
    for (int i = 0 ;  i < length; i++)
    {
      name_string = attr_string->Value(i);
      int length2 = name_string.Length();
      std::vector<char> temp_string(length2+1);
      for(int j = 1; j <= length2; j ++)
      {
        Standard_ExtCharacter c = name_string.Value(j);
        temp_string[j-1] = ToCharacter(c);
      }
      temp_string[length2] = '\0';
      const char *s = &temp_string[0];
      string = new CubitString(s);
      strings.append(string);
    }
  }
  if(lab.FindAttribute(TDataStd_RealArray::GetID(), attr_doubles))
  {
    for(int i = 0; i < attr_doubles->Length(); i++)
    {
      double d = attr_doubles->Value(i);
      doubles.append(d);
    }
  }
  if(lab.FindAttribute(TDataStd_IntegerArray::GetID(), attr_ints))
  {
    for(int i = 0; i < attr_ints->Length(); i++)
    {
      int j = attr_ints->Value(i);
      ints.append(j);
    }
  }
  CubitSimpleAttrib *tmp_attrib;
  if(strings.size() > 0 || doubles.size() || ints.size())
  {
    tmp_attrib = new CubitSimpleAttrib(&strings, &doubles, &ints);
    list.append(tmp_attrib);
  }
  for(int i = 0; i < strings.size(); i++)
    delete strings.get_and_step();
}

CubitStatus OCCAttribSet::get_attributes( TopoDS_Shape& shape,
                                          DLIList<CubitSimpleAttrib*>& list) 
{
  CubitBoolean found = CUBIT_FALSE;
  TDF_Label myLabel;

  FindShape(shape, myLabel, found);

  if(!found)
    return CUBIT_FAILURE;

  DLIList<CubitString*> strings;
  CubitString* string;
  DLIList<double> doubles;
  DLIList<int> ints;
  TDF_Label lab;
  Handle_TDataStd_Name attr_name;
  TCollection_ExtendedString name_string;
  Handle_TDataStd_ExtStringArray attr_string;
  Handle_TDataStd_RealArray attr_doubles;
  Handle_TDataStd_IntegerArray attr_ints;
  for (TDF_ChildIterator it(myLabel,CUBIT_FALSE); it.More(); it.Next())
  {
    lab = it.Value();
    get_attributes(lab, list);
  }     
  return CUBIT_SUCCESS;
}

CubitStatus OCCAttribSet::get_attributes( const CubitString& name,
                                          TopoDS_Shape& shape,
                                    DLIList<CubitSimpleAttrib*>& list ) 
{
  CubitBoolean found = CUBIT_FALSE;
  TDF_Label aLabel;

  FindShape(shape, aLabel, found);

  if(!found)
    return CUBIT_FAILURE;

  Handle_TDataStd_Name attr_name;
  TCollection_ExtendedString cstring( (Standard_CString)name.c_str() );
  TDF_Label lab;
  for (TDF_ChildIterator it(aLabel,CUBIT_FALSE); it.More(); it.Next())
  {
    lab = it.Value();

    //check if the type attribute is the same, disable it
    TCollection_ExtendedString old_string;
    if(lab.FindAttribute(TDataStd_Name::GetID(), attr_name))
    {
      old_string = attr_name->Get();

      if(old_string == cstring)
      {
        get_attributes(lab, list);
        break;
      }
    }
  }
 
  return CUBIT_SUCCESS;
}


int OCCAttribSet::attribute_count() 
{
  int count = 0;
  
  return count;
}

