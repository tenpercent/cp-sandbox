//- Filename:    CubitSimpleAttrib.cc
//- Description: Implementation of class CubitSimpleAttrib
//- Owner:       Greg Nielson
//- Checked by:
//- Version:

#include <string.h>
#include "SettingHandler.hpp"

#include "CubitSimpleAttrib.hpp"
#include "CubitString.hpp"
#include "RefEntity.hpp"
#include "RefEntityFactory.hpp"
#include "CubitFileIOWrapper.hpp"

CubitBoolean CubitSimpleAttrib::pushAttribs = CUBIT_FALSE;

//SettingHolder SettingHandler::settingObject43 = SettingHolder("Push_Attribs", CubitSimpleAttrib::set_push_attribs, CubitSimpleAttrib::get_push_attribs, CUBIT_FALSE);

CubitSimpleAttrib::CubitSimpleAttrib()
{
  characterType = "";
}

CubitSimpleAttrib::CubitSimpleAttrib( CubitString new_type )
{
  characterType = new_type;
}

CubitSimpleAttrib::~CubitSimpleAttrib()
{
  int i;
  for ( i = stringDataList.size(); i > 0; i--)
    delete stringDataList.get_and_step();
  
  for (i = doubleDataList.size(); i > 0; i--)
    delete doubleDataList.get_and_step();
  
  for (i = intDataList.size(); i > 0; i--)
    delete intDataList.get_and_step();
  
  stringDataList.clean_out();
  intDataList.clean_out();
  doubleDataList.clean_out();
}

CubitSimpleAttrib::CubitSimpleAttrib(const CubitString new_character_type,
                                     const CubitString new_string_data,
                                     const CubitString new_more_string_data,
                                     const int new_integer_data,
                                     const double new_double_data)
{
  characterType = new_character_type;

  stringDataList.clean_out();
  doubleDataList.clean_out();
  intDataList.clean_out();

  CubitString* cstemp = new CubitString(characterType);
  stringDataList.append_unique(cstemp);

  if (new_string_data != "") {
    cstemp = new CubitString(new_string_data);
    stringDataList.append_unique(cstemp);
  }
  
  if (new_more_string_data != "") {
    cstemp = new CubitString(new_more_string_data);
    stringDataList.append_unique(cstemp);
  }
  
  doubleDataList.append( new double(new_double_data) );
  intDataList.append( new int(new_integer_data) );
  characterType = "NEW_SIMPLE_ATTRIB"; 

  return;
}

CubitSimpleAttrib::CubitSimpleAttrib(CubitSimpleAttrib *csa_ptr)
{
  initialize_from_lists_of_ptrs(csa_ptr->string_data_list(),
                                csa_ptr->double_data_list(),
                                csa_ptr->int_data_list());
}

CubitSimpleAttrib::CubitSimpleAttrib(DLIList<CubitString*> *string_list,
                                     DLIList<double> *double_list,
                                     DLIList<int> *int_list)
{
  initialize_from_lists(string_list, double_list, int_list);
}


CubitSimpleAttrib &CubitSimpleAttrib::operator=(CubitSimpleAttrib &other)
{
    if (this == &other)
        return *this;

    // clean out previous data
    int i;
    for ( i = stringDataList.size(); i > 0; i--)
        delete stringDataList.get_and_step();

    for (i = doubleDataList.size(); i > 0; i--)
        delete doubleDataList.get_and_step();

    for (i = intDataList.size(); i > 0; i--)
        delete intDataList.get_and_step();

    stringDataList.clean_out();
    intDataList.clean_out();
    doubleDataList.clean_out();

    // copy the data
    initialize_from_lists_of_ptrs(other.string_data_list(),
                                  other.double_data_list(),
                                  other.int_data_list());

    return *this;
}

void CubitSimpleAttrib::initialize_from_lists_of_ptrs(DLIList<CubitString*>* string_list,
                                                      DLIList<double*>* double_list,
                                                      DLIList<int*>* int_list)
{
  if(string_list != NULL) {
    string_list->reset();
    for (int i = string_list->size(); i > 0; i--)
      stringDataList.append(new CubitString(*string_list->get_and_step()));
  }
  
  if(double_list != NULL) {
    double_list->reset();
    for (int i = double_list->size(); i > 0; i--)
      doubleDataList.append(new double(*double_list->get_and_step()));
  }
  
  if(int_list != NULL) {
    int_list->reset();
    for (int i = int_list->size(); i > 0; i--)
      intDataList.append(new int(*int_list->get_and_step()));
  }
  
  characterType = "NEW_SIMPLE_ATTRIB";
}

void CubitSimpleAttrib::initialize_from_lists(DLIList<CubitString*> *string_list,
                                              DLIList<double> *double_list,
                                              DLIList<int> *int_list)
{
  int i;
  if( string_list != NULL ) 
  {
    string_list->reset();
    for (i=string_list->size(); i--; ) 
      stringDataList.append( new CubitString( *string_list->get_and_step() ) );
  }

  if( double_list != NULL ) 
  {
    double_list->reset();
    for (i=double_list->size(); i--; ) 
      doubleDataList.append( new double( double_list->get_and_step() ) );
  }

  if( int_list != NULL ) 
  {
    int_list->reset();
    for (i=int_list->size(); i--; ) 
      intDataList.append( new int( int_list->get_and_step() ) );
  }
  
  characterType = "NEW_SIMPLE_ATTRIB";
}




CubitSimpleAttrib::CubitSimpleAttrib(const char* new_character_type,
                                     const char* new_string_data,
                                     const char* new_more_string_data,
                                     const int *new_integer_data,
                                     const double *new_double_data)
{
  if ( new_character_type != NULL )
    characterType = CubitString (new_character_type);
  else
    characterType = "UNDEFINED";

  CubitString* cstemp = new CubitString(characterType);
  stringDataList.append_unique(cstemp);

  if (new_string_data) {
    cstemp = new CubitString(new_string_data);
    stringDataList.append_unique(cstemp);
  }
  if (new_more_string_data) {
    cstemp = new CubitString(new_more_string_data);
    stringDataList.append_unique(cstemp);
  }
  
  if (new_double_data) 
    doubleDataList.append( new double(*new_double_data) );

  if (new_integer_data)
    intDataList.append( new int(*new_integer_data) );
  
  characterType = "NEW_SIMPLE_ATTRIB"; 
  return;
}

void CubitSimpleAttrib::character_type(CubitString new_type)
{
  if(characterType != "NEW_SIMPLE_ATTRIB")
  {
    characterType = new_type;
  }
  else
  {
    stringDataList.reset();
    if(stringDataList.size() >= 1)
    {
      CubitString *temp = stringDataList.get();
      *temp = new_type;
    }
    else
    {
      CubitString *temp = new CubitString(new_type);
      stringDataList.append_unique(temp);
    }
  }
}

CubitString CubitSimpleAttrib::character_type()//- returns the type of attribute
{
  if(characterType != "NEW_SIMPLE_ATTRIB")
  {
    return characterType;
  }
  else
  {
    stringDataList.reset();
    if(stringDataList.size() >= 1 && stringDataList.get() != NULL )
    {
      return *(stringDataList.get());
    }
    else
    {
      return CubitString("");
    }
  } 
}

void CubitSimpleAttrib::string_data_list( DLIList<CubitString*>* new_list)
{
  int i;
  //clean out, deleting each one
  CubitString *tmp_string;
  for(i=stringDataList.size(); i--;)
  {
    tmp_string = stringDataList.get_and_step();
    delete tmp_string; 
  } 
  stringDataList.clean_out();

  //new each one
  for( i=new_list->size(); i--; )
  {
    tmp_string = new CubitString( *(new_list->get_and_step()) ); 
    stringDataList.append( tmp_string );
  }
} 

void CubitSimpleAttrib::double_data_list(const DLIList<double*>* new_list)
{
  doubleDataList = *new_list;
}

void CubitSimpleAttrib::int_data_list(const DLIList<int*>* new_list)
{
  intDataList = *new_list;
}

CubitBoolean CubitSimpleAttrib::equivalent(CubitSimpleAttrib *first_attrib_ptr,
                                           CubitSimpleAttrib *second_attrib_ptr)
{
  CubitBoolean return_val;
  if(first_attrib_ptr == NULL && second_attrib_ptr != NULL)
    return_val = CUBIT_FALSE;
  else if (first_attrib_ptr != NULL && second_attrib_ptr == NULL)
    return_val = CUBIT_FALSE;
  else if (first_attrib_ptr == NULL && second_attrib_ptr == NULL)
    return_val = CUBIT_TRUE;
  else {
      // compare string lists
    DLIList<CubitString*> *first =
        first_attrib_ptr->string_data_list();
    DLIList<CubitString*> *second =
        second_attrib_ptr->string_data_list();
    if (first->size() != second->size()) return CUBIT_FALSE;
    first->reset();
    second->reset();
    int i;
    for (i = first->size(); i > 0; i--)
      if (*first->get_and_step() != *second->get_and_step()) return CUBIT_FALSE;

      // now compare int and double lists
    if ((first_attrib_ptr->int_data_list()->size() !=
        second_attrib_ptr->int_data_list()->size()) ||
        (first_attrib_ptr->double_data_list()->size() !=
        second_attrib_ptr->double_data_list()->size()))
      return CUBIT_FALSE;

    first_attrib_ptr->int_data_list()->reset();
    second_attrib_ptr->int_data_list()->reset();
    for (i = first_attrib_ptr->int_data_list()->size(); i > 0; i--)
      if (*first_attrib_ptr->int_data_list()->get_and_step() !=
          *second_attrib_ptr->int_data_list()->get_and_step()) return CUBIT_FALSE;
    
    first_attrib_ptr->double_data_list()->reset();
    second_attrib_ptr->double_data_list()->reset();
    for (i = first_attrib_ptr->double_data_list()->size(); i > 0; i--)
      if (*first_attrib_ptr->double_data_list()->get_and_step() !=
          *second_attrib_ptr->double_data_list()->get_and_step())
        return CUBIT_FALSE;

    return_val = CUBIT_TRUE;
  }
  return return_val;
}


//Initialize all settings in this class
void CubitSimpleAttrib::initialize_settings()
{

  SettingHandler::instance()->add_setting("Push Attribs", 
                                         CubitSimpleAttrib::set_push_attribs, 
					 CubitSimpleAttrib::get_push_attribs); 
}

void CubitSimpleAttrib::print()
{
  stringDataList.reset();
  intDataList.reset();
  doubleDataList.reset();
  
  PRINT_INFO("CSA: type = %s\n", stringDataList.get_and_step()->c_str());

  PRINT_INFO("String data: ");
  int i;
  for ( i = stringDataList.size()-1; i > 0; i--)
      PRINT_INFO("%s;", stringDataList.get_and_step()->c_str());
  if (stringDataList.size() == 0) PRINT_INFO("(none)");
  
  PRINT_INFO("\n");
  
  PRINT_INFO("Int data: ");
  for (i = intDataList.size(); i > 0; i--)
      PRINT_INFO("%d;", *intDataList.get_and_step());
  if (intDataList.size() == 0) PRINT_INFO("(none)");
  PRINT_INFO("\n");
  
  PRINT_INFO("Double data: ");
  for (i = doubleDataList.size(); i > 0; i--)
      PRINT_INFO("%f;", *doubleDataList.get_and_step());
  if (doubleDataList.size() == 0) PRINT_INFO("(none)");
  PRINT_INFO("\n");
}


void CubitSimpleAttrib::set( const CubitString& type,
                             const DLIList<int>* int_list,
                             const DLIList<double>* real_list,
                             const DLIList<CubitString*>* string_list )
{
  character_type( type );
  set( int_list ? *int_list : DLIList<int>(0) );
  set( real_list ? *real_list : DLIList<double>(0) );
  set( string_list ? *string_list : DLIList<CubitString*>(0) );
}

void CubitSimpleAttrib::set( const DLIList<int>& data )
{
  intDataList.reset();
  
  int i, count = CUBIT_MIN(intDataList.size(),data.size());
  for( i = 0; i < count; i++ )
    *intDataList.get_and_step() = data.next(i);
  for( ; i < data.size(); i++ )
    intDataList.append( new int(data.next(i)) );
  while( intDataList.size() > data.size() )
    delete intDataList.pop();
}

void CubitSimpleAttrib::set( const DLIList<double>& data )
{
  doubleDataList.reset();
  
  int i, count = CUBIT_MIN(doubleDataList.size(),data.size());
  for( i = 0; i < count; i++ )
    *doubleDataList.get_and_step() = data.next(i);
  for( ; i < data.size(); i++ )
    doubleDataList.append( new double(data.next(i)) );
  while( doubleDataList.size() > data.size() )
    delete doubleDataList.pop();
}

void CubitSimpleAttrib::set( const DLIList<CubitString*>& data )
{
  stringDataList.reset();
  stringDataList.step(); // first entry is type.
  
  int i, count = CUBIT_MIN(stringDataList.size(),data.size()) - 1;
  for( i = 0; i < count; i++ )
    *stringDataList.get_and_step() = *data.next(i);
  for( ; i < data.size(); i++ )
    stringDataList.append( new CubitString(*data.next(i)) );
  while( stringDataList.size() > (data.size()+1) )
    delete stringDataList.pop();
}

