#ifdef WIN32
#pragma warning(disable : 4786)
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <iostream>

#include "SettingHandler.hpp"
#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitUtil.hpp"
#include "CubitMessage.hpp"
#include "AppUtil.hpp"


SettingHandler* SettingHandler::instance_ = NULL;

SettingHandler* SettingHandler::instance()
{	
  if (!instance_)
    {	
      instance_ = new SettingHandler;
      if (!instance_)
	{
	  std::cerr << " *** Unable to instantiate setting_handler object ***" << std::endl;
	  exit(1);
	}
    }
  return instance_;
}

//Default Constructor
SettingHandler::SettingHandler() 
{
  debug_flags_added = CUBIT_FALSE;
} 

void SettingHandler::delete_instance() 
{
  if (instance_)
  {
    delete instance_;
    instance_ = 0;
  }
}

SettingHandler::~SettingHandler() 
{
    std::map<CubitString, SettingHolder*>::iterator iter = mSettingsList.begin();
  std::map<CubitString, SettingHolder*>::iterator end = mSettingsList.end();

  for( ; iter != end; ++iter)
    delete (*iter).second;
  mSettingsList.clear();
}

void SettingHandler::add_setting(char* name, void (*setFn) (int), 
				int (*getFn)())
{
  CubitString cs = name;
  //Check to see if a setting with the same name already exists
  if (mSettingsList.find(cs) != mSettingsList.end()) {
    std::cerr << "The " << name << " setting has the same name as another setting.\n" 
		<< "This is a bug.  Please report!\n";
    exit(-1);
  }
  SettingHolder* setting = new SettingHolder(cs, setFn, getFn);
  mSettingsList[cs] = setting;
}


void SettingHandler::add_setting(char* name, void (*setFn) (double), 
				double (*getFn)())
{
  CubitString cs = name;
  //Check to see if a setting with the same name already exists
  if (mSettingsList.find(cs) != mSettingsList.end()) {
    std::cerr << "The " << name << " setting has the same name as another setting.\n" 
		<< "This is a bug.  Please report!\n";
    exit(-1);
    
  }
  SettingHolder* setting = new SettingHolder(cs, setFn, getFn);
  mSettingsList[cs] = setting;
}


void SettingHandler::add_setting(char* name, void (*setFn) (CubitBoolean), 
				CubitBoolean (*getFn)())
{
   CubitString cs = name;
  //Check to see if a setting with the same name already exists
  if (mSettingsList.find(cs) != mSettingsList.end()) {
    std::cerr << "The " << name << " setting has the same name as another setting.\n" 
		<< "This is a bug.  Please report!\n";
    exit(-1);
    
  }
  SettingHolder* setting = new SettingHolder(cs, setFn, getFn);
  mSettingsList[cs] = setting;
}


void SettingHandler::add_setting(char* name, void (*setFn) (CubitString), 
				CubitString (*getFn)())
{
   CubitString cs = name;
  //Check to see if a setting with the same name already exists
  if (mSettingsList.find(cs) != mSettingsList.end()) {
    std::cerr << "The " << name << " setting has the same name as another setting.\n" 
		<< "This is a bug.  Please report!\n";
    exit(-1);
    
  }
  SettingHolder* setting = new SettingHolder(cs, setFn, getFn);
  mSettingsList[cs] = setting;
}

void SettingHandler::save_settings(const char *filename)
{
  
  FILE* file = fopen(filename, "w");
  
  if (file == NULL) {
    std::cerr << "File " << filename << " could not be opened.  Settings not saved." << std::endl;
    return;
   }

  //Put a header at the top of the file
  fprintf(file, "#Setting Name                                 Setting Value\n");
  fprintf(file, "#------------------------------------------------------------\n");

  // Get starting and ending iterators from the map
  std::map<CubitString, SettingHolder*>::iterator start = mSettingsList.begin();
  std::map<CubitString, SettingHolder*>::iterator end = mSettingsList.end();
  
  // Step through the map and build a file to hold the settings
  while(start != end) {
    SettingHolder* setting = (*start).second;

    //Integer settings
    if (setting->setting_type == 0) {
      fprintf(file, "%-50.50s%d \n", setting->name.c_str(), (setting->get_int_function)());      
    }
    
    //Double settings
    else if (setting->setting_type == 1) {
      fprintf(file, "%-50.50s%f \n", setting->name.c_str(), (setting->get_double_function)());
    }

    //Boolean settings
    else if (setting->setting_type == 2) {
      fprintf(file, "%-50.50s%d \n", setting->name.c_str(), (setting->get_bool_function)());
    }
    
    //String settings
    else if (setting->setting_type == 3) {
      fprintf(file, "%-50.50s%s \n", setting->name.c_str(), (setting->get_string_function)().c_str());
    }

    else
    {
      assert(false);
      PRINT_ERROR("Error with SettingHolder type!  Please report.");
    }

    ++start;
  }
  
  fclose(file);
 
}

void SettingHandler::print_settings()
{
  int rows, cols;
  if (!AppUtil::instance()->get_terminal_size(rows, cols))
    cols = 71;
  --cols;

  // Get starting and ending iterators from the map
  std::map<CubitString, SettingHolder*>::iterator iter = mSettingsList.begin();
  const std::map<CubitString, SettingHolder*>::iterator end = mSettingsList.end();
  
  // Step through the map and build a file to hold the settings
  for ( ; iter != end; ++iter ) 
  {
    SettingHolder* setting = iter->second;
    int vallen = cols - setting->name.length() - 1;
    if (vallen < 1) vallen = 1;
    switch( setting->setting_type )
    {
      case 0:
        PRINT_INFO("%s %*d\n", setting->name.c_str(), vallen,
            (setting->get_int_function)());
        break;
      case 1:
        PRINT_INFO("%s %*f\n", setting->name.c_str(), vallen,
            (setting->get_double_function)());
        break;
      case 2:
        PRINT_INFO("%s %*s\n", setting->name.c_str(), vallen, 
            (setting->get_bool_function)() ? "true" : "false");
        break;
      case 3:
        PRINT_INFO("%s %*s\n", setting->name.c_str(), vallen, 
            (setting->get_string_function)().c_str());
        break;
      default:
        PRINT_INFO("%s %*s\n", setting->name.c_str(), vallen,
            "ERROR : UNKNOWN SETTING TYPE");
    }
  }
}

void SettingHandler::save_settings()
{
  
  char* default_filename = "cubit.settings";

  FILE* file = fopen(default_filename, "w");
  
  if (file == NULL) {
    std::cerr << "File " << default_filename << " could not be opened.  Settings not saved." << std::endl;
    return;
   }

  //Put a header at the top of the file
  fprintf(file, "#Setting Name                                 Setting Value\n");
  fprintf(file, "#------------------------------------------------------------\n");

  // Get starting and ending iterators from the map
  std::map<CubitString, SettingHolder*>::iterator start = mSettingsList.begin();
  std::map<CubitString, SettingHolder*>::iterator end = mSettingsList.end();

  // Step through the map and build a file to hold the settings
  while(start != end) {
    SettingHolder* setting = (*start).second;

    //Integer settings
    if (setting->setting_type == 0) {
      fprintf(file, "%-50.50s%d \n", setting->name.c_str(), (setting->get_int_function)());      
    }
    
    //Double settings
    else if (setting->setting_type == 1) {
      fprintf(file, "%-50.50s%f \n", setting->name.c_str(), (setting->get_double_function)());
    }

    //Boolean settings
    else if (setting->setting_type == 2) {
      fprintf(file, "%-50.50s%d \n", setting->name.c_str(), (setting->get_bool_function)());
    }
    
    //String settings
    else if (setting->setting_type == 3) {
      fprintf(file, "%-50.50s%s \n", setting->name.c_str(), (setting->get_string_function)().c_str());
    }

    else
    {
      assert(false);
      PRINT_ERROR("Error with SettingHolder type!  Please report.");
    }

    ++start;
  }
  
  fclose(file);

}


void SettingHandler::restore_settings(const char* filename)
{

  //Open the file for reading
  FILE* file = fopen(filename, "r");

  if (file == NULL) {
    std::cerr << "File " << filename << " could not be opened.  Settings not restored." << std::endl;
    return;
   }

  //Read the first 2 lines of the file, we know they are not settings
  //so just get rid of them
  char junk[100];
  fgets(junk, 100, file);
  fgets(junk, 100, file);

  char name[51]; //Allocate 50 bytes for the characters and 1 for a null termination
  char new_name[50]; //A string that will hold the setting name without excess white space
  CubitString cubit_name; //A CubitString that will be used for checking if the setting
                          //is in the map.
  char value[120]; //This will be used when getting the value of a setting
  SettingHolder* setting;
  int int_value = -1;
  double double_value = 0.0;
  CubitBoolean bool_value = CUBIT_FALSE;
  CubitString string_value = "";

  while (!feof(file)) {

  //Get the setting name.  This will be in the first 50 characters of a line
  fgets(name, 51, file);

  //Create a CubitString for the setting name without all the white space
  int i;
  for (i = (strlen(name)-1); i >= 0; i--) {
    if (name[i] != 32 && name[i] != 9 && name[i] != 11) break;
  }

  //Get the first i+1 characters from name
  strncpy(new_name, name, i+1);
  new_name[i+1] = 0; //Add the null character

  //Create a CubitString
  cubit_name = new_name;

  //Find the setting in the map
  std::map<CubitString, SettingHolder*>::iterator temp;

  temp = mSettingsList.find(cubit_name);
  
    //Read in the rest of the line no matter what
  fgets(value, 120, file);
    
  if (temp == mSettingsList.end()) {
    std::cerr << "Setting " << cubit_name.c_str() << " was not found." << std::endl;  
  }
  else {
    setting = (*temp).second;
    //Integer Settings
    if (setting->setting_type == 0) {
      int_value = atoi(value);
      (setting->set_int_function)(int_value);
    }
    
    //Double Settings
    if (setting->setting_type == 1) {
      double_value = atof(value);
      (setting->set_double_function)(double_value);
    }    

    //Boolean Settings
    if (setting->setting_type == 2) {
//      bool_value = (CubitBoolean)(atoi(value));
      bool_value = (atoi(value)) ? true : false;
      (setting->set_bool_function)(bool_value);
    }    
    
    //CubitString Settings
    if (setting->setting_type == 3) {
      string_value = value;
      (setting->set_string_function)(string_value);
    }     

  }

  } //End while

  fclose(file);

}

SettingHolder* SettingHandler::get_setting_holder(CubitString name)
{ 
	if( mSettingsList.find(name) == mSettingsList.end() )
		return NULL;

	return (*(mSettingsList.find(name))).second;
}

void SettingHandler::get_settings_list(std::vector< std::pair<CubitString, SettingHolder*> > &list)
{
  // return a list of the settings
  std::map<CubitString, SettingHolder*>::iterator iter = mSettingsList.begin();
  while(iter != mSettingsList.end())
  {
    std::pair<CubitString, SettingHolder*>  tmp_pair;
    CubitString key = (*iter).first;
    tmp_pair.first = key;
    tmp_pair.second = (*iter).second;
    list.push_back( (tmp_pair) );
    iter++;
  }
}







