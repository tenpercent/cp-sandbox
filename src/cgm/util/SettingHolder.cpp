#include "SettingHolder.hpp"
#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitUtil.hpp"
#include "CubitMessage.hpp"


/////////////////////SettingHolder CONSTRUCTORS/////////////////////////


SettingHolder::SettingHolder()
{
 setting_type = SETTING_UNKNOWN_VALUE_TYPE;
  description = "";
  name = "";

  get_int_function    = NULL;
  set_int_function    = NULL; 
  get_double_function = NULL;
  set_double_function = NULL;
  get_bool_function   = NULL;
  set_bool_function   = NULL;
  get_string_function = NULL;
  set_string_function = NULL;

}

// constructor for int settings
SettingHolder::SettingHolder(CubitString n, void (*setFn) (int), 
			     int (*getFn)(), int value)
{
  setting_type = SETTING_INT_VALUE_TYPE;
  int_value    = value;
  double_value = 0.0;
  string_value = "";
  description = "";
  name = n; 

  get_int_function    = getFn;
  set_int_function    = setFn; 
  get_double_function = NULL;
  set_double_function = NULL;
  get_bool_function   = NULL;
  set_bool_function   = NULL;
  get_string_function = NULL;
  set_string_function = NULL;

}

// constructor for double settings
SettingHolder::SettingHolder(CubitString n, void (*setFn) (double), 
			     double (*getFn)(), double value)
{ 
  setting_type = SETTING_DOUBLE_VALUE_TYPE;
  int_value    = 0;
  double_value = value;
  string_value = "";
  description = "";
  name = n;

  get_int_function    = NULL;
  set_int_function    = NULL;
  get_double_function = getFn;
  set_double_function = setFn;
  get_bool_function   = NULL;
  set_bool_function   = NULL;
  get_string_function = NULL;
  set_string_function = NULL;

}


// constructor for boolean settings
SettingHolder::SettingHolder(CubitString n, void (*setFn) (CubitBoolean), 
			     CubitBoolean (*getFn)(), CubitBoolean)
{ 
  setting_type = SETTING_CUBITBOOLEAN_VALUE_TYPE;
  int_value    = 0;
  double_value = 0.0;
  description = "";
  name = n;

  get_int_function    = NULL;
  set_int_function    = NULL;
  get_double_function = NULL;
  set_double_function = NULL;
  get_bool_function   = getFn;
  set_bool_function   = setFn;
  get_string_function = NULL;
  set_string_function = NULL;
  
}

// constructor for string settings
SettingHolder::SettingHolder(CubitString n, void (*setFn) (CubitString), 
			     CubitString (*getFn)(), char* value)
{ 
  setting_type = SETTING_CUBITSTRING_VALUE_TYPE;
  int_value    = 0;
  double_value = 0.0;
  string_value = value;
  description = "";
  name = n;
  
  get_int_function   = NULL;
  set_int_function   = NULL;
  get_double_function = NULL;
  set_double_function = NULL;
  get_bool_function   = NULL;
  set_bool_function   = NULL;
  get_string_function = getFn;
  set_string_function = setFn;  

}

// This is added for settings function.
// constructor for int settings
SettingHolder::SettingHolder(CubitString n, void (*setFn) (int),
			     int (*getFn)())
{
  setting_type = SETTING_INT_VALUE_TYPE;
  description = "";
  name = n;

  get_int_function    = getFn;
  set_int_function    = setFn;
  get_double_function = NULL;
  set_double_function = NULL;
  get_bool_function   = NULL;
  set_bool_function   = NULL;
  get_string_function = NULL;
  set_string_function = NULL;

}

// constructor for double settings
SettingHolder::SettingHolder(CubitString n, void (*setFn) (double),
			     double (*getFn)())
{
  setting_type = SETTING_DOUBLE_VALUE_TYPE;
  description = "";
  name = n;

  get_int_function    = NULL;
  set_int_function    = NULL;
  get_double_function = getFn;
  set_double_function = setFn;
  get_bool_function   = NULL;
  set_bool_function   = NULL;
  get_string_function = NULL;
  set_string_function = NULL;

}


// constructor for boolean settings
SettingHolder::SettingHolder(CubitString n, void (*setFn) (CubitBoolean),
			     CubitBoolean (*getFn)())
{
  setting_type = SETTING_CUBITBOOLEAN_VALUE_TYPE;
  description = "";
  name = n;

  get_int_function    = NULL;
  set_int_function    = NULL;
  get_double_function = NULL;
  set_double_function = NULL;
  get_bool_function   = getFn;
  set_bool_function   = setFn;
  get_string_function = NULL;
  set_string_function = NULL;

}

// constructor for string settings
SettingHolder::SettingHolder(CubitString n, void (*setFn) (CubitString),
			     CubitString (*getFn)())
{
  setting_type = SETTING_CUBITSTRING_VALUE_TYPE;
  description = "";
  name = n;

  get_int_function   = NULL;
  set_int_function   = NULL;
  get_double_function = NULL;
  set_double_function = NULL;
  get_bool_function   = NULL;
  set_bool_function   = NULL;
  get_string_function = getFn;
  set_string_function = setFn;

}


SettingHolder::~SettingHolder() {}




















