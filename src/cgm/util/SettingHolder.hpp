//- Class: SettingHolder
//- Description: SettingHolder class - used for save and restore settings
//- Owner: Wesley Gill
//- Checked By:
//- Version

#ifndef SETTINGHOLDER_HPP
#define SETTINGHOLDER_HPP

#include "CubitString.hpp"
#include "CubitDefines.h"
#include "CubitUtilConfigure.h"

typedef enum SETTING_TYPE
{
  SETTING_UNKNOWN_VALUE_TYPE = -1,
  SETTING_INT_VALUE_TYPE = 0,
  SETTING_DOUBLE_VALUE_TYPE = 1,
  SETTING_CUBITBOOLEAN_VALUE_TYPE = 2,
  SETTING_CUBITSTRING_VALUE_TYPE = 3
                        // int = 0, double = 1, CubitBoolean = 2, CubitString = 3, debug = 4
} SettingType;


class CUBIT_UTIL_EXPORT SettingHolder
{
  friend class SettingHandler;
public:
  SettingHolder();
  ~SettingHolder();
 
  //- constructor for int settings
  SettingHolder(CubitString, void (*)(int), int (*)(), int);
  
  //- constructor for double settings
  SettingHolder(CubitString, void (*)(double), double (*)(), double);
 
  //- constructor for CubitBoolean settings
  SettingHolder(CubitString, void (*setFn)(CubitBoolean),
                CubitBoolean (*getFn)(), CubitBoolean);

  //- constructor for CubitString settings
  SettingHolder(CubitString n, void (*setFn) (CubitString), 
                CubitString (*getFn)(), char* value);

  //- constructor for int settings
  SettingHolder(CubitString, void (*)(int), int (*)());
  
  //- constructor for double settings
  SettingHolder(CubitString, void (*)(double), double (*)());
 
  //- constructor for CubitBoolean settings
  SettingHolder(CubitString, void (*)(CubitBoolean), CubitBoolean (*)());

  //- constructor for CubitString settings
  SettingHolder(CubitString n, void (*setFn) (CubitString),
			     CubitString (*getFn)());


private:

  //Use an enum for setting_type??

  SettingType setting_type;     // Defines the type of setting 
                        // int = 0, double = 1, CubitBoolean = 2, CubitString = 3, debug = 4
  
  CubitString name;               // Name of the setting
  CubitString description;        // Description of the setting
  int int_value;                    
  double double_value;
  CubitString string_value;

  int (*get_int_function) ();            // Pointer to get int function
  double (*get_double_function) ();      // Pointer to get double function
  CubitBoolean (*get_bool_function) ();  // Pointer to get boolean function
  CubitString (*get_string_function) (); // Pointer to get string function

  void (*set_int_function) (int);            // Pointer to set int function
  void (*set_double_function) (double);      // Pointer to set double function
  void (*set_bool_function) (CubitBoolean);  // Pointer to set bool function
  void (*set_string_function) (CubitString); // Pointer to set string function 

};

#endif





