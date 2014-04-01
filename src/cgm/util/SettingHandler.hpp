//- Class: SettingHandler
//- Description: SettingHandler class - used for save and restore settings
//- Owner: Wesley Gill
//- Checked By:
//- Version

#ifndef SETTINGHANDLER_HPP
#define SETTINGHANDLER_HPP

#include "CubitString.hpp"
#include "CubitDefines.h"
#include "SettingHolder.hpp"
#include <map>
#include <vector>
#include "CubitUtilConfigure.h"


class CUBIT_UTIL_EXPORT SettingHandler 
{
  friend class SettingHolder; 
public:
  ~SettingHandler();
  static SettingHandler* instance(); // make an instance of SettingHandler

  static void delete_instance(); //deletes instance if it exists 

  void add_setting(char* name, void (*setFn) (int), int (*getFn)());
  void add_setting(char* name, void (*setFn) (double), double (*getFn)());
  void add_setting(char* name, void (*setFn) (CubitBoolean), CubitBoolean (*getFn)());
  void add_setting(char* name, void (*setFn) (CubitString), CubitString (*getFn)());

#ifdef BOYD15
  int num_settings();
#endif
  void get_settings_list( std::vector< std::pair<CubitString, SettingHolder*> > &list);

  void print_settings(); // Output all settings using PRINT_INFO
  void save_settings(const char *filename); //Save settings into a file specified by the user
  void save_settings(); //Save settings into a default file  
  void restore_settings(const char *filename); // Restore settings stored in a file

	SettingHolder* get_setting_holder(CubitString name);
	inline SettingType get_setting_type(SettingHolder *setting);

	inline int get_setting_int(SettingHolder *setting);
	inline double get_setting_double(SettingHolder *setting);
	inline CubitBoolean get_setting_bool(SettingHolder *setting);
	inline CubitString get_setting_string(SettingHolder *setting);
	inline int get_setting_debug_value(SettingHolder *setting);

	inline void set_setting_int(SettingHolder *setting, int int_value);
	inline void set_setting_double(SettingHolder *setting, double double_value);
	inline void set_setting_bool(SettingHolder *setting, CubitBoolean bool_value);
	inline void set_setting_string(SettingHolder *setting, CubitString string_value);
	inline void set_setting_debug_value(SettingHolder *setting, int debug_value);

private:

  void (*debug_set_function)(const int index, const int value);  //Pointer to the debug set function
  int  (*debug_get_function)(const int index);  //Pointer to the debug set function

  //Private constructor for singleton pattern
  SettingHandler();  

  static SettingHandler* instance_;
  
  //The debug flags are added to the map in the constuctor of the MessageFlag class.  Although 
  //the MessageFlag constructor should only be called once(Singleton) it is called multiple times.
  //This causes an error because you cannot add the same setting to the map more than once.  This boolean
  //prevents this problem from occuring.
  CubitBoolean debug_flags_added;

  std::map<CubitString, SettingHolder*> mSettingsList;


};

inline SettingType SettingHandler::get_setting_type(SettingHolder *setting)
{ return setting->setting_type; }

inline int SettingHandler::get_setting_int(SettingHolder *setting)
{ return (setting->get_int_function)(); }

inline double SettingHandler::get_setting_double(SettingHolder *setting)
{ return (setting->get_double_function)(); }

inline CubitBoolean SettingHandler::get_setting_bool(SettingHolder *setting)
{ return (setting->get_bool_function)(); }

inline CubitString SettingHandler::get_setting_string(SettingHolder *setting)
{ return (setting->get_string_function)(); }

//inline int SettingHandler::get_setting_debug_value(SettingHolder *setting)
//{ return debug_get_function(setting->index); }

inline void SettingHandler::set_setting_int(SettingHolder *setting, int int_value)
{ (setting->set_int_function)(int_value); }

inline void SettingHandler::set_setting_double(SettingHolder *setting, double double_value)
{ (setting->set_double_function)(double_value); }

inline void SettingHandler::set_setting_bool(SettingHolder *setting, CubitBoolean bool_value)
{ (setting->set_bool_function)(bool_value); }

inline void SettingHandler::set_setting_string(SettingHolder *setting, CubitString string_value)
{ (setting->set_string_function)(string_value); }

//inline void SettingHandler::set_setting_debug_value(SettingHolder *setting, int debug_value)
//{ (*debug_set_function)(setting->index, debug_value); }

#endif







