#ifndef __CUBITATTRIBMANAGER_HPP__
#define __CUBITATTRIBMANAGER_HPP__

class RefEntity;
class CubitAttrib;
class CubitAttribUser;
class CubitSimpleAttrib;

#include "CubitDefines.h"
#include "DLIList.hpp"
#include "CubitString.hpp"
#include "CubitGeomConfigure.h"

#if defined(WIN32)
#pragma warning(push)
#pragma warning(disable : 4251)  // hide warnings about template dll exports
#endif

#include <map>

typedef CubitAttrib* (*CACreateFunction)(RefEntity*, CubitSimpleAttrib*);

class CUBIT_GEOM_EXPORT CubitAttribManager
{
public:
  CubitAttribManager();
  virtual ~CubitAttribManager();

  // stores factory function and settings for an attribute
  // returns a key that can be used to refer to the attribute type to avoid search by name
  CubitStatus register_attrib_type( const int att_type,
                                     const char* att_type_name,
                                     const char* att_internal_name,
                                     CACreateFunction p_create_function,
                                     CubitBoolean auto_actuate_flag,
                                     CubitBoolean auto_update_flag,
                                     CubitBoolean auto_write_flag,
                                     CubitBoolean auto_read_flag,
                                     CubitBoolean actuate_in_constructor,
                                     CubitBoolean actuate_after_geom_changes);
// actuateInConstructor flag: if this flag is set to true,
// this attribute will be actuated in a constructor, otherwise
// it will be actuated after a geometry operation or geometry
// import; usually, attributes which need full topological information
// must be actuated after the geometry operation is finished
//
// actuateAfterGeomChanges flag: if this flag is set to true,
// this attribute will be actuated only after all changes to geometry
// have been done.  Implies that if an attribute changes geometry,
// or if it doesn't matter, this flag should be false (the default value); 
// if it requires the geometry to be in its final state, this flag 
// should be true

  void get_registered_types(DLIList<int> &types);

  CubitAttrib *create_cubit_attrib(const int attrib_type,
                                   RefEntity *attrib_user,
                                   CubitSimpleAttrib *p_csa);

  const char * att_internal_name(int attrib_type);

  CubitBoolean auto_actuate_flag(int attrib_type);
  CubitStatus set_auto_actuate_flag(int attrib_type, CubitBoolean value);
  void set_all_auto_actuate_flags(CubitBoolean value);

  CubitBoolean auto_update_flag(int attrib_type);
  CubitStatus set_auto_update_flag(int attrib_type, CubitBoolean value);
  void set_all_auto_update_flags(CubitBoolean value);

  CubitBoolean auto_write_flag(int attrib_type);
  CubitStatus set_auto_write_flag(int attrib_type, CubitBoolean value);
  void set_all_auto_write_flags(CubitBoolean value);

  CubitBoolean auto_read_flag(int attrib_type);
  CubitStatus set_auto_read_flag(int attrib_type, CubitBoolean value);
  void set_all_auto_read_flags(CubitBoolean value);

  CubitBoolean actuate_in_constructor(int attrib_type);
  CubitBoolean actuate_after_geom_changes(int attrib_type);

  CubitStatus auto_update_attribs(RefEntity *cau);

//  CubitStatus actuate_list(DLIList<CubitAttribUser*> entity_list);

  int attrib_type(CubitSimpleAttrib *csa_ptr);
  int attrib_type(const char* name);
  int attrib_type_from_internal_name(const char* name);

//  void clear_attrib_importeds() ;
//  void report_attrib_importeds();

  int auto_flag();
  void auto_flag(int flag);

  bool silent_flag() {return silentFlag;}
  void silent_flag(const bool flag) {silentFlag = flag;}

  const char * att_name(int attrib_type);

private:
/*
  class CARegisterData
  {
  public:
    CARegisterData(const char *type_name, const char *internal_name,
                   CACreateFunction p_func, unsigned int flags)
    {
      mTypeName = type_name;
      mInternalName = internal_name;
      mCreateFunction = p_func;
      mFlags = flags;
    }

    ~CARegisterData();

    CubitString mTypeName;
    CubitString mInternalName;
    CACreateFunction mCreateFunction;
    unsigned int mFlags;
  };

  #define CA_AUTOACTUATE             0x01
  #define CA_AUTOUPDATE              0x02
  #define CA_AUTOWRITE               0x04
  #define CA_AUTOREAD                0x08
  #define CA_ACTUATEINCONSTRUCTOR    0x10
  #define CA_ACTUATEAFTERGEOMCHANGES 0x20

  std::map<int, CARegisterData*> mTypeToCAData;
//struct ltstr
//{
//  bool operator()(const char* s1, const char* s2) const
//  {
//    return strcmp(s1, s2) < 0;
//  }
//};
//  std::map<const char*, int, ltstr> mTypeNameToType;
//  std::map<const char*, int, ltstr> mInternalNameToType;
*/

  DLIList<int> mTypes;
  DLIList<char*> mTypeNames;
  DLIList<char*> mInternalNames;
  DLIList<CACreateFunction> mCreatorFunctions;
  DLIList<CubitBoolean> mAutoActuateFlags;
  DLIList<CubitBoolean> mAutoUpdateFlags;
  DLIList<CubitBoolean> mAutoWriteFlags;
  DLIList<CubitBoolean> mAutoReadFlags;
  DLIList<CubitBoolean> mActuateInConstructor;
  DLIList<CubitBoolean> mActuateAfterGeomChanges;

  // TODO - get rid of mAttribImported
//  DLIList<CubitBoolean> mAttribImported;

  bool silentFlag;

};

#if defined(WIN32)
#pragma warning(pop)
#endif

#endif //__CUBITATTRIBMANAGER_HPP__

