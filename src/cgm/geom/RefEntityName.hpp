//- Class: RefEntityName
//-
//- Description: This file maintains a map between a name and a
//- refentity. Names must be unique, but a refentity may have 
//- multiple names. This is a singleton class.
//-
//- Owner: Greg Sjaardema
//- Checked by: 
//- Version: $Id: 

#ifndef REFENTITYNAME_HPP
#define REFENTITYNAME_HPP

#include "RefEntityNameMap.hpp"
#include "SDLList.hpp"
#include "CubitGeomConfigure.h"

template <class X> class DLIList;
class CubitString;
class RefEntity;
class RefEntityNameMapList;
class RefEntity;

SDLListdeclare(RefEntityNameMapList,RefEntityNameMap*,key,CubitString)

class CUBIT_GEOM_EXPORT RefEntityName
{

public:

  //Initialize all settings in this class
  static void initialize_settings();
  static RefEntityName *instance();
  ~RefEntityName();

  void       remove_refentity_name(const CubitString &name,
                                   CubitBoolean update_attribs = CUBIT_TRUE);
  //- Remove the this name from the map. All other names for this
  //- RefEntity remain in the map

  void       remove_refentity_name(RefEntity *entity,
                                   CubitBoolean update_attribs = CUBIT_TRUE);
  //- Remove this RefEntity and all its names from the map.

  CubitStatus add_refentity_name(RefEntity *entity,
                                 CubitString &name,
                                 CubitBoolean update_attribs = CUBIT_TRUE);
  //- Add the map between this RefEntity and the name.
  //- If 'name' is invalid or duplicate, it is changed and
  //- new name is returned.
  //- Returns CUBIT_FAILURE if unsuccessful.

  CubitStatus add_refentity_name(RefEntity *entity,
                                 DLIList<CubitString*> &names,
                                 CubitBoolean update_attribs = CUBIT_TRUE);
    //- add multiple names to the name map for this entity
  
  int        get_refentity_name(const RefEntity *entity,
				DLIList<CubitString*> &names,
				int get_only_one_name = CUBIT_FALSE);
  //- Get the list of names corresponding to this RefEntity.
  //- Returns the number of names found for this RefEntity.

  int        get_refentity_name(const RefEntity *entity,
                                RefEntityNameMapList &names,
                                int get_only_one_name = CUBIT_FALSE);
  //- Get the list of names corresponding to this RefEntity.
  //- Returns the number of names found for this RefEntity.

  void       merge_refentity_names(RefEntity *retained,
                                   RefEntity *dead);
  //- Associates the names of the dead RefEntity with the retained
  //- RefEntity.

  void       switch_refentity_names(RefEntity *entity1,
                                    RefEntity *entity2);
  //- Switch names of entity1 and entity2

  RefEntity* get_refentity(const CubitString &name);
  //- Returns a pointer to the RefEntity with this name.

  void       list_refentity_names(const char *type = "all");
  //- Prints an alphabetical list of all entity names in the map.
  //- The argument is a filter on the list which restricts it to 
  //- certain entity types (body, volume, surface, curve, vertex).

  const char*      get_refentity_type(const CubitString &name);
  //- Returns a pointer to the string identifying the type of the
  //- RefEntity with this name. (Body, Volume, ...)

  int        get_refentity_id(const CubitString &name);
  //- Returns the id of the RefEntity with this name. 

  static int  get_generate_default_names();
  static void set_generate_default_names(int on_off);
  //- Specifies whether default RefEntity names should be generated (Surface1). 

  static int  get_fix_duplicate_names();
  static void set_fix_duplicate_names(int on_off);
  //- Specifies whether duplicate names should be fixed (made unique)

  void set_character(char rep, const CubitString &type);
  char get_character(const CubitString &type) const;


  static void set_suffix_setting(CubitString rep);
  static CubitString get_suffix_setting();
  static void set_replacement_setting(CubitString rep);
  static CubitString get_replacement_setting();
  
  static CubitString base_name(const CubitString *name);
  
  CubitBoolean same_base_name(const CubitString *name1,
                              const CubitString *name2);
    //- return CUBIT_TRUE if the base name (i.e. the name before the suffixChar)
    //- is the same

  static CubitBoolean get_merge_base_names() {return mergeBaseNames;};
  static void set_merge_base_names(CubitBoolean flag) {mergeBaseNames = flag;};
    //- get/set mergeBaseNames
    
  void copy_refentity_names( const RefEntity* source, 
                             RefEntity* target,
                             CubitBoolean update_attribs = CUBIT_TRUE );
  void copy_refentity_names( DLIList<RefEntity*>& source_list,
                             RefEntity* target,
                             CubitBoolean unique_base_names = CUBIT_TRUE,
                             CubitBoolean update_attribs = CUBIT_TRUE );
    //- Copy names from one entity to another.
  
protected:
  RefEntityName();
  //- Class Constructor. (Not callable by user code. Class is constructed
  //- by the {instance()} member function.
  
private:
  CubitStatus clean(CubitString &raw_name);
  //- replace invalid characters in {raw_name} with {replacementCharacter}

  CubitStatus generate_unique_name(CubitString &name);

  static RefEntityName *instance_;
  //- static pointer to unique instance of this class
  
  RefEntityNameMapList nameEntityList;
  //- A sorted DLList used to store the association between RefEntity names
  //- and the corresponding RefEntity. 

  static int fixDuplicateNames;
  static int generateDefaultNames;
  static char replacementCharacter;
  static char suffixCharacter;
 
  static CubitBoolean mergeBaseNames;
    //- if CUBIT_TRUE, when entity names are merged, names with duplicate
    //- base names are removed
};
#endif


