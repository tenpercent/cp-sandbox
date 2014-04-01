#include "RefEntityName.hpp"
#include "RefEntityNameMap.hpp"
#include "CubitString.hpp"
#include "RefEntity.hpp"
#include "CubitMessage.hpp"
#include "CubitUtil.hpp"
#include "CastTo.hpp"
#include "RefGroup.hpp"
#include "CubitAttribUser.hpp"
#include "CADefines.hpp"

#include "SettingHandler.hpp"

#include <ctype.h>
#include <assert.h>

RefEntityName* RefEntityName::instance_ = 0;
CubitBoolean RefEntityName::mergeBaseNames = CUBIT_FALSE;
int RefEntityName::generateDefaultNames = CUBIT_FALSE;
int RefEntityName::fixDuplicateNames = CUBIT_TRUE;
char RefEntityName::suffixCharacter = '@';
char RefEntityName::replacementCharacter = '_';
static int is_valid_char(char c);
static int is_valid_first_char(char c);
 
RefEntityName* RefEntityName::instance()
{
  if (instance_ == 0) {
    instance_ = new RefEntityName;
  }
  return instance_;
}

RefEntityName::RefEntityName()
{
  generateDefaultNames = CUBIT_FALSE;
  fixDuplicateNames = CUBIT_TRUE;
  replacementCharacter = '_';
  suffixCharacter = '@';
}

RefEntityName::~RefEntityName()
{
  instance_ = 0;
  for (int i=nameEntityList.size(); i > 0; i--) {
    delete nameEntityList.nullify();
    nameEntityList.step();
  }
  nameEntityList.compress();
}

void RefEntityName::remove_refentity_name(const CubitString &name,
                                          CubitBoolean update_attribs)
{
  if (nameEntityList.move_to(name))
  {
      // get the map entity's ref entity
    RefEntityNameMap *name_map = nameEntityList.remove();
    RefEntity *entity = name_map->value();
    
      // tell the ref entity to update its names
    if (update_attribs)
    {
      CubitAttrib *attrib = entity->get_cubit_attrib(CA_ENTITY_NAME);
      attrib->has_updated(CUBIT_FALSE);
      attrib->update();
    }

    delete name_map;
  }
}

void RefEntityName::remove_refentity_name(RefEntity *entity,
                                          CubitBoolean update_attribs)
{
    // NOTE: There may be multiple names for one RefEntity. Make sure to 
    //       remove all of them. 
  
  int found_match = 0;
  RefEntityNameMapList old_maps;
  
  int i;
  for ( i=nameEntityList.size(); i > 0; i--)
  {
    if (nameEntityList.get()->value() == entity)
    {
      old_maps.append(nameEntityList.nullify());
      found_match++;
    }
    nameEntityList.step();
  }
  nameEntityList.compress();
  
    // now tell the entity to update its name attribute
  if (update_attribs)
  {
    CubitAttrib *attrib = entity->get_cubit_attrib(CA_ENTITY_NAME,
                                                   CUBIT_FALSE);
    if (attrib) {
      attrib->has_updated(CUBIT_FALSE);
      attrib->update();
    }
  }

  for (i = old_maps.size(); i > 0; i--)
    delete old_maps.get_and_step();
  
}

CubitStatus RefEntityName::add_refentity_name(RefEntity *entity,
				      DLIList<CubitString*> &names,
                                              CubitBoolean update_attribs)
{
  names.reset();
  //int num_new_names = names.size();

  DLIList<CubitString*> new_names;
  
  int i;
  for ( i = names.size(); i > 0; i--)
  {
    CubitString name = *names.get();
    CubitString in_name = name;
    CubitBoolean warn_name_change = CUBIT_FALSE;
    
      // first, clean the name
    if (clean(name))
    {
//       PRINT_WARNING("Entity name '%s' is invalid. Changed to '%s'\n",
//                     in_name.c_str(), name.c_str());
      warn_name_change = CUBIT_TRUE;
    }
    
      // now, check for valid name
    CubitBoolean name_valid = CUBIT_FALSE;
    
    if (name == "")
    {
        // blank name entered - do nothing
    }
    
    else if (nameEntityList.move_to(name) &&
             nameEntityList.get()->value() == entity)
    {
        // Tried to assign same name to entity
      if ( DEBUG_FLAG(92) ) 
      {
          // check to see if it's the same as this entity's default name,
          // if so, it probably came in on an attribute, and we don't need
          // to hear about it; otherwise, write the warning
        CubitString def_name;
        entity->generate_default_name(def_name);
        if (name != def_name)
          PRINT_INFO("Entity name '%s' already assigned to %s %d\n",
                     name.c_str(), 
                     entity->class_name(), entity->id());
      }
    }
    else if (nameEntityList.move_to(name) &&
             nameEntityList.get()->value() != entity)
    {
        // Tried to assign existing name to another entity
      PRINT_DEBUG_92( "Entity name '%s' for %s %d is already used by %s %d\n",
                  name.c_str(), entity->class_name(), entity->id(),
                  nameEntityList.get()->value()->class_name(),
                  nameEntityList.get()->value()->id());
      
        // either we fix it and keep it, or we don't and get rid of it
      name_valid = CUBIT_FALSE;
      if (get_fix_duplicate_names())
      {
        if (generate_unique_name(name))
        {
          PRINT_DEBUG_92( "\t%s %d name changed to '%s'\n",
                          entity->class_name(), entity->id(), name.c_str());
          if(warn_name_change)
          {
            PRINT_WARNING("Entity name '%s' is invalid. Changed to '%s'\n",
                          in_name.c_str(), name.c_str());
          }
          
          name_valid = CUBIT_TRUE;
        }
      }
    }
    else
    {
      if(warn_name_change)
      {
        PRINT_WARNING("Entity name '%s' is invalid. Changed to '%s'\n",
                      in_name.c_str(), name.c_str());
      }
      
        // else the name must be valid
      name_valid = CUBIT_TRUE;
    }
    
    if (name_valid == CUBIT_TRUE)
    {
        // name is valid
      if (name != in_name)
          // name was changed; change in name list too
        *names.get() = name;

        // save this name to later
      new_names.append(names.get());
      
        // now step the list
      names.step();
    }
  }
  
  if (new_names.size() > 0)
  {
      // there are some valid, new names; add them, then update attribute
    new_names.reset();
    
    CubitString name;
    for (i = new_names.size(); i > 0; i--)
    {
      name = *new_names.get_and_step();
      if (nameEntityList.move_to(name) &&
          nameEntityList.get()->value() == entity) {
            PRINT_DEBUG_92("Already have name %s for %s %d.\n",
                           name.c_str(), entity->class_name(), entity->id());
      }
      
      else {
        nameEntityList.insert(new RefEntityNameMap(name, entity));
      }
    }
    
    if (update_attribs == CUBIT_TRUE)
    {
        // now tell the entity to update its name attribute
      CubitAttrib *attrib = entity->get_cubit_attrib(CA_ENTITY_NAME);
        // force update by resetting update flag
      attrib->has_updated(CUBIT_FALSE);
      attrib->update();
    }
  }
  
  return CUBIT_SUCCESS;
}

CubitStatus RefEntityName::add_refentity_name(RefEntity *entity,
				      CubitString &name,
                                              CubitBoolean update_attribs)
{
  if (name == "")
    return CUBIT_FAILURE;
  
  CubitString in_name = name;
  bool warn_name_change = false;
  
  if (clean(name))
  {
    warn_name_change = true;
  }
  
  if (nameEntityList.move_to(name))
  {
    RefEntity *old_entity = nameEntityList.get()->value();
    if (old_entity == entity)
    {
        // Tried to assign same name to entity
      if ( DEBUG_FLAG(92) ) 
      {
        PRINT_INFO("Entity name '%s' already assigned to %s %d\n",
                   name.c_str(), 
                   entity->class_name(), entity->id());
        return CUBIT_FAILURE;
      }
      return CUBIT_SUCCESS;
    }
    else
    {
        // Tried to assign existing name to another entity
      if ( DEBUG_FLAG(92) )
        PRINT_WARNING("Entity name '%s' for %s %d is already used by %s %d\n",
                      name.c_str(),
                      entity->class_name(), entity->id(),
                      old_entity->class_name(), old_entity->id());
      if (get_fix_duplicate_names())
      {
        if (generate_unique_name(name))
        {
          if (warn_name_change)
          {
            PRINT_WARNING("Entity name '%s' is invalid. Changed to '%s'\n",
              in_name.c_str(), name.c_str());
          }
          if ( DEBUG_FLAG(92) )
            PRINT_WARNING("\t%s %d name changed to '%s'\n",
            entity->class_name(), entity->id(), name.c_str());
          return add_refentity_name(entity, name);
        }
      }
      return CUBIT_FAILURE;
    }
  }

  if (warn_name_change)
  {
    PRINT_WARNING("Entity name '%s' is invalid. Changed to '%s'\n",
      in_name.c_str(), name.c_str());
  }

  RefEntityNameMap *entity_name =
    new RefEntityNameMap(name, entity);
  nameEntityList.insert(entity_name);
  
  if (update_attribs == CUBIT_TRUE)
  {
      // now tell the entity to update its name attribute
    CubitAttrib *attrib = entity->get_cubit_attrib(CA_ENTITY_NAME);
      // force update by resetting update flag
    attrib->has_updated(CUBIT_FALSE);
    attrib->update();
  }
  
  return CUBIT_SUCCESS;
}

int RefEntityName::get_refentity_name(const RefEntity *entity,
                                      DLIList<CubitString*> &names,
                                      int get_only_one)
{
    // NOTE: There may be multiple names for one RefEntity. Make sure to 
    //       access all of them. 

  RefEntityNameMapList temp_list;
  int found_match = get_refentity_name(entity, temp_list, get_only_one);
  
  for (int i=temp_list.size(); i > 0; i--)
    names.append(temp_list.get_and_step()->ptr_to_key());

  return found_match;
}

int RefEntityName::get_refentity_name(const RefEntity *entity,
                                      RefEntityNameMapList &names,
                                      int get_only_one)
{
    // NOTE: There may be multiple names for one RefEntity. Make sure to 
    //       access all of them. 
  
  nameEntityList.reset();
  int found_match = 0;
  for (int i=nameEntityList.size(); i > 0; i--)
  {
    if (nameEntityList.get()->value() == entity)
    {
      names.append(nameEntityList.get());
      found_match++;
      if (get_only_one)
        return found_match;
    }
    nameEntityList.step();
  }
  return found_match;
}

RefEntity*  RefEntityName::get_refentity(const CubitString &name) 
{
  if (nameEntityList.move_to(name)) {
    RefEntity *entity = nameEntityList.get()->value();
    return entity;
  } else {
    return NULL;
  }
}

const char* RefEntityName::get_refentity_type(const CubitString &name)
{
  RefEntity *entity = get_refentity(name);
  if (entity != NULL) {
    return entity->class_name();
  } else {
    return NULL;
  }
}

int RefEntityName::get_refentity_id(const CubitString &name)
{
  RefEntity *entity = get_refentity(name);
  if (entity != NULL) {
    return entity->id();
  } else {
    return 0;
  }
}

void RefEntityName::merge_refentity_names(RefEntity *retained,
                                          RefEntity *dead)
{
  // NOTE: There may be multiple names for one RefEntity. Make sure to 
  //       process all of them.

  //CubitBoolean was_named = CUBIT_FALSE;
  int i;

  RefEntityNameMapList combined_names;
  if (get_merge_base_names() == CUBIT_TRUE)
    get_refentity_name(retained, combined_names);
  
  get_refentity_name(dead, combined_names);

  if (combined_names.size() == 0) return;

  if (get_merge_base_names() == CUBIT_TRUE) {
      // remove names with identical base names, keeping lowest
      // first, sort the list
    combined_names.sort();
      // now, look at each name; if its previous neighbor has the same base name,
      // change the owner to the dead entity; dead entity's names will all get
      // removed afterwards
    combined_names.reset();

      // make sure we keep the first base name
    combined_names.get()->value(retained);
    combined_names.step();
    for (i = combined_names.size()-1; i > 0; i--) {
      if (same_base_name(combined_names.prev()->ptr_to_key(),
                         combined_names.get()->ptr_to_key())) {

        combined_names.get()->value(dead);
      }
      else {
        combined_names.get()->value(retained);
      }
      combined_names.step();
    }
      // now, remove the dead names by removing names for the dead entity;
      // don't update attributes, done later
    remove_refentity_name(dead);
  }
  else {
    for (i = combined_names.size(); i > 0; i--)
      combined_names.get_and_step()->value(retained);
  }
  
        // now tell the entity to update its name attribute
  CubitAttrib *attrib = retained->get_cubit_attrib(CA_ENTITY_NAME);
  attrib->has_updated(CUBIT_FALSE);
  attrib->update();
}

void RefEntityName::switch_refentity_names(RefEntity *entity1,
                                           RefEntity *entity2)
{
  // NOTE: There may be multiple names for one RefEntity. Make sure to 
  //       process all of them. 

  for (int i=nameEntityList.size(); i > 0; i--) {
    if (nameEntityList.get()->value() == entity1)
      nameEntityList.get()->value(entity2);
    else if (nameEntityList.get()->value() == entity2)
      nameEntityList.get()->value(entity1);

    nameEntityList.step();
  }

    // now tell the entities to update their name attribute; if either
    // entity now has no name, the attributes are removed in the proper manner
  CubitAttrib *attrib = entity1->get_cubit_attrib(CA_ENTITY_NAME);
  attrib->update();
  attrib = entity2->get_cubit_attrib(CA_ENTITY_NAME);
  attrib->update();
}

void RefEntityName::list_refentity_names(const char *type)
{
     // 'type' is a lowercase string specifying the type of entity to be
     // returned. (surface, body, ...)
     // if 'type' == "all", then all names will be listed.
   
   int do_all = (strcmp("all", type) == 0 || strlen(type) == 0);
   
   PRINT_INFO("\t______Name______  \t__Type__ Id\t\t_Propagated_\n");
   RefEntity *entity;
   nameEntityList.reset();

   if(nameEntityList.size() > 0){
      for(int j = nameEntityList.size(); j > 0; j--)
      {
         entity = nameEntityList.get()->value();
         if (do_all || CubitUtil::strcmp_case_insensitive(entity->class_name(),
                                                       type) == 0)
      {
         PRINT_INFO("%24s  %8s\t%-10d \tNo\n", nameEntityList.get()->key().c_str(),
                    entity->class_name(), entity->id());
         }
         nameEntityList.step();
      }
   }

}

int RefEntityName::get_generate_default_names()
{
  return generateDefaultNames;
}

void RefEntityName::set_generate_default_names(int on_off)
{
  generateDefaultNames = on_off;
}

int RefEntityName::get_fix_duplicate_names()
{
  return fixDuplicateNames;
}

void RefEntityName::set_fix_duplicate_names(int on_off)
{
  fixDuplicateNames = on_off;
}

CubitStatus RefEntityName::clean(CubitString &raw_name)
{

  if (raw_name == "") return CUBIT_FAILURE;
  
    // A valid name consists of alphanumeric characters plus '.', '_', '-', or '@'
  CubitStatus found_invalid_character = CUBIT_FAILURE;

  // Initial character must be alphabetic or "_".
  char c = raw_name.get_at(0);
  if (!is_valid_first_char(c)) {
    if (is_valid_first_char(get_character("replace")))
      raw_name.put_at(0, get_character("replace"));
    else
      raw_name.put_at(0, '_');
    found_invalid_character = CUBIT_SUCCESS;
  }
  
  for (unsigned int i = 1; i < raw_name.length(); i++) {
    c = raw_name.get_at(i);
    if (!is_valid_char(c)) {
      found_invalid_character = CUBIT_SUCCESS;
      raw_name.put_at(i, get_character("replace"));
    }
  }
  return found_invalid_character;
}

void RefEntityName::set_character(char rep, const CubitString &type)
{
  if (is_valid_char(rep)) {
    if (type.get_at(0) == 'R' || type.get_at(0) == 'r')
      replacementCharacter = rep;
    else if (type.get_at(0) == 'S' || type.get_at(0) == 's')
      suffixCharacter = rep;
    else
      PRINT_ERROR("Invalid character type '%s', must be "
		  "'replacement' or 'suffix'\n", type.c_str());
  } else {
    PRINT_ERROR("Character '%c' is not a valid entity name character\n",
		rep);
  }
}

char RefEntityName::get_character(const CubitString& type) const
{
  if (type.get_at(0) == 'R' || type.get_at(0) == 'r')
    return replacementCharacter;
  else if (type.get_at(0) == 'S' || type.get_at(0) == 's')
    return suffixCharacter;  
  else {
    PRINT_ERROR("Invalid character type '%s', must be "
		"'replacement' or 'suffix'\n", type.c_str());
    return '#';
  }
}

static int is_valid_char(char c)
{
  return (isalnum(c) || c == '.' || c == '_' || c == '@' || c == '-');
}

static int is_valid_first_char(char c)
{
  return (isalpha(c) || c == '_');
}

CubitStatus RefEntityName::generate_unique_name(CubitString &name)
{
  // The method used to generate a unique name is to append
  // 'suffixCharacter' and
  // a letter from A-Z, a-z, or 0-9 to the end of name.
  // If none of these produce a unique name, CUBIT_FALSE is returned.

  CubitString alphabet =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
  CubitString suffix("  ");
  suffix.put_at(0, suffixCharacter);
  
  CubitString internal = name;

  // See if there is an suffixCharacter sign in name already. If so, is it
  // the second to the last character in the string.
  if (name.length() < 2 ||
      name.get_at(name.length()-2) != suffixCharacter) {
    // Name does not contain suffixCharacter at correct location, add one on.
    internal += suffix;
  }
  
  CubitStatus found_unique = CUBIT_FAILURE;
  int continue_trying = CUBIT_TRUE;

  while (!found_unique && continue_trying) {
	 
    continue_trying = CUBIT_FALSE;
    int name_length = internal.length();
    unsigned int number_tested = 0;
    for (unsigned int i=0; i < alphabet.length(); i++) {
      internal.put_at(name_length-1, (char)alphabet.get_at(i));
      if (!nameEntityList.move_to(internal)) {
	found_unique = CUBIT_SUCCESS;
	break;
      }
      number_tested++;
    }

    if (number_tested == alphabet.length()) {
      // All suffixes used. Add another suffixCharacter and try again
      // Name will look like 'Name@@1' or 'Name@@@1'...
      // Find LAST suffixCharacter in name
      int ch;
      for (ch = (int)(internal.length())-1; ch >= 0; ch--) {
	if (internal.get_at(ch) == suffixCharacter) {
	  break;
	}
      }
      if (internal.get_at(ch) == suffixCharacter) {
	// Add another suffixCharacter at ch+1
	// Assured that position ch+1 exists or we wouldn't be here
	internal.put_at(ch+1, suffixCharacter);
	if (ch+2 < (int)internal.length())
	  internal.put_at(ch+2, ' ');
	else
	  internal += " ";
	continue_trying = CUBIT_TRUE;
      }
    }
  }
  
  if (found_unique)
    name = internal;

  return found_unique;
}

CubitString RefEntityName::base_name(const CubitString* name)
{
  const char *pos = strchr(name->c_str(), suffixCharacter);
  if (!pos)
    return *name;
  return name->substr(0, pos - name->c_str());
}


// returns CUBIT_TRUE if the two names have the same base name (i.e.
// the same name before the first suffixCharacter
CubitBoolean RefEntityName::same_base_name(const CubitString *name1,
                                           const CubitString *name2)
{
  const char *pos1 = strchr(name1->c_str(), suffixCharacter);
  const char *pos2 = strchr(name2->c_str(), suffixCharacter);

    // check for replacement character in one but not the other
  int length1, length2;
  if (pos1 == NULL)
    length1 = strlen(name1->c_str());
  else
    length1 = pos1 - name1->c_str();

  if (pos2 == NULL)
    length2 = strlen(name2->c_str());
  else
    length2 = pos2 - name2->c_str();

    // if the lengths are different, the base names are also different
  if (length1 != length2)
    return CUBIT_FALSE;

  if (strncmp(name1->c_str(), name2->c_str(), length1) == 0)
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}

// Copy names from one entity to another.  
// Added Aug.14,2000 by J.Kraftcheck for propogating
// names after virtual geometry operations.
void RefEntityName::copy_refentity_names( const RefEntity *source,
                                          RefEntity *target,
                                          CubitBoolean update_attribs )
{
  //No NULL pointers.
  assert( source && target );
  
  //If we can't have duplicate names, then it is not possible to 
  //copy names.
  if( ! get_fix_duplicate_names() ) return;
  
  //Assume the name is valid already, as it is attached to
  //the source entity.  Also, assume the name is unique to
  //the source entity.
  DLIList<CubitString*> names;
  get_refentity_name( source, names );
  names.reset();
  
  //For each of the names on the source entity
  for( int i = names.size(); i > 0; i-- )
  {  
    CubitString name = *names.get_and_step();
    //make the name unique
    generate_unique_name( name );
    //associate name with target
    nameEntityList.insert( new RefEntityNameMap( name, target ) );
  }
    
  if( (names.size() > 0) && (update_attribs == CUBIT_TRUE) )
  {
      // now tell the entity to update its name attribute
    CubitAttrib *attrib = target->get_cubit_attrib( CA_ENTITY_NAME );
      // force update by resetting update flag
    attrib->has_updated( CUBIT_FALSE );
    attrib->update();
  }
}



//Initialize all settings in this class
void RefEntityName::initialize_settings()
{

  SettingHandler::instance()->add_setting("Default Names", 
                                          RefEntityName::set_generate_default_names, 
					  RefEntityName::get_generate_default_names); 

  SettingHandler::instance()->add_setting("Fix Duplicate Names", 
                                          RefEntityName::set_fix_duplicate_names, 
					  RefEntityName::get_fix_duplicate_names);
   
  SettingHandler::instance()->add_setting("Replacement Character", 
					  RefEntityName::set_replacement_setting, 
					  RefEntityName::get_replacement_setting);  

  SettingHandler::instance()->add_setting("Suffix Character",
					  RefEntityName::set_suffix_setting,
					  RefEntityName::get_suffix_setting);

  SettingHandler::instance()->add_setting("Merge Base Names", 
                                          RefEntityName::set_merge_base_names, 
					  RefEntityName::get_merge_base_names); 

}


void RefEntityName::copy_refentity_names( DLIList<RefEntity*>& source_list,
                                          RefEntity* target,
                                          CubitBoolean unique_base_names,
                                          CubitBoolean update_attribs )
{
  //Validate input.
  assert( target != 0 );
  if( source_list.size() < 1 ) return;
  
  //If we can't have duplicate names, we can't copy.
  if( ! get_fix_duplicate_names() ) return;
  
  //If we need to ensure no duplicate base names...
  if( unique_base_names )
  {
    //Get a big list of all the names from the source_list
    RefEntityNameMapList name_list;
    source_list.reset();
    for( int i = source_list.size(); i > 0; i-- )
      get_refentity_name( source_list.get_and_step(), name_list );
      
    //If we don't have any names to add, we're done
    if( name_list.size() < 1 ) return;
    
    //Get the first name
    name_list.sort();
    name_list.reset();
    CubitString* prev_key = name_list.get_and_step()->ptr_to_key();
    
    //Add name to target
    CubitString name = *prev_key;
    generate_unique_name( name );
    nameEntityList.insert( new RefEntityNameMap( name, target ) );
  
    //For the rest of the names...
    for( int j = name_list.size(); j > 1; j-- )
    {
      CubitString* key_ptr = name_list.get_and_step()->ptr_to_key();
      
      //If the name has a different base name than the previous, 
      //add it to the target.
      if( !same_base_name( prev_key, key_ptr ) )
      {
        name = *key_ptr;
        generate_unique_name( name );
        nameEntityList.insert( new RefEntityNameMap( name, target ) );
      }
       
      prev_key = key_ptr;
    }
  
    //update attribute, if asked to do so
    if( update_attribs  )
    {
      CubitAttrib *attrib = target->get_cubit_attrib( CA_ENTITY_NAME );
      attrib->has_updated( CUBIT_FALSE );
      attrib->update();
    }
    
  }
  else
  {
    //If we don't care about unique base names, just add them all.
    DLIList<CubitString*> name_list;
    for( int i = source_list.size(); i > 0; i-- )
      get_refentity_name( source_list.get_and_step(), name_list );
    name_list.reset();
    add_refentity_name( target, name_list, update_attribs );
  }
}



void RefEntityName::set_suffix_setting(CubitString rep)
{
  const char* tmp = rep.c_str();
  suffixCharacter = tmp[0];
}

CubitString RefEntityName::get_suffix_setting()
{
  char tmp[] = {suffixCharacter, '\0'};
  return CubitString(tmp);
}

void RefEntityName::set_replacement_setting(CubitString rep)
{
  const char* tmp = rep.c_str();
  replacementCharacter = tmp[0];
}

CubitString RefEntityName::get_replacement_setting()
{
  char tmp[] = {replacementCharacter, '\0'};
  return CubitString(tmp);
}


















