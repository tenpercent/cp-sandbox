//- File:           CAEntityColor.cpp
//- Owner:          Joel Kopp
//- Description:    Cubit Attribute for entity colors.
//- Checked By:
//- Version:

#include "CAEntityColor.hpp"
#include "RefEntity.hpp"
#include "RefVolume.hpp"
#include "CastTo.hpp"

CubitAttrib* CAEntityColor_creator(RefEntity* entity, CubitSimpleAttrib *p_csa)
{
  CAEntityColor *new_attrib = NULL;
  if (NULL == p_csa)
  {
    new_attrib = new CAEntityColor(entity);
  }
  else
  {
    new_attrib = new CAEntityColor(entity, p_csa);
  }

  return new_attrib;
}

CAEntityColor::CAEntityColor(RefEntity* new_attrib_owner,
                       CubitSimpleAttrib *csa_ptr)
        : CubitAttrib(new_attrib_owner)
{
   assert ( csa_ptr != 0 );
   
   DLIList<int*> *i_list = csa_ptr->int_data_list();
   
   assert(i_list && i_list->size() > 0);
   i_list->reset();
   entityColor = *( i_list->get_and_step() );
}

CAEntityColor::CAEntityColor(RefEntity* new_attrib_owner)
        : CubitAttrib(new_attrib_owner)
{
  entityColor = 0;
}

CAEntityColor::~CAEntityColor()
{
}


CubitStatus CAEntityColor::actuate()
{
   if ( hasActuated)
      return CUBIT_SUCCESS;
   
   if ( !attribOwnerEntity )
      return CUBIT_FAILURE;
   
   deleteAttrib = CUBIT_FALSE;
   
   //- If actuating after import, change the color.  Else (auto actuated) color is already changed...
   if( attribOwnerEntity->color() == CUBIT_DEFAULT_COLOR || 
	   dynamic_cast<RefVolume*>(attribOwnerEntity) )
   {
      attribOwnerEntity->color(entityColor);
   }
   
   hasActuated = CUBIT_TRUE;
   return CUBIT_SUCCESS;
}

CubitStatus CAEntityColor::update()
{
      delete_attrib(CUBIT_TRUE);
      return CUBIT_SUCCESS;

/*
   if ( hasUpdated ) 
      return CUBIT_SUCCESS;
   
   if ( !attribOwnerEntity)
      return CUBIT_FAILURE;
   
   entityColor = attribOwnerEntity->color();
   
   // set the updated flag
   hasUpdated = CUBIT_TRUE;

   if( entityColor == -1 )
   {
      delete_attrib(CUBIT_TRUE);
      return CUBIT_SUCCESS;
   }
   
   return CUBIT_SUCCESS;
   */
}

void CAEntityColor::merge_owner(CubitAttrib *deletable_attrib)
{
    // take the id with the lowest value
  CAEntityColor *other_caecolor = CAST_TO(deletable_attrib, CAEntityColor);
  
  if (other_caecolor)
    entityColor = other_caecolor->color();
}

CubitSimpleAttrib* CAEntityColor::cubit_simple_attrib()
{
  DLIList<CubitString*> cs_list;
  DLIList<int> i_list;

  i_list.append ( entityColor );
    
  cs_list.append(new CubitString(att_internal_name()));

  CubitSimpleAttrib* csattrib_ptr = new CubitSimpleAttrib(&cs_list,
                                                          NULL,
                                                          &i_list);
  int i;
  for ( i = cs_list.size(); i--;) delete cs_list.get_and_step();
  
  return csattrib_ptr;
}

void CAEntityColor::print()
{
    // print info on this attribute
  
  PRINT_INFO("CAEntityColor: owner = %s %d:   color ",
             attribOwnerEntity->class_name(), attribOwnerEntity->id());
  if ( entityColor > -1 )
     PRINT_INFO("                             %d\n",
     entityColor);
  else
     PRINT_INFO("\n");
}
