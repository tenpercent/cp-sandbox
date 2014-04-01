//-------------------------------------------------------------------------
// Filename      : RefCollection.cc
//
// Purpose       : Implementation of the RefCollection class.
//                 	   
//
// Special Notes : Created to support Pro/E assemblies, parts and features
//
// Creator       : Steve Storm
//
// Creation Date : 07/16/97
//
// Owner         : Steve Storm
//-------------------------------------------------------------------------

#include "RefCollection.hpp"
#include "CollectionEntity.hpp"

RefCollection::RefCollection()
    : collectionEntityPtr(NULL)
{}

RefCollection::~RefCollection()
{}


//-------------------------------------------------------------------------
// Purpose       : Sets the associated CollectionEntity
//
// Special Notes : 
//
// Creator       : Steve Storm
//
// Creation Date : 07/31/96
//-------------------------------------------------------------------------
CubitStatus RefCollection::set_collection_entity_ptr(
                           CollectionEntity* collection_entity_ptr ) 
{
   // Make sure that a CUBIT OWNER attribute is attached to the solid
   // model entity underlying the CollectionEntity.  This attribute contains 
   // a (back)pointer to this VGI GroupingEntity.

   CubitStatus status = collection_entity_ptr->set_owner_attribute( this );

   if (status == CUBIT_SUCCESS)
   {
      collectionEntityPtr = collection_entity_ptr ;
   }
   else
   {
      PRINT_ERROR("In RefCollection::set_collection_entity_ptr\n"
                  "       Unable to set owner attribute\n" );
   }
   
   return status ;
}
