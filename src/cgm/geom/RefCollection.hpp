#ifndef REF_COLLECTION_HPP
#define REF_COLLECTION_HPP

#include "RefEntity.hpp"

class CollectionEntity;

class CUBIT_GEOM_EXPORT RefCollection : public RefEntity
{
   public:
      RefCollection();
      virtual ~RefCollection();
      CollectionEntity* get_collection_entity_ptr(){ return collectionEntityPtr; }
      CollectionEntity const* get_collection_entity_ptr() const { return collectionEntityPtr; }
      //R CollectionEntity*
      //R- A pointer to a colleciton entity associated with the grouping entity.
      //- This function returns a pointer to the collection entity associated 
      //- with this object.
      
      CubitStatus set_collection_entity_ptr(CollectionEntity* collection_entity_ptr);
      
  protected:
     CollectionEntity* collectionEntityPtr;
     
};

#endif

