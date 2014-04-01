//- Class:          CAMergePartner
//- Owner:          Tim Tautges
//- Description:    Cubit attribute for RefGroups
//- Checked by:
//- Version:

#ifndef CA_GROUP_HPP
#define CA_GROUP_HPP

#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CADefines.hpp"

class RefEntity;
class RefGroup;

class CUBIT_GEOM_EXPORT CAGroup: public CubitAttrib
{

private:

  DLIList<int> groupID;
    //- group ids containing attribOwnerEntity
  
  DLIList<int> uniqueID;
    //- unique ids of groups containing attribOwnerEntity

  DLIList<CubitString*> groupNames;
    //- names of groups containing attribOwnerEntity

  DLIList<int> sequenceNumbers;
    //- sequence numbers of this entity in the groups

  DLIList<int> numOwningGroups;
    //- for each group in groupID, number of groups owning those groups

  DLIList<int> owningGroupID;
    //- group ids containing groups containing attribOwnerEntity
  
  DLIList<int> owningUniqueID;
    //- unique ids of groups containing groups containing attribOwnerEntity

  DLIList<CubitString*> owningGroupNames;
    //- names of groups containing groups containing attribOwnerEntity

  DLIList<int> owningSequenceNumbers;
    //- sequence numbers of groups in owning groups
  
    //- for each ancestor (a group which owns only other groups, with those
    //- those groups owning only other groups), we store the group id, uid,
    //- name, and the uid of the owned group to which this is an ancestor
  DLIList<int> ancestorGroupID;
  DLIList<int> ancestorUniqueID;
  DLIList<CubitString*> ancestorGroupName;
  DLIList<int> ancestorOwnedGroupUid;
  DLIList<int> ancestorSequenceNumbers;
  
  static CubitBoolean initialize_rand;
  
public:
  CAGroup(RefEntity* = NULL);

  CAGroup(RefEntity*, CubitSimpleAttrib*);
    //- make a CAG from a simple attribute

  void initialize();
    //- initialize random number generator for this attribute

  virtual ~CAGroup();

  //HEADER- RTTI and safe casting functions.
  virtual const type_info& entity_type_info() const
     { return typeid(CAGroup);}
  //R- The geometric modeler type
  //- This function returns the type of the geometric modeler.

  CubitStatus actuate();
    //- actuate this attribute
  
  CubitStatus update();
    //- update this attribute

  CubitStatus reset();
    //- reset this attribute

  CubitSimpleAttrib* cubit_simple_attrib();
    //- return a simple attribute with this CA's data
  
  RefGroup *assign_group(RefEntity *owned_entity,
                         const int group_id, const int unique_id,
                         CubitString *group_name,
                         const int seq_num);
    //- used in actuating this CA
  
  RefGroup *assign_ancestor_group(const int ancestor_id,
                                  const int ancestor_uid,
                                  const CubitString *ancestor_name,
                                  const int owned_group_uid,
                                  const int seq_num);
  
  void build_ancestor_list(RefGroup *parent_ref_group);
  
  int int_attrib_type() {return CA_GROUP;}
    //- returns the enumerated attribute type

  virtual void has_written(CubitBoolean set_has_written);
    //- class-specific version; resets td_cage on attrib owner

  virtual CubitBoolean has_written() const {return hasWritten;};
    //- class-specific version; 

};

CubitAttrib* CAGroup_creator(RefEntity* entity, CubitSimpleAttrib *p_csa);

#endif

