//-------------------------------------------------------------------------
// Filename      : TopologyBridge.hpp
//
// Purpose       : This file contains the declarations of the base class 
//                 TopologyBridge.  Sub-classes of TopologyBridge
//                 represent the link between the Cubit representation
//                 of an entity (such as RefVertex) and the GME representation
//                 of an entity (such as an acis VERTEX).
//
// Creator       : Darryl Melander
//
// Creation Date : 01/10/99
//
// Owner         : Darryl Melander
//-------------------------------------------------------------------------

#ifndef MODEL_ENTITY_BRIDGE_HPP
#define MODEL_ENTITY_BRIDGE_HPP


#include <typeinfo>
#if !defined(NT)
using std::type_info;
#endif

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitDefines.h"
#include "CubitGeomConfigure.h"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class TopologyEntity;
class RefVolume;
class GeometryQueryEngine;
class CubitSimpleAttrib;
class CubitSimpleAttrib;
class CubitString;
class BridgeManager;
class TBOwner;
template <class X> class DLIList;
class Lump;
class ShellSM;
class Surface;
class LoopSM;
class Curve;
class CoEdgeSM;
class Point;
class TopologyBridge;
class TopologyEntity;
class BodySM;
class Lump;
class ShellSM;
class Surface;
class LoopSM;
class Curve;
class CoEdgeSM;
class Point;

// ********** END FORWARD DECLARATIONS     **********

class CUBIT_GEOM_EXPORT TopologyBridge 
{
public:
  TopologyBridge() 
      : bridgeOwner(0),
        bridgeSense(CUBIT_FORWARD)
    {}

  virtual ~TopologyBridge();
  
  const type_info& generic_entity_type_info();
  
  virtual const type_info& topology_entity_type_info() const = 0;
  
  virtual void append_simple_attribute_virt(CubitSimpleAttrib*) = 0;
    //R void
    //I name
    //I- A reference to a constant CubitString object which is the name
    //I- that is to be appended to the TopologyBridge object.
    //- The purpose of this function is to append a string
    //- attribute to the TB. The name is attached to each 
    //- of the underlying solid model entities.
  
  virtual void remove_simple_attribute_virt(CubitSimpleAttrib*) = 0;
    //R void
    //I CubitSimpleAttrib*
    //I- A reference to a CubitSimpleAttrib object which is the object
    //I- that is to be removed to this TB object.
    //- The purpose of this function is to remove a simple
    //- attribute from the TB. The attribute is attached to each of the
    //- underlying solid model entities this one points to.
  
  virtual void remove_all_simple_attribute_virt() = 0;
    //R void
    //I-
    //- The purpose of this function is to remove all simple
    //- attributes from the TB. 
  
  virtual CubitStatus get_simple_attribute(DLIList<CubitSimpleAttrib*>&) = 0;
    //R CubitSimpleAttrib*
    //R- the returned cubit simple attribute.
    //- The purpose of this function is to get the attributes
    //- of the geometry entity. The name is attached to the underlying solid
    //- model entity(ies) this one points to.
    //- The default implementation returns a NULL.
    //- MJP Note:
    //- This is the code that implements the requirement that names
    //- of VGI Entities propagate across solid model boolean
    //- operations.  The success of this relies, of course, on the underlying
    //- solid modeler being able to propagate attributes across
    //- such operations on its entities. If it cannot, then "names"
    //- of VGI entities will not propagate.
    
  virtual CubitStatus get_simple_attribute(const CubitString& name,
                                    DLIList<CubitSimpleAttrib*>& ) = 0;

  TopologyEntity* topology_entity() const;
    //R TopologyEntity*
    //R- Pointer to the TopologyEntity using this TopologyBridge.
    //- This (pure virtual) function returns a pointer to the
    //- TopologyEntity that is using this object, i.e. owner of this
    //- TopologyBridge. It is implied that a TopologyBridge is not
    //- used by more than one TopologyEntity.

  TBOwner* owner() const
    { return bridgeOwner; }
  
  void owner( TBOwner* new_owner )
    { bridgeOwner = new_owner; }
    
  BridgeManager* bridge_manager() const;
  void bridge_manager(BridgeManager* manager);
  
  CubitSense bridge_sense() 
    { return bridgeSense; }
  void reverse_bridge_sense() 
    { bridgeSense = bridgeSense == CUBIT_UNKNOWN ? CUBIT_UNKNOWN  :
                    bridgeSense == CUBIT_FORWARD ? CUBIT_REVERSED :
                                                   CUBIT_FORWARD; }
    //- Get or switch the bridge sense.
    //-
    //- The bridge sense is the sense of the TopologyBridge with
    //- respect to it's owning TopologyEntity.  The results of
    //- geometric evaluation functions on the ownging TopologyEntity
    //- may be affected by this value.  
    //-
    //- Currently, this is only set for Surfaces and Curves, as the
    //- sense of either may be opposite that of its owning RefFace
    //- or RefEdge when merging occurs.  It is not set/used on Loops
    //- or CoEdges because a) there are no geometric evaluations on
    //- the underlying LoopSMs or CoEdgeSMs and b) the relative sense
    //- can be inferred from that of the corresponding Surface and
    //- Curve.
    //-
    //- This value is saved in the MergePartner attribute.
    
  virtual GeometryQueryEngine* get_geometry_query_engine() const = 0;
    //- This function returns a pointer to the geometry query engine
    //- associated with the object.
  
  virtual int validate(const CubitString&,
                       DLIList <TopologyEntity*>&)
    { return 0; }
    //- Check that entity is valid. Returns number of problems detected.
    //- For sub-classes that don't implement this, just assume there are
    //- no problems and return 0.

  
  enum { MAX_TB_LAYER = 255 };
  virtual int layer() const { return 0; }
    //- Bridge layer at which this bridge occurs.
    //- 0   - real geometry
    //- 126 - sub-composite partition layer
    //- 127 - composite layer
    //- 128 - super-composite partition layer

  void get_parents(DLIList<TopologyBridge*> &parents);
    //- get parent topology bridges

  void get_children(DLIList<TopologyBridge*> &children,
                    bool return_hidden_entities = false,
                    int layer = MAX_TB_LAYER );
    //- get child topology bridges
  
  void get_related(const type_info& other_type,
                   DLIList<TopologyBridge*> &related);
    //- get the related entities of type other_type
  
#ifdef BOYD14
  CubitBoolean is_related(TopologyBridge *other_entity);
    //- returns CUBIT_TRUE if this and other_entity are directly related
#endif
  void bodysms(DLIList<BodySM*> &bodies,bool unique = true);
  void lumps(DLIList<Lump*> &lumps,bool unique = true);
  void shellsms(DLIList<ShellSM*> &shellsms,bool unique = true);
  void surfaces(DLIList<Surface*> &surfaces,bool unique = true);
  void loopsms(DLIList<LoopSM*> &loopsms,bool unique = true);
  void curves(DLIList<Curve*> &curves,bool unique = true);
  void coedgesms(DLIList<CoEdgeSM*> &coedgesms,bool unique = true);
  void points(DLIList<Point*> &points,bool unique = true);
    //- topology traversal of TB's; implemented based on native traversal
    //- functions in the modeler
  
  BodySM *bodysm();
  Lump *lump();
#ifdef BOYD14
  ShellSM *shellsm();
  Surface *vgi_surface();
  Curve *vgi_curve();
  CoEdgeSM *coedgesm();
  Point *point();
    //- topology traversal of TB's; implemented based on native traversal
    //- functions in the modeler
#endif
  LoopSM *loopsm();
  
  virtual void get_parents_virt(DLIList<TopologyBridge*> &parents ) = 0;
    //- Get parent topology bridges as returned by solid modeler.
    //- Derived classes must provide this method for use in the
    //- implementation of get_parents(..).  You probably want to 
    //- be using get_parents(..) rather than this method.

  virtual void get_children_virt(DLIList<TopologyBridge*> &children ) = 0;
    //- Get child topology bridges as returned by solid modeler.
    //- Derived classes must provide this method for use in the
    //- implementation of get_children(..).  You probably want to 
    //- be using get_children(..) rather than this method.
  
private:

  
  TBOwner* bridgeOwner;
  
  CubitSense bridgeSense;
    //- (See comments for bridge_sense() method.)
  
};

#endif

