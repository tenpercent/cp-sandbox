//- Class:       CubitEntity 
//- Description: CubitEntity class - the base class in the CUBIT Entity Tree.
//- Owner:       Bill Bohnhoff
//- Checked by:  Tim Tautges, 6/6/94
//- Version: $Id: 

#ifndef CUBITENTITY_HPP
#define CUBITENTITY_HPP

#include "CubitDefines.h"
#include "InvalidEntity.hpp"

#include <typeinfo>
#if !defined(NT)
using std::type_info;
#endif


class CubitBox;
class CubitVector;
class CubitString;
#ifdef CUBIT_GUI
class IGUIObservers;
#endif
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT CubitEntity
{
public:
  
    //- Heading: Constructors and Destructor
  CubitEntity() : entityId(0)
#ifdef CUBIT_GUI
                , pGUIObservers(0)
#endif
    {}
  
  virtual ~CubitEntity() ;
  
    //- Heading: Set and Inquire functions
  virtual int  id() const;
  
  virtual void set_id(int i);
    //- set the id of this entity to i
  
    //- Heading: Virtual functions
    //virtual void list()                 const = 0;           //- pure virtual
  
  virtual CubitBox bounding_box() = 0;
  virtual CubitVector center_point();
  
    //@ The overloaded members color() and is_visible() are
    //@ not pure virtual because mesh entities will NOT override these.
    //@ Mesh entities rely on the RefEntity owner for this information.
  virtual void color(int value);
  virtual int  color()                const;
  
  virtual void is_visible(int flag);
  virtual int is_visible()           const;
  virtual void is_transparent(int flag);
  virtual int is_transparent()     const;
  
  virtual const type_info&  entity_type_info()   const = 0;
  
  virtual const char* class_name() const = 0;
    //- return class name string.
  
  virtual int validate();
    //R int
    //R- number of problems detected, 0 if none (or not implemented)
  
#ifdef CUBIT_GUI
    // functions to get and set GUI event notification interface pointer
  IGUIObservers* get_GUI_observers_ptr();
  void set_GUI_observers_ptr(IGUIObservers* observers);
#endif
  
  
protected:
  int entityId;
  
#ifdef CUBIT_GUI
  IGUIObservers* pGUIObservers;
#endif
  
private:
  CubitEntity( const CubitEntity& );
  void operator=( const CubitEntity&);
};

#endif

