//- Class: Tree
//- Description: Tree is a tree class that accepts generic
//-              pointers to any object.  It is controlled by a macro that
//-              substitutes the specific pointer definitions in the class
//-              declaration. 
//+
//-              The macro {Treedeclare(name,typePtr)} is defined to
//-              create specific instances of the Tree class.
//-              {name} is the name of the Tree class created, and
//-              {typePtr} is the type of elements stored in the list.
//+
//-              The children interface is like a singlely linked list
//-              terminated by a NULL, but
//-              the underlying implementation is based on DLList
//- Assumptions: 
//-
//- Owner: Scott Mitchell
//- Checked by: 
//- Version: $Id: 

#ifndef TREE_HPP
#define TREE_HPP

#include "DLIList.hpp"
#include "MemoryManager.hpp"
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT Tree 
{
private:
  static MemoryManager memoryManager;
  //- memory management object

  DLIList<Tree*> children;   
  int endChildFlag; 
  //- have we traversed *past* the end of the children list

  Tree *parentPtr;
  
protected:
  void *data;  
  //- the node data

public:

  Tree( void *data_set = NULL, Tree *parent_set = NULL);
  //- constructor

  virtual ~Tree( );
  //- delete this and all its descendents (recursive destructor)
  
  virtual void delete_tree(){};
  //- like the destructor, only data members are recursively deleted, too.

  int is_leaf();  
  //- TRUE if no children

  int num_children();
  //- Number of child tree nodes

#ifdef BOYD15
  int depth();
  //- path length from this to root (i.e. root is at depth 0);
#endif

  void parent_item( Tree *parent_set );
  Tree *parent_item();
  //- get/set parent  

  Tree *root_item();
  //- root of whole tree of which this is a part

  //- Heading: functions for stepping through the list of children
  //- Treedeclare has these functions, minus the "_item"
  //- Try :
  //+ for (first_child(); current_child(); next_child() ) 
  //+   { current_child()->fn; }
  Tree *first_child_item();
  Tree *current_child_item();
  Tree *next_child_item();

  void empty_children();
  
  //- Heading: functions for depth first search (DFS) traversal.
  Tree *DFS_start_item();
  //- start a depth first search with this as root,
  //- DFS_next goes above this level when done.

  Tree *DFS_next_item();
  //- find next element of DFS given that this is the current element.
  //- Return NULL when done. Based on a get_and_step, so any 
  //- "current" function would be one step ahead.

  Tree *DFS_next_leaf_item();
  //- Visit the leaves in DFS order: 
  //- Keep looking at the next element until the search is done or
  //- a leaf is found.

  void add_child_item ( Tree *child );
  //- add child to the list of this's children, set child's parent to this

  SetDynamicMemoryAllocation(memoryManager)
  //- class specific new and delete operators, using memoryManager

};  

inline 
Tree::Tree( void *data_set, Tree *parent_set) {
  data = data_set;
  parentPtr = parent_set;  
  children.clean_out();
  endChildFlag = CUBIT_FALSE;
}

inline
Tree *Tree::parent_item( ) { return parentPtr; }

inline
void Tree::parent_item( Tree *parent_set ) { parentPtr = parent_set; }

inline
int Tree::is_leaf() { return children.size() == 0; }

inline
Tree *Tree::current_child_item() {
  return (children.size() && !endChildFlag) ? children.get() : NULL;
}

inline
Tree *Tree::first_child_item() {
  endChildFlag = CUBIT_FALSE;
  children.reset(); 
  return current_child_item();
}

inline
Tree *Tree::next_child_item() {
  endChildFlag = children.is_at_end();
  return (children.size() && !endChildFlag) ? children.step_and_get() : NULL;
}

inline
Tree *Tree::DFS_start_item() {
  first_child_item();
  return this;
}

inline void Tree::empty_children()
{
  children.clean_out();
}


inline
void Tree::add_child_item ( Tree *child ) {
  children.append( child );
  child->parent_item( this );
}

//- Heading : Treedeclare = how to make a tree with a specific type of node data
//+ e.g. Treedeclare(CubitFaceTree, CubitFace*); 
//+ // makes a CubitFaceTree class, with CubitFace node data.

#define Treedeclare(name,typePtr)                                          \
   class name : public Tree                                                \
{                                                                          \
public:                                                                    \
  name ( typePtr data_set = NULL, name *parent_set = NULL) :               \
    Tree ( data_set, parent_set ) {}                                       \
  name *parent() { return (name *) parent_item(); }                        \
  void parent( name *parent_set ) { parent_item( parent_set ); }           \
  name *root() { return (name *) root_item(); }                            \
  name *first_child() {  return (name *) first_child_item(); }             \
  name *next_child() { return (name *) next_child_item(); }                \
  name *child() { return (name *) current_child_item(); }                  \
  name *DFS_start() { return (name *) DFS_start_item(); }                  \
  name *DFS_next() { return (name *) DFS_next_item(); }                    \
  name *DFS_next_leaf() { return (name *) DFS_next_leaf_item(); }          \
  void add_child ( name *child ) { add_child_item( child ); }              \
  typePtr node_data() { return (typePtr) data; }                           \
  void node_data( typePtr data_set ) { data = data_set; }                  \
  void delete_tree()\
       {  for( first_child_item(); current_child_item(); next_child_item() )\
    current_child_item()->delete_tree();\
  if ( data ) delete (name *) data;\
  empty_children();\
  delete this;}\
}

#endif //- TREE_HPP

