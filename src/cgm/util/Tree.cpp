//- File Tree.cc
//- Generic tree class

#include "Tree.hpp"


Tree *Tree::root_item() {
  Tree *root_ptr = this; 
  while ( root_ptr->parent_item() != NULL )
    root_ptr = root_ptr->parent_item(); 
  return root_ptr;
  //slower, simpler, recursive =  
  // return ( parent_item() ) ? parent_item()->root() : this;
}
Tree *Tree::DFS_next_item() {  
  Tree *a_child = current_child_item();
  // get_and_step
  if (a_child) {
    next_child_item();
    return a_child->DFS_start_item();
  }
  else
    return (parent_item()) ? parent_item()->DFS_next_item() : NULL;
}

Tree *Tree::DFS_next_leaf_item() {
  Tree *next = this;
  int keep_searching;
  do {
    next = next->DFS_next_item();
    keep_searching = ( next == NULL ) ? CUBIT_FALSE 
      : (next->is_leaf() == CUBIT_FALSE); // !NULL and !leaf
  } while (keep_searching);
  return next;
}

Tree::~Tree()
{
  for( first_child_item(); current_child_item(); next_child_item() )
    delete current_child_item();
}


