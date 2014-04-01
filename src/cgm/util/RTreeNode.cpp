//---------------------------------------------------------------------------
// Class Name:  RTreeNode
// Description: Node of Rectangle tree.  Contians many of the
//              required functions for building the tree and traversing it.
// The algorithm was taken from the following paper:
//	      Guttman, A., "R-Trees: A Dynamic Index Structure for
//              Spatial Searching", Proceedings of the SIGMOD
//              Conference, Boston, June 1984, p. 47-57.
// Creation Date: 3/13/02
// Owner:  David R. White
//---------------------------------------------------------------------------

//---------------------------------
//Include Files
//---------------------------------
#include "RTreeNode.hpp"
#include "DLIList.hpp"
//---------------------------
//Initialize Static Members
//---------------------------

#ifdef INLINE_TEMPLATES
#define MY_INLINE inline
#else
#define MY_INLINE
#endif


template <class Y> MY_INLINE RTreeNode<Y>::RTreeNode (Y data, double tol,
                                                      int max_children,
                                                      int min_children)
{
  maxChildren = max_children;
  minChildren = min_children;
  myChildrenNodes = new RTreeNode<Y>* [maxChildren];
  int ii;
  for ( ii = 0; ii < maxChildren; ii++ )
    myChildrenNodes[ii] = static_cast< RTreeNode<Y>* >(NULL);
  if ( data == NULL )
  {
    PRINT_ERROR("Building RTree with null data is not allowed!\n");
    assert(data != NULL);
  }
  myData = data;
  myLevel = DATA_RNODE;
  CubitBox temp_box = data->bounding_box();
    //Check to see if any of the min/max pairs are less than the tolerance.
    //make them bigger if they are...
  CubitVector min = temp_box.minimum();
  CubitVector max = temp_box.maximum();
  if ( max.x() - min.x() < tol )
  {
    min.x(min.x()-.6*tol);
    max.x(max.x()+.6*tol);
  }
  if ( max.y() - min.y() < tol )
  {
    min.y(min.y()-.6*tol);
    max.y(max.y()+.6*tol);
  }
  if ( max.z() - min.z() < tol )
  {
    min.z(min.z()-.6*tol);
    max.z(max.z()+.6*tol);
  }
  myBoundingBox = new CubitBox(min, max);
  myParent = NULL;
  nextChildIndex = 0;
  markedFlag = 0;
  distIsBox = 1;
  myDist = CUBIT_DBL_MAX;
}
template <class Y> MY_INLINE RTreeNode<Y>::RTreeNode (CubitBox &bound_box,
                                                      int max_children,
                                                      int min_children)
{  
  maxChildren = max_children;
  minChildren = min_children;
  myBoundingBox = new CubitBox(bound_box);
  myChildrenNodes = new RTreeNode<Y>* [maxChildren];
  int ii;
  for ( ii = 0; ii < maxChildren; ii++ )
    myChildrenNodes[ii] = static_cast< RTreeNode<Y>* >(NULL);
  myData = NULL;
  myLevel = UNSET_RNODE;
  myParent = NULL;
  nextChildIndex = 0;
  markedFlag = 0;
  distIsBox = 1;
  myDist = CUBIT_DBL_MAX;
}
//-----------------------------------------------------------
// Destructor
//-----------------------------------------------------------
template <class Y> MY_INLINE RTreeNode<Y>::~RTreeNode()
{
  if ( myChildrenNodes )
    delete [] myChildrenNodes;
  if ( myBoundingBox )
    delete myBoundingBox;
}

//-----------------------------------------------------------
// Algorithm: insert
//  Insert a new index entry e into an R-Tree.
//     I1. [Find postiion for new record] Invoke choose_leaf to select
//        a leaf node in l, in which to place e.
//     I2. [Add record to leaf node].a) If l has room for
//        another entry, install E. b) Otherwise invoke split_node to
//        obtain l and ll containing e and all the old entries of l.
//     I3. [Propogate changes upward] Invoke adjust_tree on l, also passing ll
//        if a split was performed.
//     I4. [Grow Tree Taller] If node split propogation caused the root
//         to split create a new root whose children are the two resulting
//         nodes.
//-----------------------------------------------------------
template <class Y> MY_INLINE CubitStatus RTreeNode<Y>::insert(RTreeNode<Y> *e, RTreeNode<Y> *&new_root)
{
  CubitStatus stat;
  new_root = NULL;//only set this if the root node changes.  Assume
    //that this RTreeNode object is the root...

  
    // I1. Invoke choose_leaf to select a leaf node l in which to place
    //e
  RTreeNode<Y> *l = choose_leaf(this, e);
    //just test.
    // make sure l is not null.
    // make sure l is one level above e.
  if ( l ==  NULL || l->get_leaf_level() != (e->get_leaf_level() + 1) )
  {
    PRINT_ERROR("Choosing leaf for inseartion into rtree failed.\n");
    return CUBIT_FAILURE;
  }
  RTreeNode<Y> *ll = NULL;
    //I2 a) If l has room for another entry install e.
  if ( l->can_add() )
  {
    l->add_child(e, CUBIT_TRUE);
  }
  else
  {
      //I2 b)
      //Otherwise invoke split_node to obtain l and ll containing e and
      //all the old entries of l.
    stat = quadratic_split_node(l,e,ll);
    if ( stat != CUBIT_SUCCESS || ll == NULL )
    {
      PRINT_ERROR("Problems splitting node during insertion to RTree.\n");
      return CUBIT_FAILURE;
    }
  }
    //I3.  Propagate changes upward.
    //I4, grow tree taller (do both inside adjust_tree...)
  stat = adjust_tree(l, ll, this, new_root);
  if ( stat!= CUBIT_SUCCESS )
    return stat;
  return CUBIT_SUCCESS;
}

//--------------------------------------------
// Algorithm: ChooseLeaf: Select a leaf node in which to place
// a new index entry e.  Recursive search the subtrees of n
// until n is a leaf node.
//----------------------------------------------
template <class Y> MY_INLINE RTreeNode<Y>* RTreeNode<Y>::choose_leaf( RTreeNode<Y>* n,
                                                            RTreeNode<Y>* e )
{   
    //If n is a leaf node, or one level greater than e,
    //we are done.
  if ( n->get_leaf_level() == (e->get_leaf_level() + 1) )
    return n;
    //Now choose the entry f in n (children of n that is)
    //whose bounding box f_box needs
    //least enlargement to include e_box. Resolve ties by
    //choosing the entry with the bounding box of smallest volume.
  double min_enlargement = CUBIT_DBL_MAX, curr_enlargement;
  RTreeNode<Y> *curr_node;
  int child_index = -1;
  int ii;
  for(ii = 0; (ii < maxChildren) && (n->myChildrenNodes[ii] != NULL); ii++  )
  {

    curr_node = n->myChildrenNodes[ii];
    curr_enlargement = calc_enlargement(curr_node, e);
    if ( curr_enlargement <= min_enlargement )
    {
      if ( curr_enlargement == min_enlargement && child_index >= 0 )
      {
          //only reset if the curr_node has a smaller volume.
        double curr_vol = volume(curr_node);
        double old_vol = volume(n->myChildrenNodes[child_index]);
        if ( old_vol < curr_vol )
          continue;
      }
      child_index = ii;
      min_enlargement = curr_enlargement;
    }  
  }
    //do error checking...
  if ( child_index == -1 || child_index >= maxChildren )
    return static_cast< RTreeNode<Y>* >(NULL);
  RTreeNode<Y> *f = n->myChildrenNodes[child_index];
    //Now continue on...
  curr_node = choose_leaf(f,e);
  return curr_node;
}
//----------------------------------------------------------------------
// calc_enlargement:  Calculate the enlargement required for increasing
// the bounding box of current so that it would encapsulate the bounding
// box of add_to.  So to do that, create the union of the two bounding
// boxes, then of that supper box subtrace the volume of the current.
// The result should be the volumetric difference between how much
// current has and how much it would need be or the enlargement.
//----------------------------------------------------------------------
template <class Y> MY_INLINE double RTreeNode<Y>::calc_enlargement(RTreeNode<Y> *current, RTreeNode<Y> *add_to )
{
    //The enlargement area is the volume of the box that would
    //be the union of current and add_to minus the volume of the current.
  CubitBox curr_box = current->bounding_box();
  CubitBox add_to_box = add_to->bounding_box();
  CubitBox supper = curr_box;
    //Unite add_to_box to the curr_box.
  supper|= add_to_box;
  double area_big = volume(supper);
  return area_big - volume(current);
}
template <class Y> MY_INLINE double RTreeNode<Y>::calc_enlargement(CubitBox &current, CubitBox &add_to )
{
    //The enlargement area is the volume of the box that would
    //be the union of current and add_to minus the volume of the current.
  CubitBox supper = current;
    // unite the add_to box.
  supper |= add_to;
  double area_big = volume(supper);
  return area_big - volume(current);
}
//------------------------------------------------------------------
// Algorithm: adjust_tree
// Description:  Ascend from a leaf node L to the root, adjusting covering
// bounding boxes and propagating nodes splits as necesary.
//------------------------------------------------------------------
template <class Y> MY_INLINE CubitStatus RTreeNode<Y>::adjust_tree(RTreeNode<Y> *l, RTreeNode<Y> *ll,
                                   RTreeNode<Y> *root_node,
                                   RTreeNode<Y> *&new_root)
{
  CubitStatus stat;
  //we need to move up the tree and correct things that have changed.
  if ( l == root_node )
  {
    if ( ll == NULL )
      return CUBIT_SUCCESS;
    else
    {
        //Create a new root node and store l and ll there
	  CubitBox root_box = l->bounding_box();
      root_box |= ll->bounding_box();
      new_root = new RTreeNode<Y>(root_box, maxChildren, minChildren);
      int new_level = l->get_leaf_level() + 1;
      new_root->set_leaf_level(new_level);
      new_root->add_child(l, CUBIT_TRUE);
      new_root->add_child(ll, CUBIT_TRUE);
      return CUBIT_SUCCESS;
    }
  }
  RTreeNode<Y> *parent_node = l->get_parent();
  RTreeNode<Y> *new_group = NULL;
  if ( ll != NULL )
  {
      //We need to add ll to the parent if we can,
      //and then we need to update the parent's bounding box...
    if ( parent_node->can_add() )
    {
      parent_node->add_child(ll, CUBIT_FALSE);
        //we need to recalculate the bounding box for the
        //entire set since both l and ll were modified...
      parent_node->recalc_b_box();
    }
    else
    {
        //Now we must split the children of the parent. l should
        //already be in the chilren list of the paretn.  So send
        //to split node the parent_node and ll.
        //parent node during this process will have its b_box recalced.
      stat = quadratic_split_node(parent_node, ll, new_group);
      if ( stat != CUBIT_SUCCESS || new_group == NULL )
      {
        PRINT_ERROR("Problems splitting node during insertion to RTree.\n");
        return CUBIT_FAILURE;
      }
    }
  }
  else
  {
      //just recalulate the b_box for the parent_node.
    parent_node->recalc_b_box();
  }
  stat = adjust_tree(parent_node, new_group, root_node, new_root);
  if ( stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
  return CUBIT_SUCCESS;
}
//------------------------------------------------------------------
// Algorithm: quadratic_split_node
//  Description:  For the RTree, this is where most of the variations
//  on implementations have occured.  The best method proposed by Guttman
//  was a quadratic split, which I'll implement here.  The Rstar tree
//  did some slightly more complicated things which I might try later.
// Description of function:
//  e is the node to which we want to add l, but l's children's list
//  is full.  Split the children of l and e up into two groups.  The
//  groups are stored in l and ll.
// Assume:  assume that l's myChildrenNodes are maxChildren in size.
//
//  Divide a set of maxChildren + 1 index entries into two groups.
//  QS1  [Pick first entry for each group]  Apply algorithm pick_seeds,
//       to choose two entries to be the first elements of the groups.
//       Assign each to a group.
//  QS2  [Check if done]  If all entries have been assigned, stop.
//       If one group has so few entries that all the rest must
//       be assigned to it in order for it to have the minimum number m,
//       assign them and stop.
//  QS3  [Select entry to assign]  Invoke Algorithm pick_next to choose
//       the next entry to assign.  Add it to the group whose covering
//       rectangle will have to be enlarged least to accommodate it.
//       Resolve ties by adding the entry to the group with the smaller
//       volume, then to the one with fewer entries, then to either. Repeat
//       QS2.
//
//------------------------------------------------------------------
template <class Y> MY_INLINE CubitStatus RTreeNode<Y>::quadratic_split_node( RTreeNode<Y> *l,
                                             RTreeNode<Y> *e,
                                             RTreeNode<Y> *&ll )
{
  int ii;
    //create a new list containing all the nodes we want to split.
  RTreeNode<Y> **input_list = new RTreeNode<Y>* [maxChildren + 1];
  DLIList <RTreeNode<Y>*> nodes_remaining;
  for ( ii = 0; ii < maxChildren; ii++ )
  {
    input_list[ii] = l->myChildrenNodes[ii];
  }
  input_list[maxChildren] = e;
  
    //QS1, pick first entry for each group.
  RTreeNode<Y> *seed_1, *seed_2;
  CubitStatus stat = pick_seeds(input_list, maxChildren+1,seed_1,
                                seed_2);
  if ( stat != CUBIT_SUCCESS )
  {
    delete [] input_list;
    return stat;
  }
    //now flush out l. This cleans out the bounding box and
    //chindrenNodes and resets the bounding box.  Also
    //create ll, make l and ll non-leaf nodes and add
    //seed_1 and seed_2 to l and ll. (doesn't matter which...)
  l->flush(seed_1->bounding_box());
  l->add_child(seed_1, CUBIT_FALSE);
    //this is still a leaf node.
    //this will change if necessary the parent...
  l->set_leaf_level(e->get_leaf_level() + 1);
  ll = new RTreeNode<Y>(seed_2->bounding_box(), maxChildren, minChildren );
  ll->add_child(seed_2, CUBIT_FALSE);
  ll->set_leaf_level(l->get_leaf_level());
    //build the nodes remaining list...
  for ( ii = 0; ii < maxChildren+1; ii++ )
  {
    if ( input_list[ii] != seed_1 &&
         input_list[ii] != seed_2 )
      nodes_remaining.append(input_list[ii]);
  }
  delete [] input_list;
  RTreeNode<Y> *next_node;
  CubitBoolean add_to_group_1;
    //Q2
  while (nodes_remaining.size() > 0 )
  {
      //Q2 continued.
    if ( l->space_left() < minChildren &&
         minChildren - l->space_left() >= nodes_remaining.size() )
    {
        //just add the rest of the nodes to l.
      for ( ii = nodes_remaining.size(); ii > 0; ii-- )
        l->add_child(nodes_remaining.get_and_step(), CUBIT_TRUE);
      nodes_remaining.clean_out();
      break;
    }
    else if ( ll->space_left() < minChildren &&
         minChildren - ll->space_left() >= nodes_remaining.size() )
    {
        //just add the rest of the nodes to l.
      for ( ii = nodes_remaining.size(); ii > 0; ii-- )
        ll->add_child(nodes_remaining.get_and_step(), CUBIT_TRUE);
      nodes_remaining.clean_out();
      break;
    }
      //Q3
      //pick next selects the next node and the group to
      //put it in.  It also removes the node from nodes remaining.
      //Some of these steps were added to pick next for efficiency...
    stat = pick_next(nodes_remaining, l, ll, next_node,
                     add_to_group_1 );
    if ( stat != CUBIT_SUCCESS )
      return stat;
    if ( add_to_group_1 )
      l->add_child(next_node, CUBIT_TRUE);
    else
      ll->add_child(next_node, CUBIT_TRUE);
  }
  return CUBIT_SUCCESS;
}
template <class Y> MY_INLINE void RTreeNode<Y>::flush( CubitBox &new_box )
{
  int ii;
  nextChildIndex = 0;
  for ( ii = 0; ii < maxChildren; ii++ )
    myChildrenNodes[ii] = NULL;
  delete myBoundingBox;
  myBoundingBox = new CubitBox(new_box);
}
template <class Y> MY_INLINE void RTreeNode<Y>::add_child(RTreeNode<Y> *child_node,
							  CubitBoolean recalc_bound_box)
{
  assert(nextChildIndex < maxChildren && child_node != NULL );
  myChildrenNodes[nextChildIndex] = child_node;
    //update the bounding box. by uniting with child node...
  if ( recalc_bound_box )
  {
    CubitBox *old_box = myBoundingBox;
    myBoundingBox = new CubitBox( *old_box |= child_node->bounding_box());
    delete old_box;
  }
  nextChildIndex++;
  child_node->set_parent(this);
}
template <class Y> MY_INLINE CubitBoolean RTreeNode<Y>::can_add()
{
  if (nextChildIndex >= maxChildren )
    return CUBIT_FALSE;
  else
    return CUBIT_TRUE;
}
template <class Y> MY_INLINE int RTreeNode<Y>::space_left()
{
  return maxChildren - nextChildIndex;
}
//------------------------------------------------------------------
// Algorithm pick_seeds: 
//  Select two entries to be the first elements of the groups.
//  PS1 [Calculate inefficiency of grouping entries together] For
//      each pair of entries, e1 and e2, compose a rectangle (bounding_box)
//      j, including e1 and e2, calculate:
//      d = volume(j) - volume(e1) - volume(e2)
//  PS2 [Choose the most wasteful pair] Choose the pair with the largest
//      d.
//------------------------------------------------------------------
template <class Y> MY_INLINE CubitStatus RTreeNode<Y>::pick_seeds(RTreeNode<Y> **input_list,
                                                        const int input_list_size,
                                                        RTreeNode<Y> *&seed_1,
                                                        RTreeNode<Y> *&seed_2)
{
  int ii, jj;
  RTreeNode<Y> *e_1, *e_2;
  CubitBox e_box_1, e_box_2, j;
  double d, max_d = -CUBIT_DBL_MAX;
  seed_1 = static_cast< RTreeNode<Y>* >(NULL);
  seed_2 = static_cast< RTreeNode<Y>* >(NULL);
  
  for(ii = 0; ii < input_list_size; ii++ )
  {
    e_1 = input_list[ii];
    e_box_1 = e_1->bounding_box();
    for ( jj = ii+1; jj < input_list_size; jj++ )
    {
      e_2 = input_list[jj];
      e_box_2 = e_2->bounding_box();
        //unite the boxes.
	  j = e_box_1;
      j |= e_box_2;
        //find the most wastefull boxes to separate the groups.
      d = volume(j) - volume(e_box_1) - volume(e_box_2);
      if ( d > max_d )
      {
        seed_1 = e_1;
        seed_2 = e_2;
        max_d = d;
      }
    }
  }
  return CUBIT_SUCCESS;
}
//------------------------------------------------------------------
// Algorithm pick_next: 
//  Select one remaining entry for classification in a group.
//  PN1 [Determine cost of putting each entry in group] For each entry
//      e not yet in a group, calculate d1 = area increase required in the
//      covering rectangle of group_1 to include E.  Calculate d2 similarly
//      for group_2
//  PN2 [Find entry with greatest preference for one group]  Choose
//      any entry with the maximum difference between d1 and d2.
//------------------------------------------------------------------
template <class Y> MY_INLINE CubitStatus RTreeNode<Y>::pick_next(DLIList <RTreeNode<Y>*> &remaining_nodes,
                                                       RTreeNode<Y>* group_1,
                                                       RTreeNode<Y>* group_2,
                                                       RTreeNode<Y>*& next,
                                                       CubitBoolean &add_to_group_1)
{
  int ii, next_index = 0;
  double d1, d2, max_diff = -CUBIT_DBL_MAX;
  add_to_group_1 = CUBIT_TRUE;
  RTreeNode<Y> *max_diff_node = static_cast< RTreeNode<Y>* >(NULL);
  RTreeNode<Y> *curr_node;
  CubitBox group_1_box = group_1->bounding_box();
  CubitBox group_2_box = group_2->bounding_box();
  CubitBox curr_box;
  double v1 = volume(group_1_box);
  double v2 = volume(group_2_box);
  remaining_nodes.reset();
  for ( ii = 0; ii < remaining_nodes.size(); ii++ )
  {
    curr_node = remaining_nodes.get_and_step();
    curr_box = curr_node->bounding_box();
    d1 = calc_enlargement(group_1_box, curr_box);
    d2 = calc_enlargement(group_2_box, curr_box);
    if ( d1 > d2 )
    {
      if ( max_diff < d1 - d2 )
      {
        //add to group whose covering area would have
        //to be enlarged least.
        add_to_group_1 = CUBIT_FALSE;
        max_diff = d1 - d2;
        max_diff_node = curr_node;
        next_index = ii;
      }
    }
    else if ( d2 > d1 )
    {
      if ( max_diff < d2 - d1 )
      {
        //add to group whose covering area would have
        //to be enlarged least.
        add_to_group_1 = CUBIT_TRUE;
        max_diff = d2 - d1;
        max_diff_node = curr_node;
        next_index = ii;
      }
    }
    else {
      if ( max_diff < 0.0 )
      {
          //Add to group with smaller area.
        if ( v1 < v2 )
          add_to_group_1 = CUBIT_TRUE;
        else if ( v2 < v1 )
          add_to_group_1 = CUBIT_FALSE;
        else
        {
            //add to group with fewest entries.
          int num_left_1 = group_1->space_left();
          int num_left_2 = group_2->space_left();
          if ( num_left_1 > num_left_2 )
            add_to_group_1 = CUBIT_TRUE;
          else
            add_to_group_1 = CUBIT_FALSE;
        }
        max_diff = 0.0;
        max_diff_node = curr_node;
        next_index = ii;
      }
    }
  }
  next = NULL;
  if ( max_diff_node == NULL )
    return CUBIT_FAILURE;
  else
    next = max_diff_node;
    //remove next from the remaining_nodes list.
  remaining_nodes.reset();
  remaining_nodes.step(next_index);
  RTreeNode<Y> *check_node = remaining_nodes.remove();
  if ( check_node != max_diff_node )
  {
    PRINT_ERROR("Error in pick next algorithm logic of rtree...");
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}
template <class Y> MY_INLINE void RTreeNode<Y>::recalc_b_box()
{
  if(myLevel == DATA_RNODE )
    return;
  int ii;
  CubitBox temp_box;
  CubitBoolean first_box = CUBIT_TRUE;
  for ( ii = 0; ii < nextChildIndex; ii++ )
  {
	if ( first_box )
	{
	  temp_box = myChildrenNodes[ii]->bounding_box();
	  first_box = CUBIT_FALSE;
	}
	else
      temp_box |= myChildrenNodes[ii]->bounding_box();
  }
  delete myBoundingBox;
  myBoundingBox = new CubitBox(temp_box);
  return;
}
//-------------------------------------------------------------
// Algorithm: remove.  Remove index record e from an R-tree.
//   D1)  [Find node containing record].  Invoke find_leaf to locate
//        the leaf node l containing e.  Stop if the record was not
//        found.
//   D2)  [Delete record.]  Remove e from l.
//   D3)  [Propagate changes.]  Invoke CondenseTree, passing L.
//   D4)  [Shorten tree.]  If the root node has only one child
//        after the tree has been adjusted, make the child the new
//        root.
//-------------------------------------------------------------
template <class Y> MY_INLINE CubitBoolean RTreeNode<Y>::remove( Y e,
                                                      RTreeNode<Y> *&new_root,
                                                      CubitBoolean &delete_root)
{
    //D1) Find node containting record.
  if(this->is_data()){
    if(this->get_data() == e){
      myData = NULL;
      myLevel = UNSET_RNODE;
      return CUBIT_TRUE;
    }
    else if(this->num_children() == 0){
      return CUBIT_FALSE;
    }
  }
  RTreeNode<Y> *l = NULL;
  CubitBox my_box = e->bounding_box();
  CubitStatus stat = find_leaf(e, my_box, this, l);
  if ( l == NULL || stat != CUBIT_SUCCESS )
    return CUBIT_FALSE;
    //Now l is the RTreeNode that holds the actual data (a DATA_RNODE)
    //not a leaf.  This was done for efficiency.
  RTreeNode<Y> *data_node = l;
  l = data_node->get_parent();
    //D2) [Delete record]  Remove e from l.
  
    //remove the data node from the children and delete
    //the node.
  l->remove_child(data_node);
  delete data_node;

    //D3) [Propogate Changes].
  stat = condense_tree(l, this, new_root);
    //D4) [Shorten the tree].
  RTreeNode<Y> *root = this;
  if ( new_root != NULL )
    root = new_root;
  if ( root->num_children() == 1 )
  {
    new_root = root->get_child(0);
    new_root->set_parent(static_cast< RTreeNode<Y>* >(NULL));
    delete_root = CUBIT_TRUE;
  }
  return CUBIT_TRUE;
}
template <class Y> MY_INLINE CubitStatus RTreeNode<Y>::find_leaf( Y e,
                                                        CubitBox &e_box,
                                                        RTreeNode<Y> *t,
                                                        RTreeNode<Y> *&l )
{
  int ii;
  l = NULL;
  int loop_size = t->num_children();
  RTreeNode<Y> *curr_node;
  if ( t->get_leaf_level() > LEAF_RNODE )
  {
    for ( ii = 0; ii < loop_size; ii++ )
    {
      curr_node = t->get_child(ii);
      if ( curr_node == NULL )
      {
        PRINT_ERROR("Problems finding boxes in range.\n");
        assert(curr_node != NULL);
        return CUBIT_FAILURE;
      }
      if ( e_box.overlap(GEOMETRY_RESABS, curr_node->bounding_box()) )
      {
          //okay now search through this now.
        find_leaf(e, e_box,curr_node,l);
        if ( l != NULL )
          return CUBIT_SUCCESS;
      }
    }
  }
  else if ( t->is_leaf() )
  {
      //search through the children for e.
    for ( ii = 0; ii < loop_size; ii++ )
    {
      curr_node = t->get_child(ii);
      if ( curr_node == NULL )
      {
        PRINT_ERROR("Problems finding boxes in range.\n");
        assert(curr_node != NULL);
        return CUBIT_FAILURE;
      }
      if ( curr_node->get_data() == e )
      {
        l = curr_node;
        return CUBIT_SUCCESS;
      }
    }
  }
  return CUBIT_SUCCESS;
}
template <class Y> MY_INLINE CubitBoolean RTreeNode<Y>::remove_child( RTreeNode<Y> *child )
{
    //first find which item this child is at.
  int ii;
  int index_child = -1;
  int loop_size = this->num_children();
  for ( ii = 0; ii < loop_size; ii++ )
  {
    if ( myChildrenNodes[ii] == child )
      index_child = ii;
  }
  if ( index_child == -1 )
    return CUBIT_FALSE;
    //Now we need to bubble the array from this point
    //upward.
  for ( ii = index_child; ii < loop_size-1; ii++ )
  {
    myChildrenNodes[ii] = myChildrenNodes[ii+1];
  }
    //decrement the position of the next available child.
  nextChildIndex--;
    //now go from nextChildIndex to the end and make sure it is
    //null.
  for (ii = nextChildIndex; ii < maxChildren; ii++ )
    myChildrenNodes[ii] = NULL;
  
  return CUBIT_TRUE;
}
//--------------------------------------------------------------------------
// Algorithm: condense_tree
//  Given a leaf node l from which an entry has been deleted, eliminate
//  the node if it has too few entries and relocate its entries.  Propagate
//  node elimination upaward as necessary.  Adjust all covering rectangles
//  on the path to the root, making them smaller if possible.
//  CT1) [Initialize]  Set n=l, Set q, the set of eliminated nodes, to be
//       empty.
//  CT2) [Find parent entry]  If n is the root, go to CT6.  Otherwise
//       let p be the parent of n, and let en be n's entry in p.
//  CT3) [Eliminate under-full node.] If n has fewer than minChildren,
//       delete en from p and add n to set q.
//  CT4) [Adjust covering rectangle]  If n has not been eliminated, adjust
//       en's bounding box to tightly contain all entries in n.
//  CT5) [Move up one level in tree] Set n=p, and repeat from CT2.
//  CT6) [Reinsert orphaned entries].  Reinsert all entries of nodes in set q.
//       entries from eliminated leaf nodes are re-inserted in tree leaves
//       as described in algorithm insert, but entries from higher-level
//       nodes must be placed higher in the tree so that leaves of their
//       dependent subtrees will be on the same level as leaves of the
//       main tree.
//--------------------------------------------------------------------------
template <class Y> MY_INLINE CubitStatus RTreeNode<Y>::condense_tree(RTreeNode<Y> *l,
                                                           RTreeNode<Y> *root,
                                                           RTreeNode<Y> *&new_root )
{
  int ii;
  new_root = NULL;
    //CT1)
  RTreeNode<Y> *n = l, *p;
  DLIList <RTreeNode<Y>*> set_q;
    //CT2)
  while ( n != root )
  {
    p = n->get_parent();
    if ( n->num_children() < minChildren )
    {
        //CT3
        //take these children and add them to set_q.
      for ( ii = 0;ii < n->num_children(); ii++ )
        set_q.append(n->get_child(ii));
        //remove n from p.
      p->remove_child(n);
        //delete n.
      delete n;
        //now continue on.
    }
    else
    {
        //CT4
      n->recalc_b_box();
    }
      //CT5
    n = p;
  }
    //now reinsert all nodes in set_q.
  RTreeNode<Y> *curr_node, *temp_root;
  temp_root = root;
  for (ii = set_q.size(); ii > 0; ii-- )
  {
    curr_node = set_q.get_and_step();
    temp_root->insert(curr_node, new_root);
    if ( new_root != NULL )
      temp_root = new_root;
  }
  if ( temp_root != root )
    new_root = temp_root;
  return CUBIT_SUCCESS;
}



  
