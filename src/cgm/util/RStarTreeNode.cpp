//---------------------------------------------------------------------------
// Class Name:  RStarTreeNode
// Description: Node of Rectangle tree.  Contians many of the
//              required functions for building the tree and traversing it.
// The algorithm was taken from the following paper:
//	      Norbert Beckmann, H. Kriegel, R. Schnieder, and B. Seegar,
//              "The R*-tree: An Efficient and Robust Access Method
//              for Points and Rectangles", Proceedings of ACM SIGMOD
//              Int'l. Conf. on Management of Data, pp. 322-331, 1990.

// Creation Date: 7/21/02
// Owner:  David R. White
//---------------------------------------------------------------------------

//---------------------------------
//Include Files
//---------------------------------
#include "RStarTreeNode.hpp"
#include "DLIList.hpp"
#include "CpuTimer.hpp"
//---------------------------
//Initialize Static Members
//---------------------------

#ifdef INLINE_TEMPLATES
#define MY_INLINE inline
#else
#define MY_INLINE
#endif
static int id = 0;

template <class Y> MY_INLINE RStarTreeNode<Y>::RStarTreeNode (Y data, double tol,
                                                      int max_children,
                                                      int min_children)
{
  myId = id++;
  maxChildren = max_children;
  minChildren = min_children;
  myChildrenNodes = new RStarTreeNode<Y>* [maxChildren];
  int ii;
  for ( ii = 0; ii < maxChildren; ii++ )
    myChildrenNodes[ii] = (RStarTreeNode<Y>*) NULL;
  if ( data == NULL )
  {
    PRINT_ERROR("Building RTree with null data is not allowed!\n");
    assert(data != NULL);
  }
  myData = data;
  myLevel = DATA_RSTARNODE;
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
template <class Y> MY_INLINE RStarTreeNode<Y>::RStarTreeNode (CubitBox &bounding_box,
                                                      int max_children,
                                                      int min_children)
{  
  myId = id++;
  maxChildren = max_children;
  minChildren = min_children;
  myBoundingBox = new CubitBox(bounding_box);
  myChildrenNodes = new RStarTreeNode<Y>* [maxChildren];
  int ii;
  for ( ii = 0; ii < maxChildren; ii++ )
    myChildrenNodes[ii] = (RStarTreeNode<Y>*) NULL;
  myData = NULL;
  myLevel = UNSET_RSTARNODE;
  myParent = NULL;
  nextChildIndex = 0;
  markedFlag = 0;
  distIsBox = 1;
  myDist = CUBIT_DBL_MAX;
}
//-----------------------------------------------------------
// Destructor
//-----------------------------------------------------------
template <class Y> MY_INLINE RStarTreeNode<Y>::~RStarTreeNode()
{
  if ( myChildrenNodes )
    delete [] myChildrenNodes;
  if ( myBoundingBox )
    delete myBoundingBox;
}
template <class Y> MY_INLINE void RStarTreeNode<Y>::validate_tree(int print)
{
	int ii;
   if (print )
          {
            PRINT_INFO("Parent %d: Children: ", myId);
            for ( ii = 0; ii < num_children(); ii++ )
            {
              RStarTreeNode<Y> *curr_node = myChildrenNodes[ii];
              PRINT_INFO("%d ", curr_node->myId);
            }
            PRINT_INFO("\n");
          }
	for ( ii = 0; ii < num_children(); ii++ )
	{
		RStarTreeNode<Y> *curr_node = myChildrenNodes[ii];
		assert (curr_node->get_parent() == this );
		curr_node->validate_tree(print);
	}
    return;
}
//-----------------------------------------------------------
// Algorithm: insert
//  Insert a new index entry e into an R-Tree.
//     I1. [Find postiion for new record] Invoke choose_sub_tree to select
//        a leaf node in l, in which to place e.
//     I2. [Add record to leaf node].a) If l has room for
//        another entry, install E. b) Otherwise invoke overflow_treatment to
//        insert e by reinserting in a different order or spliting l.
//     I3. [Propogate changes upward] Invoke adjust_tree on l, also passing ll
//        if a split was performed.
//     I4. [Grow Tree Taller] If node split propogation caused the root
//         to split create a new root whose children are the two resulting
//         nodes.
//-----------------------------------------------------------
template <class Y> MY_INLINE CubitStatus RStarTreeNode<Y>::insert(RStarTreeNode<Y> *e,
                                                                  RStarTreeNode<Y> *&new_root,
                                                                  int *overflow_flags,
                                                                  int levels)
{
  int print1=0;
  if ( print1 )
    this->validate_tree(print1);
  CubitStatus stat;
  new_root = NULL;//only set this if the root node changes.  Assume
    //that this RStarTreeNode object is the root...
  RStarTreeNode<Y> *root = this;
  
    // I1. Invoke choose_sub_tree to select a leaf node l in which to place
    //e
  RStarTreeNode<Y> *l = choose_sub_tree(this, e);
  assert(l->get_parent() != NULL || l == this );
  
    //just test.
    // make sure l is not null.
    // make sure l is one level above e.
  if ( l ==  NULL || l->get_leaf_level() != (e->get_leaf_level() + 1) )
  {
    PRINT_ERROR("Choosing leaf for inseartion into rtree failed.\n");
    return CUBIT_FAILURE;
  }
  RStarTreeNode<Y> *ll = NULL;
    //I2 a) If l has room for another entry install e.
  if ( l->can_add() )
  {
    l->add_child(e, CUBIT_TRUE);
  }
  else
  {
      //Call the overflow.
    stat = overflow_treatment(l, e, ll, root, new_root,
                              overflow_flags, levels);
    if ( stat != CUBIT_SUCCESS )
      return stat;
  }
    //adjust the bounding boxes and if needed
    //create a new root...
  assert(l->get_parent() != NULL || l == root || l == new_root );
  int print = 0;
  if ( new_root == NULL && print )
	  this->validate_tree(print);
  else if ( new_root != NULL && print )
	  new_root->validate_tree(print);

  stat = adjust_tree(l, ll, root, new_root,
                     overflow_flags, levels);
  if ( stat!= CUBIT_SUCCESS )
    return stat;
  return CUBIT_SUCCESS;
}
//--------------------------------------------
// Algorithm: overflow_treatment
// Decides whether or not to do a reinsert or
// a split.  Basically e should go into l, but
// there is no more room for it...
//--------------------------------------------
template <class Y> MY_INLINE
CubitStatus RStarTreeNode<Y>::overflow_treatment( RStarTreeNode<Y>* l,
                                                  RStarTreeNode<Y>* e,
                                                  RStarTreeNode<Y> *&ll,
                                                  RStarTreeNode<Y> *root,
                                                  RStarTreeNode<Y> *&new_root,
                                                  int *overflow_flags, int levels)
{
  assert(l->get_parent() != NULL || l == root || l == new_root );
  
  CubitStatus stat;
    //Test is this level is not the root level.
  if ( l->get_leaf_level() != (levels-1) && overflow_flags[l->get_leaf_level()] == 0)
  {
      //mark this level as having been reinserted...
    overflow_flags[l->get_leaf_level()] = 1;
    stat = reinsert(l,e,root,new_root,overflow_flags,levels);
    if ( stat != CUBIT_SUCCESS )
      return stat;
  }
  else
  {
    stat = split_node(l, e, ll);
    if ( stat != CUBIT_SUCCESS )
      return stat;
  }
  return stat;
}
//--------------------------------------------
// Private Algorithm: sort_center_distance
// Used for sorting in decreasing order (max first)
// rtree nodes based on their distance value.
//--------------------------------------------
template <class Y> MY_INLINE
int RStarTreeNode<Y>::sort_center_distance( RStarTreeNode<Y> *&n_1,
                                            RStarTreeNode<Y> *&n_2 )
{
  if ( n_1->get_dist() > n_2->get_dist() )
    return -1;
  else if ( n_1->get_dist() < n_2->get_dist() )
    return 1;
  else
    return 0;
}
//--------------------------------------------
// Private Algorithm: reinsert
// This algorithm chooses p nodes to remove
// from l and reinsert them into the tree.
//--------------------------------------------
template <class Y> MY_INLINE
CubitStatus RStarTreeNode<Y>::reinsert(RStarTreeNode<Y>* l,
                                   RStarTreeNode<Y>* e,
                                   RStarTreeNode<Y> *root,
                                   RStarTreeNode<Y> *&new_root,
                                   int *overflow_flags, int levels)
{
  DLIList <RStarTreeNode<Y>*> ordered_entries;
  RStarTreeNode<Y> *curr_node;
  CubitBox big_bound = l->bounding_box();
  big_bound |= e->bounding_box();
  CubitVector center_big = big_bound.center();
  CubitVector center_curr;
  double dist_sq;
  int ii;
  for ( ii = 0; ii < maxChildren; ii++)
  {
    curr_node = l->myChildrenNodes[ii];
    center_curr = curr_node->bounding_box().center();
    dist_sq = (center_curr-center_big).length_squared();
    curr_node->set_dist(dist_sq);
    ordered_entries.append(curr_node);
  }
  center_curr = e->bounding_box().center();
  dist_sq = (center_curr-center_big).length_squared();
  e->set_dist(dist_sq);
  ordered_entries.append(e);
  ordered_entries.sort( sort_center_distance );
    //Make sure the sorting worked...
  if (ordered_entries.get()->get_dist() < ordered_entries.next()->get_dist())
  {
    PRINT_ERROR("Sorting failed in R*Tree.\n");
    assert(0);
    return CUBIT_FAILURE;
  }
    //Calculate P.  The rstar tree says to use 30% of M.
    //I'll round up...
  double P = .3*maxChildren;
  int p = (int) (P+0.5);
  DLIList <RStarTreeNode<Y>*> reinsert_nodes;
  for ( ii = 0; ii < p; ii++ )
  {
    reinsert_nodes.append(ordered_entries.get_and_step());
  }
    //Now reverse the reinsert nodes, inorder to reinsert
    //the minimum distance ones first as the paper says
    //this far outperforms the max ones.
  reinsert_nodes.reverse();

    //remove these nodes from l.
  CubitBoolean e_reinserted = CUBIT_FALSE;
  
  for ( ii = 0; ii < p; ii++ )
  {
    curr_node = reinsert_nodes.get_and_step();
      //remember e wasn't part of l anyways...
    if ( curr_node == e )
    {
      e_reinserted = CUBIT_TRUE;
      continue;
    }
    l->remove_child(curr_node);
    curr_node->set_parent(NULL);
  }
    //ressize the bounding box.
  if ( !e_reinserted )
  {
    l->add_child(e, CUBIT_FALSE);
  }
  l->recalc_b_box();
  CubitStatus stat;
  RStarTreeNode<Y> *changed_root = NULL;
  for ( ii = 0; ii < p; ii++ )
  {
    curr_node = reinsert_nodes.get_and_step();
    stat = root->insert(curr_node, new_root,
                        overflow_flags, levels);
    if ( stat != CUBIT_SUCCESS || curr_node->get_parent() == NULL)
    {
      PRINT_ERROR("RStarTree::reinsert insertion failed.\n");
      return stat;
    }
    if ( new_root != NULL )
    {
      changed_root = new_root;
      root = new_root;
    }
  }
    //if the root was split during this, like at one of the middle nodes,
    //new root would get reset to null again.  Soo, luckily we saved that
    //change!  Reassign changed_root to new_root.
  if ( changed_root != NULL )
    new_root = changed_root;
  return CUBIT_SUCCESS;
}

//--------------------------------------------
// Algorithm: choose_sub_tree: Select a leaf node in which to place
// a new index entry e.  Recursive search the subtrees of n
// until n is a leaf node.
//----------------------------------------------
template <class Y> MY_INLINE
RStarTreeNode<Y>* RStarTreeNode<Y>::choose_sub_tree( RStarTreeNode<Y>* n,
                                                     RStarTreeNode<Y>* e )
{
    //If n is a leaf node, or one level greater than e,
    //we are done.
  if ( n->get_leaf_level() == (e->get_leaf_level() + 1) )
    return n;
    //Now choose the entry f in n (children of n that is)
    //If the children of n are leaf nodes, then find the entry f in n
    //  whose rectangle needs least overlap enalargement to include the new data
    //  rectangle.  Resolve ties by choosing the entry whose rectangle needs least
    //  are enlargement, then the entry with the rectangle of smallest area.
    //Else Choose the entry f in n whose rectangle needs least area enlargment to include the new
    //data rectangle.  Resolve ties by choosing the entry with the rectangle of smallest area.
  double min_enlargement = CUBIT_DBL_MAX, curr_enlargement;
  double min_overlap = CUBIT_DBL_MAX, curr_overlap;
  RStarTreeNode<Y> *curr_node;
  int child_index = -1;
  int ii;
  CubitBox bounding_box;
  for(ii = 0; (ii < maxChildren) && (n->myChildrenNodes[ii] != NULL); ii++  )
  {

    curr_node = n->myChildrenNodes[ii];
	assert(curr_node->get_parent() != NULL );
    if ( curr_node->get_leaf_level() == (e->get_leaf_level() + 1) )
    {
      curr_overlap = calc_overlap(curr_node, e);
      if ( curr_overlap <= min_overlap )
      {
        if ( curr_overlap == min_overlap && child_index >= 0 )
        {
          double curr_enl = calc_enlargement(curr_node, e);
          double best_enl = calc_enlargement(n->get_child(ii), e);
          if ( curr_enl > best_enl )
            continue;
          if ( curr_enl == best_enl )
          {
              //only reset if the curr_node has a smaller volume.
            double curr_vol = volume(curr_node);
            double old_vol = volume(n->myChildrenNodes[child_index]);
            if ( old_vol < curr_vol )
              continue;
          }
        }
        child_index = ii;
        min_overlap = curr_overlap;
      }
    }
    else
    {
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
  }
    //do error checking...
  if ( child_index == -1 || child_index >= maxChildren )
    return (RStarTreeNode<Y>*)NULL;
  RStarTreeNode<Y> *f = n->myChildrenNodes[child_index];
    //Now continue on...
  curr_node = choose_sub_tree(f,e);
  return curr_node;
}
//----------------------------------------------------------------------
// calc_overlap:  Calculate the total overlap between the add_to and the
//                children of current.
//----------------------------------------------------------------------
template <class Y> MY_INLINE double RStarTreeNode<Y>::calc_overlap(RStarTreeNode<Y> *current,
                                                                   RStarTreeNode<Y> *add_to)
{
  int ii, jj;
  CubitBox add_to_box = add_to->bounding_box();
  double total_volume = 0.0;
    //calculate the total overlap currently.
  CubitBox curr_child_box, other_child_box, temp_box;
  for ( ii = 0; ii < current->num_children(); ii++ )
  {
    curr_child_box = current->get_child(ii)->bounding_box();
    for ( jj = 0; jj < current->num_children(); jj++ )
    {
      if ( ii == jj )
        continue;
      temp_box = curr_child_box;
      other_child_box.reset(current->get_child(jj)->bounding_box());
      temp_box &= other_child_box;
      total_volume += volume(temp_box);
    }
  }
  double prev_total = total_volume;
    //add to it the overlap that would occur.
  for ( ii = 0; ii < current->num_children(); ii++ )
  {
    curr_child_box.reset( current->get_child(ii)->bounding_box());
    curr_child_box &= add_to_box;
    total_volume += volume(curr_child_box);
  }
    //now find the overlap enlargment, total - prev_total...
  return (total_volume-prev_total);
}

//----------------------------------------------------------------------
// calc_enlargement:  Calculate the enlargement required for increasing
// the bounding box of current so that it would encapsulate the bounding
// box of add_to.  So to do that, create the union of the two bounding
// boxes, then of that supper box subtrace the volume of the current.
// The result should be the volumetric difference between how much
// current has and how much it would need be or the enlargement.
//----------------------------------------------------------------------
template <class Y> MY_INLINE
double RStarTreeNode<Y>::calc_enlargement(RStarTreeNode<Y> *current, RStarTreeNode<Y> *add_to )
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
template <class Y> MY_INLINE
double RStarTreeNode<Y>::calc_enlargement(CubitBox &current, CubitBox &add_to )
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
template <class Y> MY_INLINE
CubitStatus RStarTreeNode<Y>::adjust_tree(RStarTreeNode<Y> *l, RStarTreeNode<Y> *ll,
                                          RStarTreeNode<Y> *root_node,
                                          RStarTreeNode<Y> *&new_root,
                                          int *overflow_flags,
                                          int levels)
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
      new_root = new RStarTreeNode<Y>(root_box, maxChildren, minChildren);
      int new_level = l->get_leaf_level() + 1;
      new_root->set_leaf_level(new_level);
      new_root->add_child(l, CUBIT_TRUE);
      new_root->add_child(ll, CUBIT_TRUE);
      return CUBIT_SUCCESS;
    }
  }
  else if ( l == new_root && ll == NULL )
  {
	  return CUBIT_SUCCESS;
  }
  else if ( l == new_root && ll != NULL )
  {
   //Create a new root node and store l and ll there
    CubitBox root_box = l->bounding_box();
    root_box |= ll->bounding_box();
    new_root = new RStarTreeNode<Y>(root_box, maxChildren, minChildren);
    int new_level = l->get_leaf_level() + 1;
    new_root->set_leaf_level(new_level);
    new_root->add_child(l, CUBIT_TRUE);
    new_root->add_child(ll, CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }

  RStarTreeNode<Y> *parent_node = l->get_parent();
  RStarTreeNode<Y> *new_group = NULL;
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
      stat = overflow_treatment(parent_node, ll, new_group, root_node, new_root,
                                overflow_flags, levels);
      if ( stat != CUBIT_SUCCESS )
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
  if ( parent_node->get_parent() == NULL && 
       parent_node != root_node &&
       parent_node != new_root )
  {
    PRINT_INFO("level = %d\n", parent_node->get_leaf_level());
    PRINT_INFO("levels = %d\n", levels);
    PRINT_ERROR("parent_node (%d) == NULL\n", parent_node->myId);
    PRINT_ERROR("And l (%d) ", l->myId);
    assert(0);
  }
  stat = adjust_tree(parent_node, new_group, root_node, new_root,
                     overflow_flags, levels);
  if ( stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
  return CUBIT_SUCCESS;
}
template <class Y> MY_INLINE
int RStarTreeNode<Y>::sort_high_x(RStarTreeNode<Y> *&n_1,
                                  RStarTreeNode<Y> *&n_2 )
{
  CubitVector n_1_high = n_1->bounding_box().maximum();
  CubitVector n_2_high = n_2->bounding_box().maximum();

  if ( n_1_high.x() < n_2_high.x() )
    return -1;
  else if ( n_1_high.x() == n_2_high.x() )
    return 0;
  else
    return 1;
}
template <class Y> MY_INLINE
int RStarTreeNode<Y>::sort_high_y(RStarTreeNode<Y> *&n_1,
                                  RStarTreeNode<Y> *&n_2 )
{
  CubitVector n_1_high = n_1->bounding_box().maximum();
  CubitVector n_2_high = n_2->bounding_box().maximum();

  if ( n_1_high.y() < n_2_high.y() )
    return -1;
  else if ( n_1_high.y() == n_2_high.y() )
    return 0;
  else
    return 1;
}
template <class Y> MY_INLINE
int RStarTreeNode<Y>::sort_high_z(RStarTreeNode<Y> *&n_1,
                                  RStarTreeNode<Y> *&n_2 )
{
  CubitVector n_1_high = n_1->bounding_box().maximum();
  CubitVector n_2_high = n_2->bounding_box().maximum();

  if ( n_1_high.z() < n_2_high.z() )
    return -1;
  else if ( n_1_high.z() == n_2_high.z() )
    return 0;
  else
    return 1;
}
template <class Y> MY_INLINE
int RStarTreeNode<Y>::sort_low_x(RStarTreeNode<Y> *&n_1,
                                 RStarTreeNode<Y> *&n_2 )
{
  CubitVector n_1_low = n_1->bounding_box().minimum();
  CubitVector n_2_low = n_2->bounding_box().minimum();

  if ( n_1_low.x() < n_2_low.x() )
    return -1;
  else if ( n_1_low.x() == n_2_low.x() )
    return 0;
  else
    return 1;
}
template <class Y> MY_INLINE
int RStarTreeNode<Y>::sort_low_y(RStarTreeNode<Y> *&n_1,
                                 RStarTreeNode<Y> *&n_2 )
{
  CubitVector n_1_low = n_1->bounding_box().minimum();
  CubitVector n_2_low = n_2->bounding_box().minimum();

  if ( n_1_low.y() < n_2_low.y() )
    return -1;
  else if ( n_1_low.y() == n_2_low.y() )
    return 0;
  else
    return 1;
}
template <class Y> MY_INLINE
int RStarTreeNode<Y>::sort_low_z(RStarTreeNode<Y> *&n_1,
                                 RStarTreeNode<Y> *&n_2 )
{
  CubitVector n_1_low = n_1->bounding_box().minimum();
  CubitVector n_2_low = n_2->bounding_box().minimum();

  if ( n_1_low.z() < n_2_low.z() )
    return -1;
  else if ( n_1_low.z() == n_2_low.z() )
    return 0;
  else
    return 1;
}
  
//------------------------------------------------------------------
// Algorithm: split_node
// This function is rather tricky since it really isn't well
// described Beckmann's paper very well.  I looked at other online
// docs and descriptions and came to the current implementation.
// As I understand it the current function does the following:
// First descide which axis the nodes should be split along.
// To accomplish this the nodes that are going to be split (the
// children of l and the node e), are added to two lists.  The lists
// are then sorted according to their high and low values along
// the three axis.  Then for each  each high and low,
// d distributions are created with possible groupings.
// Where d = (maxChildren -2*minChildren +2). These distributions are
// then used to calculate the total margin for each axis.  The axis
// with the minimum margin is selected.  While calculating the margins
// for each distribution, the best "distribution" for each axis is also
// selected.  The best distribution will be the one that has the minimum
// overlap over the entire set of distributions, and high and low sets.
// When the axis is chosen, the correct distribution is then also stored
// or known.  The function then splits l into l and ll.
//------------------------------------------------------------------
template <class Y> MY_INLINE
CubitStatus RStarTreeNode<Y>::split_node( RStarTreeNode<Y> *l,
                                          RStarTreeNode<Y> *e,
                                          RStarTreeNode<Y> *&ll )
{
    int ii;
    //create a new list containing all the nodes we want to split.
      //create two lists.
  DLIList <RStarTreeNode<Y>*> ordered_low, ordered_high;
  for ( ii = 0; ii < maxChildren; ii++)
  {
    ordered_low.append(l->myChildrenNodes[ii]);
    ordered_high.append(l->myChildrenNodes[ii]);
  }
  ordered_low.append(e);
  ordered_high.append(e);
    //the input list contains all of the nodes.

  int d = maxChildren - 2*minChildren + 2;

    //Now do the first step, choose the split axis.
    //loop over each dimension.
  double local_margin, min_margin = CUBIT_DBL_MAX;
  DLIList<RStarTreeNode<Y>*> best_group_1, best_group_2;
  
  for ( ii = 0; ii < 3; ii++ )
  {
      //Sort the lists according to the high and low dimension.
      //Both lists are ordered lowest to highest however.
    switch(ii)
    {
      case(0):
          //this is the x dimension.
        ordered_low.sort(sort_low_x);
        ordered_high.sort(sort_high_x);
        break;
      case(1):
          //this is the y dimension.
        ordered_low.sort(sort_low_y);
        ordered_high.sort(sort_high_y);
        break;
      case(2):
          //this is the z dimension.
        ordered_low.sort(sort_low_z);
        ordered_high.sort(sort_high_z);
        break;
    }
      //Now loop over the distributions and sum the margins for the
      //different distributions.  The axis where the sum of its margins
      //is minimal is the correct axis.
    int k;
    local_margin = 0.0;
    double min_overlap = CUBIT_DBL_MAX;
    double min_volume = CUBIT_DBL_MAX;
    DLIList<RStarTreeNode<Y>*> group_1_low, group_1_high, local_best_1;
    DLIList<RStarTreeNode<Y>*> group_2_low, group_2_high, local_best_2;
      //just do this so that the code can look familar with the paper.
    int m = minChildren;
    int M = maxChildren;
      //Also determine with distribution is best among these in this axis.
      //Store those groups in case this axis is optimum.
    for ( k = 0; k < d; k++ )
    {
        //build the 4 groups...
      int jj;
      group_1_low.clean_out();
      group_2_low.clean_out();
      group_1_high.clean_out();
      group_2_high.clean_out();
      for ( jj = 0; jj < (m-1+k); jj++ )
      {
        group_1_low.append(ordered_low.next(jj));
        group_1_high.append(ordered_high.next(jj));
      }
      for ( jj = (m-1+k); jj < (M+1); jj++ )
      {
        group_2_low.append(ordered_low.next(jj));
        group_2_high.append(ordered_high.next(jj));
      }
      assert(group_1_low.size() + group_2_low.size() == M+1 );
        //Okay we have the groups.  Now calculate the metrics.
        //First find the bounding boxes for the groups.
      CubitBox group_1_low_box = super_box(group_1_low);
      CubitBox group_2_low_box = super_box(group_2_low);
      CubitBox group_1_high_box = super_box(group_1_high);
      CubitBox group_2_high_box = super_box(group_2_high);
      
      local_margin += margin(group_1_low_box);
      local_margin += margin(group_2_low_box);
      local_margin += margin(group_1_high_box);
      local_margin += margin(group_2_high_box);
        //Okay now that we have the margin, that is the portion of the
        //code for choosing the correct axis.  Now make sure if this axis
        //is the right one, we find the right distribution.

      double overlap_low, overlap_high;
        //remember &= is the overlap or intersection and the volume calculates
        //the volume of the overlap or intersection.
      overlap_low = volume(group_1_low_box &= group_2_low_box);
      overlap_high = volume(group_1_high_box &= group_2_high_box);
      CubitBoolean use_low = (overlap_low < overlap_high)? CUBIT_TRUE : CUBIT_FALSE;
      double temp_overlap = use_low ? overlap_low : overlap_high;
        //Choose the best distribution based on the mininum distribution
      if ( temp_overlap < min_overlap )
      {
        min_overlap = temp_overlap;
        if ( use_low )
        {
          local_best_1 = group_1_low;
          local_best_2 = group_2_low;
        }
        else
        {
          local_best_1 = group_1_high;
          local_best_2 = group_2_high;
        }
      }
        //break ties based on the smallest volumes.
      else if ( temp_overlap == min_overlap )
      {
          //supposed to resolve this by choosing the one with the minimum area.
        double tmp_vol;
        if ( use_low ){
          tmp_vol = volume(group_1_low_box);
          tmp_vol += volume(group_2_low_box);
        }
        else {
          tmp_vol = volume(group_1_high_box);
          tmp_vol += volume(group_2_high_box);
        }
        if ( tmp_vol < min_volume )
        {
          min_volume = tmp_vol;
          if ( use_low )
          {
            local_best_1 = group_1_low;
            local_best_2 = group_2_low;
          }
          else
          {
            local_best_1 = group_1_high;
            local_best_2 = group_2_high;
          }
        }
      }
    }
      //After the margin has been sumed for the entire distributions,
      //choose the axis with the min margin.  Note I'm not storing the
      //axis because for each distribution I'm also chosing the local
      //best based on overlap.  Store that local best as the overal all
      //best.  It only gets stored if the axis is optimum...
    if ( local_margin < min_margin )
    {
      min_margin = local_margin;
      best_group_1 = local_best_1;
      best_group_2 = local_best_2;
    }
  }
    //Okay now we have the groups.  Clean out l, create ll.
  l->flush(best_group_1.get()->bounding_box());
  l->add_child(best_group_1.get_and_step(), CUBIT_FALSE);
  l->set_leaf_level(e->get_leaf_level() + 1);
  for ( ii = 1; ii < best_group_1.size(); ii++ )
    l->add_child(best_group_1.get_and_step(), CUBIT_TRUE);
  ll = new RStarTreeNode<Y>(best_group_2.get()->bounding_box(),
                        maxChildren, minChildren);
  ll->add_child(best_group_2.get_and_step(), CUBIT_FALSE);
  ll->set_leaf_level(l->get_leaf_level());
  for ( ii = 1; ii < best_group_2.size(); ii++ )
    ll->add_child(best_group_2.get_and_step(), CUBIT_TRUE);

  return CUBIT_SUCCESS;
}
//-----------------------------------------------
//Private Function: Margin
//  Calculates the margin of bounding box.
//-----------------------------------------------
template <class Y> MY_INLINE
double RStarTreeNode<Y>::margin(CubitBox &bounding_box)
{
  double margin = 4*(bounding_box.x_range() + bounding_box.y_range()
                     + bounding_box.z_range());
  return margin;
}
//-----------------------------------------------
//Private Function: super_box
//  Calculates the overall bounding box of the rtree
//  nodes in the list.
//-----------------------------------------------
template <class Y> MY_INLINE
CubitBox RStarTreeNode<Y>::super_box(DLIList<RStarTreeNode<Y>*> &node_list)
{
  int ii;
  CubitBox bounding_box = node_list.get_and_step()->bounding_box();
  for ( ii = 1; ii < node_list.size(); ii++ )
  {
    bounding_box |= node_list.get_and_step()->bounding_box();
  }
  return bounding_box;
}

template <class Y> MY_INLINE void RStarTreeNode<Y>::flush( CubitBox &new_box )
{
  int ii;
  nextChildIndex = 0;
  for ( ii = 0; ii < maxChildren; ii++ )
    myChildrenNodes[ii] = NULL;
  delete myBoundingBox;
  myBoundingBox = new CubitBox(new_box);
}
template <class Y> MY_INLINE void RStarTreeNode<Y>::add_child(RStarTreeNode<Y> *child_node,
                                                CubitBoolean recalc_b_box)
{
  assert(nextChildIndex < maxChildren && child_node != NULL );
  myChildrenNodes[nextChildIndex] = child_node;
    //update the bounding box. by uniting with child node...
  if ( recalc_b_box )
  {
    CubitBox *old_box = myBoundingBox;
    myBoundingBox = new CubitBox( *old_box |= child_node->bounding_box());
    delete old_box;
  }
  nextChildIndex++;
  child_node->set_parent(this);
}
template <class Y> MY_INLINE CubitBoolean RStarTreeNode<Y>::can_add()
{
  if (nextChildIndex >= maxChildren )
    return CUBIT_FALSE;
  else
    return CUBIT_TRUE;
}
template <class Y> MY_INLINE int RStarTreeNode<Y>::space_left()
{
  return maxChildren - nextChildIndex;
}
template <class Y> MY_INLINE void RStarTreeNode<Y>::recalc_b_box()
{
  if(myLevel == DATA_RSTARNODE )
    return;
  int ii;
  CubitBox temp_box;
  CubitBoolean is_first = CUBIT_TRUE;
  for ( ii = 0; ii < nextChildIndex; ii++ )
  {
    if ( is_first )
    {
      is_first = CUBIT_FALSE;
      temp_box = myChildrenNodes[ii]->bounding_box();
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
template <class Y> MY_INLINE CubitBoolean RStarTreeNode<Y>::remove( Y e,
                                                      RStarTreeNode<Y> *&new_root,
                                                      CubitBoolean &delete_root)
{
    //D1) Find node containting record.
  RStarTreeNode<Y> *l = NULL;
  CubitBox my_box = e->bounding_box();
  CubitStatus stat = find_leaf(e, my_box, this, l);
  if ( l == NULL || stat != CUBIT_SUCCESS )
    return CUBIT_FALSE;
    //Now l is the RStarTreeNode that holds the actual data (a DATA_RSTARNODE)
    //not a leaf.  This was done for efficiency.
  RStarTreeNode<Y> *data_node = l;
  l = data_node->get_parent();
    //D2) [Delete record]  Remove e from l.
  
    //remove the data node from the children and delete
    //the node.
  l->remove_child(data_node);
  delete data_node;

    //D3) [Propogate Changes].
  stat = condense_tree(l, this, new_root);
    //D4) [Shorten the tree].
  RStarTreeNode<Y> *root = this;
  if ( new_root != NULL )
    root = new_root;
  if ( root->num_children() == 1 )
  {
    new_root = root->get_child(0);
    new_root->set_parent((RStarTreeNode<Y>*)NULL);
    delete_root = CUBIT_TRUE;
  }
  return CUBIT_TRUE;
}
template <class Y> MY_INLINE CubitStatus RStarTreeNode<Y>::find_leaf( Y e,
                                                        CubitBox &e_box,
                                                        RStarTreeNode<Y> *t,
                                                        RStarTreeNode<Y> *&l )
{
  int ii;
  CubitStatus stat;
  l = NULL;
  int loop_size = t->num_children();
  RStarTreeNode<Y> *curr_node;
  if ( t->get_leaf_level() > LEAF_RSTARNODE )
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
        stat = find_leaf(e, e_box,curr_node,l);
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
template <class Y> MY_INLINE CubitBoolean RStarTreeNode<Y>::remove_child( RStarTreeNode<Y> *child )
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
template <class Y> MY_INLINE CubitStatus RStarTreeNode<Y>::condense_tree(RStarTreeNode<Y> *l,
                                                           RStarTreeNode<Y> *root,
                                                           RStarTreeNode<Y> *&new_root )
{
  int ii;
  new_root = NULL;
    //CT1)
  RStarTreeNode<Y> *n = l, *p;
  DLIList <RStarTreeNode<Y>*> set_q;
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
  RStarTreeNode<Y> *curr_node, *temp_root;
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



  
