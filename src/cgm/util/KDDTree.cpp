//-----------------------------------------------------------------
//- Class:   KDDTree
//- Author:  Kevin Albrecht 
//- Created: 13 May 2003
//- Updated: 8 Feb 2004
//-
//- Description:
//-   Dynamic version of the k-d tree, where k=3.
//-
//- References:
//-
//-   Hanan Samet.  Design and Analysis of Spatial Data Structures.
//-      Addison-Wesley, Reading, MA, 1990.
//-
//-   Jon Louis Bentley. Multidimensional binary search trees used
//-     for associative searching. In Communications of the ACM,
//-     18(9), pages 509-517, September 1975.
//-----------------------------------------------------------------

#if !defined(TEMPLATE_DEFS_INCLUDED) || defined(INCLUDED_FROM_KDD_TREE_HEADER)

//---------------------------------
// Include Files
//---------------------------------

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "KDDTree.hpp"
#include "KDDTreeNode.hpp"
#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "PriorityQueue.hpp"

//---------------------------------
// Define Methods
//---------------------------------

//- Constructor
//-  * The following only applies if self-balancing is turned on: if
//-    selfBalancingDeletionTolerance is set to 0, then there is no limit
//-    to the number of nodes that can be marked for deletion at a time;
//-    otherwise, the tree will rebalance itself whenever the percentage
//-    of nodes on the tree marked for deletion is greater than the
//-    tolerance.
template <class Z> KDDTree<Z>::KDDTree (double tol, CubitBoolean selfBalancingOn,
                                        double selfBalancingDeletionTolerance,
                                        CubitBoolean randomOn)
{
  root = NULL;

  myTolerance = tol;
  mySelfBalancingOn = selfBalancingOn;
  myDeletionTolerance = selfBalancingDeletionTolerance;
  myMarkedNodes = 0;
  myRandomOn = randomOn;

  if (myRandomOn)
  {
    //// seed the random number generator
    srand( static_cast<unsigned>(time(NULL)) );
  }
}

//- Destructor
template <class Z> KDDTree<Z>::~KDDTree()
{
  int i;
  for (i = myAddList.size(); i > 0; i--)
  {
    delete myAddList.pop();
  }
  for (i = myNodeList.size(); i > 0; i--)
  {
    KDDTreeNode<Z> *node = myNodeList.pop();
    delete node;
  }
}

//- Immediately put all nodes on list onto the tree
template <class Z> CubitStatus KDDTree<Z>::dump_list ()
{
  while (myAddList.size() > 0)
  {
    KDDTreeNode<Z> *node = myAddList.pop();
    insert_node (node);
  }

  return CUBIT_SUCCESS;
}

//- "insert_node"
//- Dynamically insert the data into the k-d tree 
//- 
//- Algorithm INSERT (From Bentley):
//-   This algoritm is passed an object "data" of class "Z",
//-   which has a bounding_box() method.  If there is already
//-   a node in the tree with equal bounding box center point,
//-   it is put in the right subtree.
//-   I0. [Create new node] Create a node P with the bounding box
//-       specified, and set P.LEFT <- null, P.RIGHT <- null, and
//-       P.DISC <- null.
//-   I1. [Check for null tree] If ROOT = null, then set ROOT <- P
//-       and return CUBIT_SUCCESS; otherwise, set Q <- ROOT (Q
//-       will move down the tree).
//-   I2. [Compare] Compare the nodes and set the child in the
//-       correct direction to T.
//-   I3. [Move down] Set Q <- child of Q and go to I2.
//-   I4. [Insert new node in tree] Set the child of Q to P, then
//-       set the children of P to null. Set the discriminator of
//-       P to be the discriminator after that in Q.
//-
template <class Z> CubitStatus KDDTree<Z>::insert_node (KDDTreeNode<Z>* P)
{
  KDDTreeNode<Z> *F = NULL;      // father node
  KDDTreeNode<Z> *T;             // temp node

  if (root == NULL)
  {
    root = P;
    P->set_disc (DIMX);
  }
  else
  {
    T = root;
    DIRECTION direction = DIR_NULL;

    while (T != NULL)
    {
      F = T; // remember the father
      direction = P->compare (T);

      CubitVector tmin = T->safetyBox.minimum();
      CubitVector tmax = T->safetyBox.maximum();
      CubitVector pmin = P->safetyBox.minimum();
      CubitVector pmax = P->safetyBox.maximum();

      if (pmin.x() < tmin.x()) tmin.x (pmin.x());
      if (pmin.y() < tmin.y()) tmin.y (pmin.y());
      if (pmin.z() < tmin.z()) tmin.z (pmin.z());
      if (pmax.x() > tmax.x()) tmax.x (pmax.x());
      if (pmax.y() > tmax.y()) tmax.y (pmax.y());
      if (pmax.z() > tmax.z()) tmax.z (pmax.z());

      T->safetyBox.reset (tmin, tmax);

      T = T->get_child (direction);
    }

    F->set_child (P, direction);
    P->set_disc (F->next_disc());
    P->parent = F;
  }
  myNodeList.push (P);

  return CUBIT_SUCCESS;      
}

//- "modifind"
//- Rearrange the array around the median point.
//-
//- Description:
//-   This is the MODIFIND algorithm of V. Zabrodsky, but modified to choose a random
//-   pivot point.  This is in turn a modified version of the Hoare FIND algorithm for
//-   finding the median point.  Running time is O(n^2) in the worst case and O(n) in
//-   the average case.
//-
//-   Zabrodsky's algorithm:
//-     http://www.geocities.com/zabrodskyvlada/aat/a_modi.html
//-     http://www.geocities.com/zabrodskyvlada/3alg.html
//-
//- Results:
//-   Reordering of input array such that A[K] has the value it would have if A were
//-   sorted, L<=I<=K will imply A[I]<=A[K], and K<=I<=R will imply A[I]>=A[K]
template <class Z> int KDDTree<Z>::modifind (DIMENSION dim, int left, int right,
                                             KDDTreeNode<Z>* array[])
{
  int K = ((right - left + 1) / 2) + left + 1;
  int L = left + 1;
  int R = right + 1;
  int I, J;
  KDDTreeNode<Z>* node; // "X" in MODIFIND
  KDDTreeNode<Z>* temp; // temp for swapping; "W" in MODIFIND

  //// Choose pivot between left and right
  if (myRandomOn == CUBIT_TRUE)
  {
    int pivot = ( rand() % (right - left) ) + left; // create random pivot between left and right

    //// Swap array[pivot] and array[K]
    temp = array[pivot];
    array[pivot] = array[K - 1];
    array[K - 1] = temp;
  }

  while (L < R)
  {
    node = array[K - 1];
    I = L;
    J = R;
    while (! ((J < K) || (K < I)) )
    {
      if (dim == DIMX)
      {
        while (array[I - 1]->x < node->x) {
          I++;
        }
        while (node->x < array[J - 1]->x) {
          J--;
        }
      }
      else if (dim == DIMY)
      {
        while (array[I - 1]->y < node->y) {
          I++;
        }
        while (node->y < array[J - 1]->y) {
          J--;
        }
      }
      else
      {
        while (array[I - 1]->z < node->z) {
          I++;
        }
        while (node->z < array[J - 1]->z) { 
          J--; 
        }
      }

      //// Swap array[I] and array[J]
      temp = array[I - 1];
      array[I - 1] = array[J - 1];
      array[J - 1] = temp;

      I++;
      J--;
    }

    if (J < K)
    {
      L = I;
    }
    if (K < I)
    {
      R = J;
    }
  }

  return K - 1;
}

//- "balance" and "recursive_balance"
//- Create a balanced tree out of the nodes on the tree and in the Add List.
//- This is used to balance the tree manually; it is also called by the
//- find method when self-balancing is on.
//-
//- Description:
//-   This is the OPTIMIZE algorithm of Bentley.  Total running time is
//-   O(n lg n).  It uses the MODIFIND algorithm to find the median.
//-
//- Results:
//-   The tree produced will be balanced so that all leaf nodes occur on
//-   at most two adjacent levels.
//-
template <class Z> CubitStatus KDDTree<Z>::balance ()
{
  int arraypos = 0;
  KDDTreeNode<Z> ** array = new KDDTreeNode<Z>*[myAddList.size () + myNodeList.size ()];

  root = NULL;

  while (myAddList.size () > 0)
  {
    array[arraypos] = myAddList.pop();
    arraypos++;
  }
  while (myNodeList.size () > 0)
  {
    array[arraypos] = myNodeList.pop();

    if (array[arraypos]->valid == CUBIT_FALSE)
    {
      array[arraypos]->left = NULL;
      array[arraypos]->right = NULL;
      delete array[arraypos];
    }
    else
    {
      arraypos++;
    }
  }

  int left = 0;
  int right = (arraypos - 1);
  root = recursive_balance (DIMX, left, right, array, NULL);

  myMarkedNodes = 0;

  delete [] array;
  return CUBIT_SUCCESS;
}

template <class Z> KDDTreeNode<Z> *KDDTree<Z>::recursive_balance
 (DIMENSION dim, int left, int right, KDDTreeNode<Z>** array, KDDTreeNode<Z>* parent)
{
  if (left > right)
  {
    return NULL;
  }
  else
  {
    KDDTreeNode<Z>* P;
    int K;

    if (left != right)
    {
      K = modifind (dim, left, right, array);
      P = array[K];
    }
    else
    {
      K = left;
      P = array[left];
    }

    myNodeList.push (P);

    P->safetyBox.reset (P->boundingBox);
    for (int i = left; i <= right; i++)
    {
      CubitVector imin = array[i]->safetyBox.minimum();
      CubitVector imax = array[i]->safetyBox.maximum();
      CubitVector pmin = P->safetyBox.minimum();
      CubitVector pmax = P->safetyBox.maximum();

      if (imin.x() < pmin.x()) pmin.x (imin.x());
      if (imin.y() < pmin.y()) pmin.y (imin.y());
      if (imin.z() < pmin.z()) pmin.z (imin.z());
      if (imax.x() > pmax.x()) pmax.x (imax.x());
      if (imax.y() > pmax.y()) pmax.y (imax.y());
      if (imax.z() > pmax.z()) pmax.z (imax.z());

      P->safetyBox.reset (pmin, pmax);
    }

    DIMENSION nextDim;
    switch (dim)
    {
      case DIMX: nextDim = DIMY; break;
      case DIMY: nextDim = DIMZ; break;
      default:   nextDim = DIMX;
    }

    P->set_disc (dim);
    P->parent = parent;

    P->left = recursive_balance (nextDim, left, K - 1, array, P);
    P->right = recursive_balance (nextDim, K + 1, right, array, P);

    return P;
  }
}

//- Find the depth of the tree
template <class Z> int KDDTree<Z>::find_max_height ()
{
  int depth = 0, maxdepth = 0;
  recursive_find_max_height (root, depth, &maxdepth);

  return maxdepth;
}

//- Find the depth of the tree
template <class Z> void KDDTree<Z>::recursive_find_max_height
 (KDDTreeNode<Z> *the_root, int depth, int *maxdepth)
{
  if (the_root)
  {
    depth++;
    if (depth > *maxdepth)
    {
      *maxdepth = depth;
    }
    recursive_find_max_height (the_root->left, depth, maxdepth);
    recursive_find_max_height (the_root->right, depth, maxdepth);
  }
}

//- Add a node with the data to the list
template <class Z> CubitStatus KDDTree<Z>::add (Z data)
{
  KDDTreeNode<Z> *P = new KDDTreeNode<Z> (data);
  P->safetyBox.reset (data->bounding_box());

  myAddList.push (P);

  return CUBIT_SUCCESS;
}

//- Return a pointer to the node containing the specified data
template <class Z> KDDTreeNode<Z> *KDDTree<Z>::find_node_containing_data (KDDTreeNode<Z> *subtreeRoot, Z data)
{
  KDDTreeNode<Z> *T = new KDDTreeNode<Z>(data); // temp node to use in searching
  KDDTreeNode<Z> *P = subtreeRoot;              // node to delete
  DIRECTION D;

  //// Find the node
  while (P != NULL)
  {
    if ((P->boundingBox.minimum() == T->boundingBox.minimum()) &&
        (P->boundingBox.maximum() == T->boundingBox.maximum()))
    {
      if (P->valid == CUBIT_TRUE)
      {
        break; // the bounding boxes match and this node has not been deleted
      }
    }

    D = T->compare_with_equality (P);

    if (D == DIR_EITHER)
    {
      KDDTreeNode<Z> *leftResult = find_node_containing_data (P->get_child (DIR_LEFT), data);

      if (leftResult != NULL)
      {
        P = leftResult;
        break;
      }

      KDDTreeNode<Z> *rightResult = find_node_containing_data (P->get_child (DIR_RIGHT), data);
      
      if (rightResult != NULL)
      {
        P = rightResult;
        break;
      }

      P = NULL;
    }
    else
    {
      P = P->get_child (D);
    }
  }

  delete T;

  return P;
}

//- Remove the data member's entry in the tree. Returns CUBIT_TRUE
//- if item removed, CUBIT_FALSE if item not in tree.
template <class Z> CubitBoolean KDDTree<Z>::remove (Z data)
{
  //// If the Add List is not empty, action must be taken
  if (myAddList.size() > 0)
  {
    if (mySelfBalancingOn == CUBIT_TRUE) // self-balancing is on, so rebalance the tree
    {
      balance ();
    }
    else // self-balancing is off, so put everything in the Add List onto the tree
    {
      dump_list ();
    }
  }

  //// Tree is empty
  if (root == NULL)
  {
    return CUBIT_FALSE;
  }
  //// Tree is not empty
  else
  {
    KDDTreeNode<Z> *P = find_node_containing_data (root, data);

    if (P == NULL) // no matching node was found
    {
      return CUBIT_FALSE;
    }
    else // mark the matching node for deletion
    {
      if (P->valid == CUBIT_FALSE)
      {
        return CUBIT_FALSE; // this node was already deleted
      }

      P->valid = CUBIT_FALSE; // set the node to be deleted

      myMarkedNodes++;
      if (myDeletionTolerance != 0)
      {
        if ( (( static_cast<double>(myMarkedNodes) / myNodeList.size()) > myDeletionTolerance) &&
             (myMarkedNodes > 1)
           )
        {
          balance ();
        }
      }

      return CUBIT_TRUE;
    }
  }
}

//- Find members intersecting this range box
template <class Z> CubitStatus KDDTree<Z>::find (const CubitBox &range_box, DLIList <Z> &range_members)
{
  //// If the Add List is not empty, action must be taken
  if (myAddList.size() > 0)
  {
    if (mySelfBalancingOn == CUBIT_TRUE) // self-balancing is on, so rebalance the tree
    {
      balance ();
    }
    else // self-balancing is off, so put everything in the Add List onto the tree
    {
      dump_list ();
    }
  }

  //// Find all of the members of the tree that intersect this range_box
  if (root != NULL)
  {
    recursive_find (root, range_box, range_members);
  }

  return CUBIT_SUCCESS;
}

//- Recursively find members intersecting this range box (called by "find")
template <class Z> void KDDTree<Z>::recursive_find
  ( KDDTreeNode<Z> *rect_tree,
    const CubitBox &range_box,
    DLIList <Z> &range_members
  )
{
  //// Check for overlap with the safety box
  if ( ! range_box.overlap (myTolerance, rect_tree->safetyBox) )
  {
    return; // no overlap, return
  }

  //// Check for overlap with the bounding box
  if ( range_box.overlap (myTolerance, rect_tree->boundingBox) )
  {
    if (rect_tree->valid == CUBIT_TRUE)
    {
      range_members.append (rect_tree->data); // append the data to the list.
    }
  }

  if (rect_tree->left != NULL)
  {
    recursive_find (rect_tree->left, range_box, range_members);
  }
  if (rect_tree->right != NULL)
  {
    recursive_find (rect_tree->right, range_box, range_members);
  }

  return;
}

//- Finds the minimum distance squared between the given point and the box. If
//-  the point is on or in the box, the min distance is zero.
template <class Z> double KDDTree<Z>::min_dist_sq (CubitVector &q, CubitBox &b_box)
{
  CubitVector b_min = b_box.minimum();
  CubitVector b_max = b_box.maximum();
  CubitVector r;

  //// set "r" in the x-dim
  if (q.x () < b_min.x ())
  {
    r.x (b_min.x ());
  }
  else if (q.x () > b_max.x ())
  {
    r.x (b_max.x ());
  }
  else
  {
    r.x (q.x ());
  }
  
  //// set "r" in the y-dim
  if (q.y () < b_min.y ())
  {
    r.y (b_min.y ());
  }
  else if (q.y () > b_max.y ())
  {
    r.y (b_max.y ());
  }
  else
  {
    r.y (q.y ());
  }
  
  //// set "r" in the z-dim
  if (q.z () < b_min.z ())
  {
    r.z (b_min.z ());
  }
  else if (q.z () > b_max.z ())
  {
    r.z (b_max.z ());
  }
  else
  {
    r.z (q.z ());
  }
  
  double dist = (q-r).length_squared();

  return dist;
}

template <class Z> bool KDDTree<Z>::less_than_func (KDDTreeNode<Z> *&node_a,
                                                    KDDTreeNode<Z> *&node_b)
{
  if (node_a->get_dist() < node_b->get_dist ())
  {
    return true;
  }
  else
  {
    return false;
  }
}

//- "k_nearest_neighbor"
//- Find the K nearest neighbors to a point.
//-
//- Description:
//-   This algorithm is based on the best-first search.  The goal of this
//-   algorithm is to minimize the number of nodes visited by using the
//-   distance to each subtree's bounding box to avoid visiting subtrees
//-   which could not possibly contain one of the k nearest objects.
//-
template <class Z> CubitStatus KDDTree<Z>::k_nearest_neighbor
  (CubitVector &q, int k, double &closest_dist, DLIList<Z> &nearest_neighbors,
   typename KDDTree<Z>::DistSqFunc dist_sq_point_data
  )  
{
  //// Create the priority queues
  PriorityQueue<KDDTreeNode<Z>*> *queue = new PriorityQueue<KDDTreeNode<Z>*> (KDDTree<Z>::less_than_func);
  PriorityQueue<KDDTreeNode<Z>*> *queueTemp = new PriorityQueue<KDDTreeNode<Z>*> (KDDTree<Z>::less_than_func);


  KDDTreeNode<Z> *element = root;

  // push this node on the queue
  element->set_dist (min_dist_sq (q, element->safetyBox));
  element->set_dist_data (DD_SAFETY);
  queue->push (element);
  

  // if the k closest nodes on the tree are not leaf-nodes, expand the closest
  //   non-leaf node
  while ( !queue->empty() )
  {
    element = queue->top();
    queue->pop();

    if (element->get_dist_data() == DD_LEAF)
    {
      // this node is a leaf, so it can be pushed onto the temporary queue
      queueTemp->push (element);
    }
    else
    {
      // one of the top k nodes is a non-leaf node, so expand it
      if (element->left)
      {
        element->left->set_dist (min_dist_sq (q, element->left->safetyBox));
        element->left->set_dist_data (DD_SAFETY);
        queue->push (element->left);
      }
      if (element->right)
      {
        element->right->set_dist (min_dist_sq (q, element->right->safetyBox));
        element->right->set_dist_data (DD_SAFETY);
        queue->push (element->right);
      }
      element->set_dist (dist_sq_point_data (q, element->data));
      element->set_dist_data (DD_LEAF);
      queue->push (element);

      // take all the elements in the temporary queue and reinsert them into
      //   the actual queue
      while ( !queueTemp->empty() )
      {
        queue->push ( queueTemp->top() );
        queueTemp->pop ();
      }
    }

    if (queueTemp->size() == k)
    {
      // success-- place the k nodes into the nearest_neighbors list
      element = queueTemp->top();
      queueTemp->pop();
      closest_dist = element->get_dist();
      nearest_neighbors.append (element->data);

      while ( !queueTemp->empty() )
      {
        nearest_neighbors.append ( queueTemp->top()->data );
        queueTemp->pop();
      }
     
      return CUBIT_SUCCESS;
    }
  }
  return CUBIT_FAILURE;
}

#endif
