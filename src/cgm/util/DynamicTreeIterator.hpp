
#ifndef DYNAMIC_TREE_ITERATOR_HPP
#define DYNAMIC_TREE_ITERATOR_HPP

#include <functional>
#include "CubitDefines.h"
#include "DynamicDLIIterator.hpp"

// use this functor to iterate over first items in a map
template <class Pair, class First>
struct MapSelect1st : public std::unary_function<Pair, First>
{
  const First& operator()(const Pair& x) const 
  {
    return x.first;
  }
};

// use this functor to iterate over second items in a map
template <class Pair, class Second>
struct MapSelect2nd : public std::unary_function<Pair, Second>
{
  const Second& operator()(const Pair& x) const 
  {
    return x.second;
  }
};

// use this functor to iterate over items in a set
template <class Item>
struct SetSelect : public std::unary_function<Item, Item>
{
  const Item& operator()(const Item& x) const 
  {
    return x;
  }
};

template <class TreeType, class Y, class Which> struct TreeSelectFind 
: public std::unary_function<typename TreeType::value_type, bool>
{
  Y* mToFind;
  TreeSelectFind(Y* y) : mToFind(y) {}
  bool operator()(const typename TreeType::value_type& comp) const
  {
    return mToFind == Which()(comp);
  }
};

// TreeType is a set, map, or multimap, OutType is whatever type DynamicDLIIterator needs to return
// and Which is one of MapSelect1st, MapSelect2nd, SetSelect, or whatever you want
// for example
// typedef std::set<CubitEntity*> MySet;
// MySet the_set;
// DynamicTreeIterator<MySet, CubitEntity, SetSelect<MySet::value_type> >* iter = 
//        new DynamicTreeIterator<MySet, CubitEntity, SetSelect<MySet::value_type> >;
// iter->watch(the_set);

template < class TreeType, class OutType, class Which > class DynamicTreeIterator
  : public DynamicDLIIterator<OutType>
{
  public:
    DynamicTreeIterator() : mTree(0) {}
    ~DynamicTreeIterator(){}

    void watch(TreeType& the_tree)
    {
      mTree = &the_tree;
      mIter = mTree->begin();
    }

    virtual void reset()
    {
      mIter = mTree->begin();
    }

    virtual void step(int n = 1)
    {
      if(mTree->empty())
        return;
      while(n--)
      {
        if(++mIter == mTree->end())
          mIter = mTree->begin();
      }
    }

    virtual unsigned int size() const
    {
      return mTree->size();
    }
    virtual OutType* get() const
    {
      if(mTree->empty())
        return NULL;
      return Which()(*mIter);
    }
    virtual OutType* next(int n = 1) const
    {
      if(mTree->empty())
        return NULL;
      typename TreeType::iterator iter = mIter;
      while(n--)
      {
        if(++iter == mTree->end())
          iter = mTree->begin();
      }
      return Which()(*iter);
    }
    virtual OutType* get_and_step(int n = 1)
    {
      OutType* tmp = get();
      step(n);
      return tmp;
    }

    virtual CubitBoolean move_to(OutType* item)
    {
      typename TreeType::iterator iter = find_item(item);
      if(iter != mTree->end())
      {
        mIter = iter;
        return CUBIT_TRUE;
      }
      return CUBIT_FALSE;
    }

    virtual CubitBoolean is_in_list(OutType* item) const
    {
      return find_item(item) != mTree->end();
    }
    virtual CubitBoolean is_at_end() const
    {
      if(mTree->empty())
        return CUBIT_TRUE;
      typename TreeType::iterator iter = mIter;
      ++iter;
      return iter == mTree->end();
    }
    virtual CubitBoolean is_at_beginning() const
    {
      return mTree->begin() == mIter;
    }
    virtual DynamicDLIIterator<OutType>* clone() const
    {
      DynamicTreeIterator<TreeType, OutType, Which>* new_iterator = new DynamicTreeIterator<TreeType, OutType, Which>;
      *new_iterator = *this;
      return new_iterator;
    }

  private:
    typename TreeType::iterator find_item(OutType* item) const
    {
      typename TreeType::iterator iter = std::find_if(mIter, mTree->end(), 
                TreeSelectFind<TreeType, OutType, Which>(item));
      if(iter == mTree->end())
      {
        iter = std::find_if(mTree->begin(), mIter, 
                    TreeSelectFind<TreeType, OutType, Which>(item));
        if(iter == mIter)
          iter = mTree->end();
      }
      return iter;
    }
    TreeType* mTree;
    typename TreeType::iterator mIter;
};

#endif // DYNAMIC_TREE_ITERATOR_HPP

