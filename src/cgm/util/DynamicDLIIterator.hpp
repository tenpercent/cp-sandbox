#ifndef DYNAMIC_DLI_ITERATOR_HPP
#define DYNAMIC_DLI_ITERATOR_HPP

#include "DLIList.hpp"

template <class ListType, class OutType> class CastingDLIIterator
{
public:
  typedef ListType* internal_type;
  typedef OutType* output_type;
  
  CastingDLIIterator() 
    {}
  virtual ~CastingDLIIterator()
    {}
  
  void watch(const DLIList<internal_type> *dl_i_list)
    { mInternal.watch(dl_i_list); }
  
  void reset()
    { mInternal.reset(); }
  
  void step( int n = 1 ) 
    { mInternal.step(n); }
  
  int size() const
    { return mInternal.size(); }
  
  output_type get() const
    { return dynamic_cast<output_type>(mInternal.next(0)); }
  
  output_type next(int n = 1) const
    { return dynamic_cast<output_type>(mInternal.next(n)); }
  
  output_type get_and_step(int n = 1) 
    { return dynamic_cast<output_type>(mInternal.get_and_step(n)); }
  
  CubitBoolean move_to(output_type item)
    { return mInternal.move_to(dynamic_cast<internal_type>(item)); }
  
  CubitBoolean is_in_list(output_type item) const
    { return mInternal.is_in_list(dynamic_cast<internal_type>(item)); }

  CubitBoolean is_at_end() const
    { return mInternal.is_at_end(); }
  
  CubitBoolean is_at_beginning() const
    { return mInternal.is_at_beginning(); }
  
private:
  DLIListIterator<internal_type> mInternal;
};

template <class X> class DynamicDLIIterator
{
public:
  DynamicDLIIterator() 
    {}
  virtual ~DynamicDLIIterator()
    {}
  
  virtual void reset() = 0;
  virtual void step(int n = 1) = 0;
  virtual unsigned int size() const = 0;
  virtual X* get() const = 0;
  virtual X* next(int n = 1) const = 0;
  virtual X* get_and_step(int n = 1) = 0; 
  virtual CubitBoolean move_to(X* item) = 0;
  virtual CubitBoolean is_in_list(X* item) const = 0;
  virtual CubitBoolean is_at_end() const = 0;
  virtual CubitBoolean is_at_beginning() const = 0;
  virtual DynamicDLIIterator<X>* clone() const = 0;
};

template<class ListType, class OutType> class DynamicDLIIteratorImpl :
  public DynamicDLIIterator<OutType>,
  public CastingDLIIterator<ListType, OutType>
{
public:
    DynamicDLIIteratorImpl() 
    {}
  virtual ~DynamicDLIIteratorImpl()
    {}
  
  virtual void reset()
    { CastingDLIIterator<ListType, OutType>::reset(); }
  
  virtual void step(int n = 1)
    { CastingDLIIterator<ListType, OutType>::step(n); }
  
  virtual unsigned int size() const
    { return CastingDLIIterator<ListType, OutType>::size(); }
  
  virtual OutType* get() const
    { return CastingDLIIterator<ListType, OutType>::get(); }
  
  virtual OutType* next(int n = 1) const
    { return CastingDLIIterator<ListType, OutType>::next(n); }
  
  virtual OutType* get_and_step(int n = 1)
    { return CastingDLIIterator<ListType, OutType>::get_and_step(n); }
  
  virtual CubitBoolean move_to(OutType* item)
    { return CastingDLIIterator<ListType, OutType>::move_to(item); }
  
  virtual CubitBoolean is_in_list(OutType* item) const
    { return CastingDLIIterator<ListType, OutType>::is_in_list(item); }
  
  virtual CubitBoolean is_at_end() const
    { return CastingDLIIterator<ListType, OutType>::is_at_end(); }
  
  virtual CubitBoolean is_at_beginning() const
    { return CastingDLIIterator<ListType, OutType>::is_at_beginning(); }

  virtual DynamicDLIIterator<OutType>* clone() const
    {
      DynamicDLIIteratorImpl<ListType, OutType>* rv = new DynamicDLIIteratorImpl<ListType, OutType>;
      *rv = *this;
      return rv;
    }
};

template <class ListType, class OutType> class CastingDynamicDLIIterator :
  public DynamicDLIIterator<OutType>
{
public:
  typedef ListType* internal_type;
  typedef OutType* output_type;

  CastingDynamicDLIIterator(DynamicDLIIterator<ListType>* internal) 
          : mInternal(internal)
    {}
  virtual ~CastingDynamicDLIIterator()
    { delete mInternal; }
  
  void reset()
    { mInternal->reset(); }
  
  void step( int n = 1 ) 
    { mInternal->step(n); }
  
  unsigned int size() const
    { return mInternal->size(); }
  
  output_type get() const
    { return dynamic_cast<output_type>(mInternal->next(0)); }
  
  output_type next(int n = 1) const
    { return dynamic_cast<output_type>(mInternal->next(n)); }
  
  output_type get_and_step(int n = 1) 
    { return dynamic_cast<output_type>(mInternal->get_and_step(n)); }
  
  CubitBoolean move_to(output_type item)
    { return mInternal->move_to(dynamic_cast<internal_type>(item)); }
  
  CubitBoolean is_in_list(output_type item) const
    { return mInternal->is_in_list(dynamic_cast<internal_type>(item)); }

  CubitBoolean is_at_end() const
    { return mInternal->is_at_end(); }
  
  CubitBoolean is_at_beginning() const
    { return mInternal->is_at_beginning(); }

  DynamicDLIIterator<OutType>* clone() const
    { return new CastingDynamicDLIIterator(mInternal->clone()); }
  
  
private:
  DynamicDLIIterator<ListType>* mInternal;
};


#endif

