//class DoubleListItem
// for creating a (sorted) list of pointers to doubles.
// samitch

#ifndef DOUBLE_LIST_ITEM_HPP
#define DOUBLE_LIST_ITEM_HPP

class DoubleListItem {
public:
  DoubleListItem (double set_value) { myValue = set_value; }
  double myValue;
  double double_value() {return myValue;}
};

#endif

