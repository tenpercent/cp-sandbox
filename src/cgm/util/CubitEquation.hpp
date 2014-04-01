#ifndef CUBITEQUATION_HPP
#define CUBITEQUATION_HPP

class CubitEquation
{
public:
  CubitEquation()
    {}
  virtual ~CubitEquation()
    {}

  virtual double evaluate(double x, double y, double z, int n) = 0;
  virtual CubitEquation* clone() const = 0;
};

#endif
