//- Class: TextProgressTool
//- Description: interface for displaying progress on a process
//- Created: 11/6/01
//- Checked By:
//- Version:

#ifndef TEXTPROGRESSTOOL_HPP
#define TEXTPROGRESSTOOL_HPP

#include "ProgressTool.hpp"
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT TextProgressTool : public ProgressTool
{
public:
  TextProgressTool();
  virtual ~TextProgressTool();

  virtual void start(int nLower, int nUpper,
               const char* title = NULL,
               const char* info_string = NULL,
               CubitBoolean bSmooth = CUBIT_TRUE,
               CubitBoolean bHasCancelButton = CUBIT_FALSE);

  virtual void end();

  virtual void step();
  // move the control bar one integer step based on nLower and nUpper

  virtual void percent( double pcnt );
  // move the control bar to a specified percent
  // pcnt should be between 0 and 1

private:
  int nTotal;
  int nCurrent;
  int nShown;
};

#endif

//EOF

