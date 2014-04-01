//- Class: ProgressTool
//- Description: interface for displaying progress on a process
//- Created: 11/6/01
//- Checked By:
//- Version:

#ifndef PROGRESSTOOL_HPP
#define PROGRESSTOOL_HPP

#include "CubitDefines.h"

class ProgressTool
{
public:
  virtual void start(int nLower, int nUpper,
               const char* title = NULL,
               const char* info_string = NULL,
               CubitBoolean bSmooth = CUBIT_TRUE,
               CubitBoolean bHasCancelButton = CUBIT_FALSE) = 0;

  virtual void end() = 0;

  virtual void step() = 0;
  // move the control bar one integer step based on nLower and nUpper

  virtual void percent( double pcnt ) = 0;
  // move the control bar to a specified percent
  // pcnt should be between 0 and 1
};

#endif

//EOF

