//- Class: StubProgressTool
//- Description: default progress tool implementation that does nothing
//- Created: 04/29/02 by Byron Hanks
//- Checked By:
//- Version:

#ifndef STUBPROGRESSTOOL_HPP
#define STUBPROGRESSTOOL_HPP

#include "CubitDefines.h"
#include "ProgressTool.hpp"
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT StubProgressTool : public ProgressTool
{
protected:

public:
  StubProgressTool();
  virtual ~StubProgressTool();

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
};

#endif

//EOF

