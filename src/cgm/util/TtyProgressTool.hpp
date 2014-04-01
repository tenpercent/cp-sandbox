//-------------------------------------------------------------------------
// Filename      : TtyProgressTool.hpp
//
// Purpose       : Progress bar for use in text terminals (ttys)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/16/03
//-------------------------------------------------------------------------
#ifndef TTY_PROGRESS_TOOL_HPP
#define TTY_PROGRESS_TOOL_HPP

#include "ProgressTool.hpp"
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT TtyProgressTool : public ProgressTool
{
public:

  inline TtyProgressTool();
  virtual ~TtyProgressTool();
  
  virtual void start( int lower, int upper, 
                      const char* title = 0,
                      const char* info = 0,
                      CubitBoolean smooth = CUBIT_TRUE,
                      CubitBoolean cancelButton = CUBIT_FALSE );
  
  virtual void end();
  
  virtual void step();
  
  virtual void percent( double p );
  
  void clear_all();  // clear all object data (inverse of start(..)).
  
private:

  bool update();     // check for change in tty width and update if necessary
  void display(int step, int pct);    // display progress line

  int myRange;       // number of steps to reach 100% 
  int currentVal;    // current step count
  int currPercent;   // current percent
  int prevWidth;     // result of most recent query for terminal width
  
  char* lineBuffer;   // buffer holding image of entire line of output
  char* barStart;     // start of progress bar (pointer into lineBuffer)
  char* bufferEnd;    // one past last buffer character

  char* firstMessage; // copy of first text message passed into start(..)
  char* secondMessage;// copy of second text message passed into start(..)
};

inline TtyProgressTool::TtyProgressTool()
  : prevWidth(0), lineBuffer(0), firstMessage(0), secondMessage(0) {}

#endif

  
  
