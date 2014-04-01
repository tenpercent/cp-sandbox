//-------------------------------------------------------------------------
// Filename      : TtyProgressTool.cpp
//
// Purpose       : Progress bar for use in text terminals (ttys)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/16/03
//-------------------------------------------------------------------------
#ifndef NT
#  include <unistd.h>
#endif
#include <string.h>
#include "AppUtil.hpp"
#include "TtyProgressTool.hpp"
#include "CubitUtil.hpp"

  // Enable/disable display of progress bar and/or numeric percent.
const bool DISPLAY_PROGRESS_BAR = true;
const bool DISPLAY_NUM_PERCENT  = true;

  // Min size for progress bar.  Will not be shown if available
  // space is less than this value.
const int MIN_PROGRESS_BAR_WIDTH = 10;
  // Maximum size for progress bar.  If 0, limited only by available
  // space.
const int MAX_PROGRESS_BAR_WIDTH =  0;

  // Characters composing progress bar.
const char PROGRESS_BAR_START   = '|';
const char PROGRESS_BAR_END     = '|';
const char PROGRESS_BAR_FILLED  = '=';
const char PROGRESS_BAR_CURRENT = '>';
const char PROGRESS_BAR_EMPTY   = ' ';


/* buffer layout:
 +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
 |\r|m |e |s |s |a |g |e |: |  || |= |= |= |= |= |> |  |  || |9 |9 |9 |% |
 +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
  ^                             ^                                         ^
  |                             |                                         |
  +--- lineBuffer               +--- barStart                bufferEnd ---+
*/


TtyProgressTool::~TtyProgressTool()
{ 
  clear_all();
}

void TtyProgressTool::clear_all()
{
    // free allocated buffers 
    // (messages were created with strdup (CubitUtil::util_strdup),
    // which uses malloc so use free() (CubitUtil::util_strdup_free())).
  delete [] lineBuffer; 
  if (firstMessage) 
    CubitUtil::util_strdup_free(firstMessage);
  if (secondMessage)
    CubitUtil::util_strdup_free(secondMessage);
  
    // zero data so calls to step or percent fail.
  myRange = 0;
  prevWidth = 0;
  lineBuffer = barStart = bufferEnd = 0;
  firstMessage = secondMessage = 0;
}


bool TtyProgressTool::update()
{
    // Minimum/maximum space for progress bar.  (Add 2 to constants
    // for the leading and trailing '|' character.) 
  const int MIN_PROGRESS_WIDTH = MIN_PROGRESS_BAR_WIDTH + 2;
  const int MAX_PROGRESS_WIDTH = MAX_PROGRESS_BAR_WIDTH + 2;
  
    // Get terminal width.  Assume 80 if query fails.
  int width, height;
  if ( !AppUtil::instance()->get_terminal_size(height, width) )
    width = 80;

#ifdef NT
    // Windows cmd.exe wraps the line when the last character
    // of the line is output, rather than the first character
    // after the end of the line.  Stop one short so to avoid
    // newlines.
  width--;
#endif

    // If terminal width hasn't changed then no update is
    // required.  return.
  if ( width && width == prevWidth )
    return true;
  prevWidth = width;
  
    // If width isn't at least 6 don't try to display anything.
    // Clear data and return.
  if ( width < 4 )
  {
    delete [] lineBuffer;
    lineBuffer = barStart = bufferEnd = 0;
    return false;
  }

    // Re-create lineBuffer for the new terminal width if necessary.
  if ( !lineBuffer || (bufferEnd - lineBuffer - 1 < width) )
  {
    delete [] lineBuffer;
#ifdef NT
    lineBuffer = new char [width+2];
    lineBuffer[width+1] = '\r';
#else
    lineBuffer = new char [width+1];
#endif
    if ( !lineBuffer )
      return false;
  }
  
  bufferEnd = lineBuffer + width + 1;    // one past end of buffer
  barStart = 0;              // start of progress bar (null if no progress bar)
  char* bar_end = 0;         // end of progress bar

    // Current start and end (fill lineBuffer until start == end)
  char* start = lineBuffer;
  char* end = bufferEnd;
  
    // Use \r to move cursor to beginning of line without advancing a line.
  *(start++) = '\r'; // leading \r
  
    // If displaying numeric percent on end of line, allocate 4 spaces
    // at end of line and place '%' char in the fourth one.
  if (DISPLAY_NUM_PERCENT && (end - start > 3))
  {
    end[-1] = '%';
    end -= 4;
  }

    // Copy first message into line
  if (firstMessage && (end - start > 2))
  {
    int len = strlen(firstMessage);
      // need to truncate message?
    if (len + 2 > end - start)
      len = end - start - 2;
    memcpy(start, firstMessage, len);
    start += len;
      // add trailing ": "
    *(start++) = ':';
    *(start++) = ' ';
  }
  
    // Allocate space for minimum progress bar if sufficient space
  if (DISPLAY_PROGRESS_BAR && (end - start > MIN_PROGRESS_WIDTH))
  {
    bar_end = end;
    end -= MIN_PROGRESS_WIDTH;
  }
  
    // Copy second message into line buffer
  if (secondMessage && (end - start > 2))
  {
    int len = strlen(secondMessage);
      // need to truncate message?
    if (len + 2 > end - start) 
      len = end - start - 2;
    memcpy(start, secondMessage, len);
    start += len;
      // add trailing ": "
    *(start++) = ':';
    *(start++) = ' ';
  }
  
    // Finish progress bar setup (if progress bar is displayed at all)
  if (bar_end)
  {
      // Grow progress bar up to MAX_PROGRESS_WIDTH (or all the 
      // remaining space if MAX_PROGRESS_BAR_WIDTH is unset)
    int diff = bar_end - start - MAX_PROGRESS_WIDTH;
    end = start;
    if (MAX_PROGRESS_BAR_WIDTH && diff >  0)
      end += diff;
    
      // Put in bounding '|' chars for progress bar and
      // initialize barStart to point *after* the leading '|'
    barStart = end;
    barStart[0] = PROGRESS_BAR_START;
    bar_end[-1] = PROGRESS_BAR_END;
    barStart++;
  }
  
    // fill any remaining space with spaces
  if (start < end)
    memset( start, ' ', end - start );
      
  return true;
}

void TtyProgressTool::start( int lower, int upper, 
                             const char* s1, const char* s2, 
                             CubitBoolean, CubitBoolean )
{
    // get rid of any old state if someone forgot to call end().
  clear_all();
  
    // store passed range
  myRange = upper - lower;
  currentVal = currPercent = 0;
  prevWidth = 0;
  
    // copy passed messages
  const char* empty = "";
  if (s1)
    firstMessage = CubitUtil::util_strdup((char*)s1);
  else
    firstMessage = CubitUtil::util_strdup((char*)empty);
  if (s2)
    secondMessage = CubitUtil::util_strdup((char*)s2);

    // generate line buffer
  display(0,0);
}

void TtyProgressTool::step()
{
  int new_count = currentVal + 1;
  display( new_count, 100 * new_count / myRange );
}

void TtyProgressTool::percent( double p )
{
  display( (int)(p * myRange), (int)(100 * p) );
}

void TtyProgressTool::end()
{
    // clear terminal line
  if (lineBuffer)
  {
    memset( lineBuffer + 1, ' ', bufferEnd - lineBuffer - 1 );
    fwrite( lineBuffer, bufferEnd - lineBuffer, 1, stdout );
    putchar('\r');
    fflush( stdout );
  }

    // clear internal data  
  clear_all();
}

void TtyProgressTool::display( int new_count, int new_percent )
{
  int prev_width = prevWidth;
  bool do_output = false;
  
    // check for change in terminal width and update if necessary.
  if (!update())
    return;
  
    // if tty size changed or percent changed, need to update  
  if (prev_width != prevWidth || new_percent != currPercent)
  {
    currPercent = new_percent;
    do_output = true;
  }

    // calculate start of numeric percent in lineBuffer
  char* pctStart = bufferEnd;
  if (DISPLAY_NUM_PERCENT)
    pctStart -= 4; // start of numeric percent
  
    // check if the progress bar output is changing, and if
    // so set do_output to true
  size_t bar_width = 0, bar_drawn = 0, old_drawn;
  if (barStart)
  {
    bar_width = pctStart - barStart - 1;
    bar_drawn = new_count >= myRange ? bar_width :
                new_count <= 0       ? 0         :
                bar_width * new_count / myRange;
    old_drawn = currentVal >= myRange ? bar_width :
                currentVal <= 0       ? 0         :
                bar_width * currentVal / myRange;
    if (bar_drawn != old_drawn)
      do_output = true;
  }
  currentVal = new_count;
  
  if (DISPLAY_NUM_PERCENT && do_output)
  {
      // normalize numeric percent to values that can be displayed
    int pct = currPercent;
    if (pct < 0)
      pct = 0;
    else if(pct > 999)
      pct = 999;

      // display numeric percent in last three chars of line
    pctStart[2] = (char)('0' + pct % 10);
    pctStart[1] = pctStart[0] = ' ';
    for ( int i = 1; i >= 0 && pct > 9; i-- )
    {
      pct /= 10;
      pctStart[i] = (char)('0' + pct % 10);
    }
  }
  
    // if progress bar is to be displayed, update it.
  if (barStart && do_output)
  {
      // write '=' for filled section of bar
    if (bar_drawn > 1)
      memset(barStart, PROGRESS_BAR_FILLED, bar_drawn - 1);
      // write '>' between filled and unfilled sections of bar
    if (bar_drawn > 0)
      barStart[bar_drawn - 1] = PROGRESS_BAR_CURRENT;
      // write spaces in unfilled part of bar
    if (bar_drawn < bar_width)
      memset(barStart + bar_drawn, PROGRESS_BAR_EMPTY, bar_width - bar_drawn);
  }
  
    // write the buffer to the terminal
  if (do_output)
  {
    fwrite( lineBuffer, bufferEnd - lineBuffer, 1, stdout );
    fflush( stdout );
  }
}
