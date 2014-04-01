//#include "string.h"
#include "TextProgressTool.hpp"
#include "CubitMessage.hpp"

//=============================================================================
// Description: constructor
// Notes:  This will initialize the progress bar and display it
//         The difference between nLower and nUpper will be the number of
//         integral steps the progress bar will before the process is complete
// Author: sjowen
// Date: 11/7/01
//=============================================================================
TextProgressTool::TextProgressTool()
{
  nTotal = 0;
  nCurrent = 0;
  nShown = 0;
}

TextProgressTool::~TextProgressTool()
{
}

void TextProgressTool::start(int nLower, int nUpper,
                           const char* title           /* = NULL */,
                           const char* info_string     /* = NULL      */,
                           CubitBoolean ,
                           CubitBoolean )
{
  const char* progress_title = title ? title : "Progress";

  //bSmooth; bHasCancelButton;
  PRINT_INFO("%s\n",progress_title);
  if (info_string != NULL)
  {
    PRINT_INFO("%s\n",info_string);
  }
  PRINT_INFO("0   |    |    |    |    50   |    |    |    |  100\n");

  nTotal = nUpper - nLower;
  nCurrent = 0;
  nShown = 0;
}

//=============================================================================
// Description: end
// Notes:  closes/finishes the prgress bar window
// Author: sjowen
// Date: 11/7/01
//=============================================================================
void TextProgressTool::end()
{
  PRINT_INFO("\n");
}

//=============================================================================
// Description: makes one integer step of the progress bar
// Notes:
// Author: sjowen
// Date: 11/7/01
//=============================================================================
void TextProgressTool::step()
{
  nCurrent++;
  int n = 50 * nCurrent / nTotal;
  for( ; nShown < n; nShown++ )
    PRINT_INFO("*");
}

//=============================================================================
// Description: moves the progress bar to the specified percentage
// Notes:  only moves if percentage is greater the the current step
//         pcnt should be between 0 and 1
// Author: sjowen
// Date: 11/7/01
//=============================================================================
void TextProgressTool::percent( double pcnt )
{
  int ncur = (int)(pcnt * (double)nTotal + 0.5);
  if(ncur > nTotal) 
    ncur = nTotal;
  if (ncur > nCurrent)
  {
    int ii;
    for (ii=nCurrent; ii<ncur; ii++)
    {
      step();
    }
    nCurrent = ncur;
  }
}
