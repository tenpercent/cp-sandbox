#include "StubProgressTool.hpp"

StubProgressTool::StubProgressTool()
{
}

StubProgressTool::~StubProgressTool()
{
}

//=============================================================================
// Description: 
// Notes:  This will initialize the progress bar and display it
//         The difference between nLower and nUpper will be the number of
//         integral steps the progress bar will before the process is complete
// Author:bwhanks
// Date: 4/29/02
//=============================================================================
void
StubProgressTool::start(int , int ,
                           const char* ,
                           const char* ,
                           CubitBoolean ,
                           CubitBoolean )
{
}

//=============================================================================
// Description: end
// Notes:  closes/finishes the prgress bar window
// Author:bwhanks
// Date: 4/29/02
//=============================================================================
void
StubProgressTool::end()
{
}

//=============================================================================
// Description: makes one integer step of the progress bar
// Notes:
// Author:bwhanks
// Date: 4/29/02
//=============================================================================
void StubProgressTool::step()
{
}

//=============================================================================
// Description: moves the progress bar to the specified percentage
// Notes:  only moves if percentage is greater the the current step
//         pcnt should be between 0 and 1
// Author:bwhanks
// Date: 4/29/02
//=============================================================================
void StubProgressTool::percent( double )
{
}
