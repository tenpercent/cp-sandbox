//to call floating-point check from the Fortran routine
#include <math.h>
#ifdef _MSC_VER
#	include <windows.h> // MSVC Workaround
#	define isnan(x) ((x) != (x))
#	define isinf(x) isnan(x-x) 
#endif
#include "Ani3D_FC.h"

void fpcheck( double *a, int *flag ) 
{ 
   *flag = isnan(*a) && !isinf(*a); 
}

