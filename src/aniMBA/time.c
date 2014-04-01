#ifndef _MSC_VER
#	include <unistd.h>             /* sysconf(_SC_CLK_TCK) */
#	include <sys/times.h>          /* struct tms, times() */
#	define SECS_PER_CLOCK (1./sysconf(_SC_CLK_TCK))
#else
#	include <time.h>
#	define SECS_PER_CLOCK  (1./CLOCKS_PER_SEC)
#endif

#include "Ani3D_FC.h"



//void seconds( double *usecs, double *ssecs ) 
void seconds( double *usecs ) 
{
#  ifndef _MSC_VER
   struct tms t;
#  endif
   clock_t utime;
   clock_t stime;

#  ifdef _MSC_VER
   utime = clock();
   stime = clock();
#  else
   times(&t);
   utime = t.tms_utime;
   stime = t.tms_stime;
#  endif

   *usecs = utime * SECS_PER_CLOCK;
   // *ssecs = stime * SECS_PER_CLOCK;
}


