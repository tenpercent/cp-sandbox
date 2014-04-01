#include <stdio.h>
#include "CpuTimer.hpp"

#include "CubitMessage.hpp"

// Added by Cat for NT port
#ifndef NT
  #include <sys/param.h>
  #include <sys/times.h>
#endif

#include <time.h>	
#ifdef SOLARIS
#include <unistd.h>
#endif

#ifndef NT
#ifndef HZ
#ifdef CLK_TCK
#define HZ CLK_TCK
#else
#define HZ 60
#endif
#endif
#else
#define HZ CLOCKS_PER_SEC
#endif

// Added by Cat for NT port
#ifdef NT
struct tms {
	   clock_t tms_utime;		/* user time */
	   clock_t tms_stime;		/* system time */
	   clock_t tms_cutime;		/* user time, children */
	   clock_t tms_cstime;		/* system time, children */
      };
#endif


CpuTimer::CpuTimer()
{
  tms current;

  // Added by Cat for NT port
#ifdef NT
  nt_times(&current);    
#else
  times( &current );
#endif
  
  // Store current value of cpu time
  cpu = current.tms_utime +
        current.tms_stime +
        current.tms_cutime +
        current.tms_cstime;
  cpuInitial = cpu;
}

// Return time values

double
CpuTimer::cpu_secs()
{
  tms current;

  // Added by Cat for NT port
#ifdef NT
  nt_times(&current);    
#else
  times( &current );
#endif
  
  time_t cpu_now = current.tms_utime +
    current.tms_stime +
    current.tms_cutime +
    current.tms_cstime;
  time_t delta = cpu_now - cpu;
  cpu   = cpu_now;
#ifdef NT
  if( delta == 0 )
     delta = 1;
#endif
  return (double) delta / HZ;
}

double
CpuTimer::elapsed()
{
  tms current;

// Added by Cat for NT port
#ifdef NT
  nt_times(&current);    
#else
  times( &current );
#endif
// Store totals
  
  time_t cpu_now = current.tms_utime +
    current.tms_stime +
    current.tms_cutime +
    current.tms_cstime;
  time_t elapsed = cpu_now - cpuInitial;
  return (double) elapsed / HZ;
}


// Added by Cat for NT port
#ifdef NT

void CpuTimer::nt_times(tms *sys)
{
	sys->tms_utime =0;	
	sys->tms_stime = clock();		
	sys->tms_cutime = 0;		
	sys->tms_cstime = 0;
        //PRINT_DEBUG( CUBIT_DEBUG_3,"clock returned %d\n", sys->tms_stime );
}
#endif
