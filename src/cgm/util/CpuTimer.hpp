//- Class: CpuTimer

#if !defined(CPU_TIMER)
#define CPU_TIMER
#include <sys/types.h>
#include "CubitUtilConfigure.h"


class CUBIT_UTIL_EXPORT CpuTimer {
public:
  CpuTimer();			//- initialise to current system time
  double cpu_secs();		//- return CPU time in seconds since last
                                //- to cpu_secs();
  double elapsed();             //- return CPU time in seconds since 'birth'
  
private:
  time_t cpu;			//- cpu time in 1/HZ units since last call
                                //- to cpu_secs()
  time_t cpuInitial;             //- cpu time in 1/HZ units since construction.

 // Added by Cat for NT port
  #ifdef NT
  void  nt_times(struct tms *);
  #endif


};

#endif

