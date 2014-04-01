//-------------------------------------------------------------------------
// Filename      : AppUtil.cc
//
// Purpose       : This file represents the Cubit application itself.
//
// Special Notes :
//
// Creator       : Darryl Melander
//
// Date          : 06/08/98
//
// Owner         : Darryl Melander
//-------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#ifndef NT
#  include <unistd.h>
#  include <termios.h>
#  include <sys/ioctl.h>
#endif
#include <ctype.h>
#include <time.h>
//#include <new.h>

#include "AppUtil.hpp"
#include "CubitMessage.hpp"
#include "CubitString.hpp"
#include "CubitObserver.hpp"
#include "SettingHandler.hpp"
#include "StubProgressTool.hpp"

// Different platforms follow different conventions
#ifndef NT
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif
#ifdef SOLARIS
#include <procfs.h>
#include <sys/syscall.h>
extern "C" int getrusage(int, struct rusage *);
#ifndef RUSAGE_SELF
#include </usr/ucbinclude/sys/rusage.h>
#endif
#endif

#ifdef HP
#include <sys/syscall.h>
#if 0
#if OS_VERSION != 1020
#define getrusage(a, b)  syscall(SYS_GETRUSAGE, a, b)
#endif
#endif
#endif

AppUtil* AppUtil::instance_ = NULL;

#ifdef HP
#ifndef UNIX_PORT
extern "C" int Error;
#else
int Error = 0;
#endif //unix_port
#endif
#ifdef NT
#  include <direct.h>
#  include <windows.h>
//#  define strdup _strdup
#endif
#include "CubitUtil.hpp"



extern "C" void cubit_update_terminal_size(int)
{
#if defined(NT)

  // http://msdn.microsoft.com/library/default.asp?url=/library/en-us/dllproc/base/getconsolescreenbufferinfo.asp
  // http://msdn.microsoft.com/library/default.asp?url=/library/en-us/dllproc/base/scrolling_a_screen_buffer_s_window.asp
  HANDLE h_stdout = GetStdHandle(STD_OUTPUT_HANDLE);
  CONSOLE_SCREEN_BUFFER_INFO scr_info;
  if (GetConsoleScreenBufferInfo(h_stdout, &scr_info))
  {
	  AppUtil::instance()->set_terminal_size(
        // dwSize.Y is will return the height of the scroll-back buffer.
        // srWindow contains the position of the visible rect of the
        // window in the buffer.  That is the height of the window.
      scr_info.srWindow.Bottom - scr_info.srWindow.Top,
        // use width of scroll buffer rather than size of window
        // for number of columns.  The window size is useless as
        // caller would need to know the offset in the X direction
        // to make use of it and we aren't returning that.  Besides,
        // the user presumably set the scroll buffer to this width
        // because that's the line length (s)he wants.
      scr_info.dwSize.X
    );
  }

#elif defined(TIOCGWINSZ)

    // On UNIX platforms, register this function as a handler for
    // SIGWINCH.  The system will then call this function any time
    // the user changes the size of the terminal window.
#ifdef SIGWINCH
  signal( SIGWINCH, &cubit_update_terminal_size );
#endif

  const int file = fileno(stdout);
  struct winsize size;
  if( ioctl( file, TIOCGWINSZ, (char*)&size ) == 0 )
  {
    AppUtil::instance()->set_terminal_size( size.ws_row, size.ws_col );
  }

#endif
}


CubitBoolean AppUtil::catching_sigint_ = CUBIT_FALSE;

// Returns the singleton instance of the app object.
// If it is being created for the first time, it
// does application initialization.
AppUtil* AppUtil::instance()
{
   if (instance_ == NULL)
   {
      instance_ = new AppUtil();

      if ( !instance_ )
      {
         PRINT_ERROR(" *** Unable to initialize application ***\n");
         exit(1);
      }

        // Initialize interrupted flag to false.
      cubit_intr = CUBIT_FALSE;

        // initialize the static observer list
      CubitObserver::init_static_observers();

      cubit_update_terminal_size(0);

      initialize_settings();
   }
   return instance_;
}

AppUtil::AppUtil()
{
    term_width = 0;
    term_height = 0;
    mAppStarted = CUBIT_FALSE;

      // reset catching_sigint_ here, even though it's static, 'cuz it used
      // to be reset to false on construction before going static

    catching_sigint_ = CUBIT_FALSE;

    // create a default progress tool that does nothing
    // an app can replace this with an app specific progress tool
    mProgressTool = new StubProgressTool();
}

AppUtil::~AppUtil()
{
  instance_ = NULL;
  catch_interrupt( CUBIT_FALSE );
  if (mProgressTool)
      delete mProgressTool;
  mProgressTool = NULL;
}

#ifdef CAT
// Returns how many days until the app expires.  A negative
// number inicates it is already expired.  If there is
// an extension, has_extension is set to CUBIT_TRUE.
int AppUtil::days_to_expire(int year, int month, int day ) const
{
     // Setup variables and constants
   struct tm expiration_time;
   time_t raw_expire;
   time_t raw_current;
   int days;
   const int SECONDS_PER_DAY = 24 * 60 * 60;

     // Get the current time, raw
   time(&raw_current);

     // Setup a dissected expiration date
   expiration_time.tm_year  = year - 1900;
   expiration_time.tm_mon   = month - 1;
   expiration_time.tm_mday  = day;
   expiration_time.tm_hour  = 0;
   expiration_time.tm_min   = 0;
   expiration_time.tm_sec   = 0;
   expiration_time.tm_isdst = 0;

     // Convert expiration date from dissected to raw
   raw_expire = mktime(&expiration_time);
     // Get difference between dates
   days = (int)difftime(raw_expire, raw_current) / SECONDS_PER_DAY;
   //days = INT_MAX;//(int)difftime(raw_expire, raw_current) / SECONDS_PER_DAY;
//changed back by aga@cat|11/2/04

     // Return days left
   return days;
}
#endif

// Prints out information about the session, like
// execution time and memory problems.
void AppUtil::report_resource_usage() const
{

#ifndef JANUS
#ifndef NT
   struct rusage r_usage;
   float utime, stime;
   apputil_getrusage(r_usage);
   utime = (float)r_usage.ru_utime.tv_sec +
      ((float)r_usage.ru_utime.tv_usec/1.e6);
   stime = (float)r_usage.ru_stime.tv_sec +
      ((float)r_usage.ru_stime.tv_usec/1.e6);
   static int pagesize = getpagesize();
   PRINT_INFO("\nEnd of Execution\n");
   PRINT_INFO("User time             = %f\n", utime);
   PRINT_INFO("System time           = %f\n", stime);
   PRINT_INFO("Total time            = %f\n", utime+stime);
   PRINT_INFO("Max resident set size = %ld bytes\n", r_usage.ru_maxrss*pagesize);
   PRINT_INFO("Int resident set size = %ld\n", r_usage.ru_idrss);
   PRINT_INFO("Minor Page Faults     = %ld\n", r_usage.ru_minflt);
   PRINT_INFO("Major Page Faults     = %ld\n", r_usage.ru_majflt);
   PRINT_INFO("Swaps                 = %ld\n", r_usage.ru_nswap);
#endif// NT
#endif//JANUS
}


// Interrupt handling defintions.
volatile CubitBoolean cubit_intr = CUBIT_FALSE;

// Default handler for SIGINT provided by AppUtil.
extern "C" void sigint_handler(int)
{
#ifndef CUBIT_NO_SIGNAL
  if( signal( SIGINT, sigint_handler ) == SIG_ERR )
    PRINT_ERROR("Cannot continue to catch SIGINT!\n");
  cubit_intr = CUBIT_TRUE;
#endif
}

// Set signal handler for SIGINT.
void AppUtil::catch_interrupt( CubitBoolean yesno )
{
#ifndef CUBIT_NO_SIGNAL
  if( yesno )
  {
    if( signal( SIGINT, sigint_handler ) == SIG_ERR )
      PRINT_ERROR("Can't catch SIGINT!\n");
    else
      catching_sigint_ = CUBIT_TRUE;
  }
  else
  {
    if( signal( SIGINT, SIG_DFL ) == SIG_ERR )
      PRINT_ERROR("Can't reset SIGINT handler!\n");
    else
      catching_sigint_ = CUBIT_FALSE;
  }
#endif
}

void AppUtil::startup(int argc, char ** argv)
{

    if (mAppStarted)
        return;

    if ( argc > 0 ) {
      char* path = CubitUtil::util_strdup(argv[0]);

      char* end;
      for (end = path; *end; end++);

      for ( ; end > path; end--) {
        if ( *end == '/'
#ifdef NT
             || *end == '\\'
#endif
           ) {
          break;
        }
      }
      *end = '\0';
      end = path;

#ifndef NT
      bool local = *end != '/';
#else
      bool local = true;
      if (isalpha(*end) && (*(end+1) == ':'))
      {
        if (*(end+2) == '/' || *(end+2) == '\\')
          local = false;
        else
          end += 2;
      }
      else if ( ((*end == '/') || (*end == '\\')) && (*end == *(end+1)) )
        local = false;
#endif

      if (local)
      {
#ifndef NT
        char buffer[PATH_MAX];
        getcwd(buffer, PATH_MAX);
#else
        char buffer[_MAX_PATH];
        _getcwd(buffer, _MAX_PATH);
#endif
        cubitDir = buffer;
        cubitDir += "/";
        cubitDir += end;
      }
      else
        cubitDir = path;

      CubitUtil::util_strdup_free(path);
    }

    mAppStarted = CUBIT_TRUE;
    // add startup code here
    // nothing to do for now

    return;
}

int AppUtil::shutdown()
{
    // When the user is done, print out additional info if requested
   if (DEBUG_FLAG(3))
      report_resource_usage();

     // Return the number of errors in the session
   int ret_val = ( CubitMessage::instance()->error_count() );
   if ( ret_val > 0 )
   {
     PRINT_ERROR("Errors found during CUBIT session.\n");
   }

     // Delete the singletons
   delete CubitMessage::instance();

   CubitObserver::term_static_observers();

   SettingHandler::delete_instance();

   mAppStarted = CUBIT_FALSE;

   return ret_val;
}
void AppUtil::apputil_getrusage(struct rusage &r_usage) const
{
    // get the resource usage as defined by getrusage
#ifndef JANUS
#ifndef NT
   getrusage(RUSAGE_SELF, &r_usage);

   if (r_usage.ru_maxrss == 0) {
       // this machine doesn't return rss - try going to /proc
       // print the file name to open
#ifdef SOLARIS
     char buffer[120];
     strcpy(buffer, "/proc/self/psinfo");
     int file_ptr = -1;
     file_ptr = open(buffer, O_RDONLY);
     if (file_ptr < 0) return;
     struct psinfo myps_info;
     read(file_ptr, &myps_info, sizeof(myps_info));
     static int page_size = getpagesize();
     r_usage.ru_maxrss = myps_info.pr_rssize*1024/page_size;
     r_usage.ru_idrss = myps_info.pr_size*1024/page_size;
     close (file_ptr);
#endif // SOLARIS
#ifdef CUBIT_LINUX
     char file_str[4096], dum_str[4096];
     int file_ptr = -1, file_len;
     file_ptr = open("/proc/self/stat", O_RDONLY);
     file_len = read(file_ptr, file_str, sizeof(file_str)-1);
     if (file_len == 0) return;
     close(file_ptr);
     file_str[file_len] = '\0';
       // read the preceeding fields and the ones we really want...
     int dum_int;
     unsigned int dum_uint, vm_size, rss;
     static int page_size = getpagesize();
     int num_fields = sscanf(file_str,
                             "%d " // pid
                             "%s " // comm
                             "%c " // state
                             "%d %d %d %d %d " // ppid, pgrp, session, tty, tpgid
                             "%u %u %u %u %u " // flags, minflt, cminflt, majflt, cmajflt
                             "%d %d %d %d %d %d " // utime, stime, cutime, cstime, counter, priority
                             "%u %u " // timeout, itrealvalue
                             "%d " // starttime
                             "%u %u", // vsize, rss
                             &dum_int,
                             dum_str,
                             dum_str,
                             &dum_int, &dum_int, &dum_int, &dum_int, &dum_int,
                             &dum_uint, &dum_uint, &dum_uint, &dum_uint, &dum_uint,
                             &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, &dum_int,
                             &dum_uint, &dum_uint,
                             &dum_int,
                             &vm_size, &rss);
     if (num_fields == 24) {
       r_usage.ru_maxrss = rss/page_size;
       r_usage.ru_idrss = vm_size/page_size;
     }
#endif // CUBIT_LINUX
   }
#endif // NT
#endif // JANUS
}

bool AppUtil::get_terminal_size( int& rows, int& cols )
{
  cols = term_width;
  rows = term_height;
  return rows || cols;
}

void AppUtil::set_terminal_size( int rows, int cols )
{
  term_width = cols;
  term_height = rows;
}

ProgressTool *AppUtil::progress_tool()
{
  assert(mProgressTool != NULL);

  return mProgressTool;
}

void AppUtil::progress_tool(ProgressTool* pTool)
{
  assert(pTool != NULL);  // existence is assumed elsewhere in the code
  if (pTool)
  {
      if (mProgressTool)
      {
          delete mProgressTool;
      }
      mProgressTool = pTool;
  }
}

void AppUtil::initialize_settings()
{
  CubitMessage::initialize_settings();
}
