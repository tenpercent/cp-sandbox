
#include <fstream>
#include <iomanip>

#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitUtil.hpp"
#include <assert.h>
#include <string.h>
#include <vector>

#include "SettingHandler.hpp"

#ifdef NT
#define vsnprintf _vsnprintf
//#define strdup _strdup
#endif

int CubitMessage::errorCount = 0;
int CubitMessage::warningCount = 0;
//
// Message type usage:
// PRINT_ERROR:    Message to tell user why the task did not complete.
// PRINT_WARNING:  Message to tell user why completed task may not be what
//                 was requested.
// PRINT_INFO:     Message to tell user about status and progress.
// PRINT_DEBUG:    Message to developer tied to a global debug flag.
// DIAGNOSTIC:     Message to developer.

CubitMessage* CubitMessage::instance_ = NULL;
CubitMessageHandler* CubitMessage::mHandler = NULL;
int CubitMessage::infoFlag = CUBIT_TRUE;
int CubitMessage::diagnosticFlag = CUBIT_FALSE;
int CubitMessage::warningFlag = CUBIT_TRUE;
std::ofstream* CubitMessage::loggingStream = NULL;
CubitString* CubitMessage::loggingFile = NULL;
std::ofstream* CubitMessage::loggingErrorStream = NULL;
CubitString* CubitMessage::loggingErrorFile = NULL;

CubitMessage* CubitMessage::instance()
{
  if (!instance_)
  {
    instance_ = new CubitMessage;
    if (!instance_)
    {
	  std::cerr << " *** Unable to instantiate message object ***" << std::endl;
      exit(1);
    }
  }
  return instance_;
}

CubitMessage::CubitMessage()
{
  infoFlag      = CUBIT_TRUE;
  warningFlag   = CUBIT_TRUE;
  diagnosticFlag= CUBIT_FALSE;
  loggingStream = NULL;
  loggingFile   = NULL;
  loggingErrorStream = NULL;
  loggingErrorFile   = NULL;
  currentDebugFlag = CUBIT_DEBUG_1;

    // Initialize the debugFlag array
  static MessageFlag staticDebugFlag[] =
  {
    MessageFlag( 0, "UNUSED"),
    MessageFlag( 1, "Debug Graphics toggle for some debug options."),
    MessageFlag( 2, "Whisker weaving information"),
    MessageFlag( 3, "Timing information for 3D Meshing routines."),
    MessageFlag( 4, "Graphics Debugging (DrawingTool)"),
    MessageFlag( 5, "FastQ debugging"),
    MessageFlag( 6, "Submapping graphics debugging"),
    MessageFlag( 7, "Knife progress whisker weaving information"),
    MessageFlag( 8, "Mapping Face debug / Linear Programing debug "),
    MessageFlag( 9, "Paver Debugging"),
    MessageFlag(10, "WW: removed hex seam flag"),
    MessageFlag(11, "Nodeset Associativity debugging"),
    MessageFlag(12, "Fastq activity"),
    MessageFlag(13, "Mesh entities"),
    MessageFlag(14, "Model activity"),
    MessageFlag(15, "2.5D Debugging (Project, Translate, Rotate)"),
    MessageFlag(16, "RefFace activity"),
    MessageFlag(17, "Use Count debugging"),
    MessageFlag(18, "Webcut debugging"),
    MessageFlag(19, "Feature Merge / Unmerge debugging"),
    MessageFlag(20, "Parallel meshing activity"),
    MessageFlag(21, "Boundary Layer Tool Debugging"),
    MessageFlag(22, "ExodusMesh sizing function debugging"),
    MessageFlag(23, "Draw after joining chords in WW"),
    MessageFlag(24, "SelfCrossingLoop (and derivatives) debug info"),
    MessageFlag(25, "Extra invalidity checking in WW"),
    MessageFlag(26, "Surface Smoothing"),
    MessageFlag(27, "Primal Construction debugging, see also flag 70"),
    MessageFlag(28, "Plastering debugging"),
    MessageFlag(29, "Volume SubMapping"),
    MessageFlag(30, "Volume Mapping"),
    MessageFlag(31, "CleanUp debugging"),
    MessageFlag(32, "Pyramid debugging"),
    MessageFlag(33, "Whisker Weaving inside chord list face drawing"),
    MessageFlag(34, "If on Whisker Weaving doesn't merge sheets"),
    MessageFlag(35, "If on WW query displays sheets before joining chords"),
    MessageFlag(36, "Enable/Disable idr_keyword_debugger function"),
    MessageFlag(37, "Superdrive debugging"),
    MessageFlag(38, "WW hex formation messages"),
    MessageFlag(39, "Doublet Pillower graphics output"),
    MessageFlag(40, "Element Quality debugging output"),
    MessageFlag(41, "Check_join adjacency check"),
    MessageFlag(42, "Auto vertex type and sweep verification"),
    MessageFlag(43, "Programmer Errors for SubMapping"),
    MessageFlag(44, "Submapping Graphics Debugging"),
    MessageFlag(45, "Pillow Sheet debugging"),
    MessageFlag(46, "Paver breakout detection (expensive)"),
    MessageFlag(47, "Extra LP debugging (see flag 8 also)"),
    MessageFlag(48, "Geometry sizing function."),
    MessageFlag(49, "Draws Face by Face Creation in Paving"),
    MessageFlag(50, "Debugging for AutoSchemeSelect"),
    MessageFlag(51, "Modified Paver (row by row, more intersection checking)"),
    MessageFlag(52, "User Interface: If flag is enabled, filenames being\n"
                "\t\t\tused for input will be echoed and each input line\n"
                "\t\t\twill be echoed prior to being parsed."),
    MessageFlag(53, "Surface Morpher debugging"),
    MessageFlag(54, "Parser debugging"),
    MessageFlag(55, "Stairtool general debugging"),
    MessageFlag(56, "Stairtool face intersection debugging"),
    MessageFlag(57, "Relative Interval/Length setting"),
    MessageFlag(58, "StcVertex debugging of Whisker Weaving" ),
    MessageFlag(59, "VGI developer error"),
    MessageFlag(60, "StcVertex debugging of Looping" ),
    MessageFlag(61, "List number of points used in curve faceting" ),
    MessageFlag(62, "Print verbose information on group operations"),
    MessageFlag(63, "Label Whisker Weaving diagrams tersely"),
    MessageFlag(64, "No label on Whisker Weaving diagrams"),
    MessageFlag(65, "Volume Morpher debugging"),
    MessageFlag(66, "Print debug information on importing Pro/E geometry"),
    MessageFlag(67, "List number of triangles used in surface faceting" ),
    MessageFlag(68, "Debug information on tetrizing volumes" ),
    MessageFlag(69, "Debug information on tet mesher geometry access" ),
    MessageFlag(70, "STC Pillowing, see also flag 27" ),
    MessageFlag(71, "Hex-Tet Debug Flag"),
    MessageFlag(72, "DoubletPillower text messages"),
    MessageFlag(73, "Auto Surface debugging (use new auto surf select)"),
    MessageFlag(74, "Feature-based decomposition info"),
    MessageFlag(75, "Many-to-many sweep imprint debugging"),
    MessageFlag(76, "Virtual point and partition curve"),
    MessageFlag(77, "Volume interval matching"),
    MessageFlag(78, "Tipton Smoother jacobian modification enabler"),
    MessageFlag(79, "VoidCleanupTool debugging"),
    MessageFlag(80, "Hex-Tet informational messages"),
    MessageFlag(81, "Curve Morpher Debugging"),
    MessageFlag(82, "Diamond debugging"),
    MessageFlag(83, "AutoSizer debugging"),
    MessageFlag(84, "Surface auto decomposition"),
    MessageFlag(85, "U-SubMapping debugging"),
    MessageFlag(86, "Virtual curve and partition surface"),
    MessageFlag(87, "Composite curve and composite surface"),
    MessageFlag(88, "Volume partitioning"),
    MessageFlag(89, "Tet meshing warning and debug messages"),
    MessageFlag(90, "Geometry attributes"),
    MessageFlag(91, "Smoothing Debug Output"),
    MessageFlag(92, "Print name changed warnings"),
    MessageFlag(93, "Hex Fix Up"),
    MessageFlag(94, "Entity name attribute"),
    MessageFlag(95, "Group imprint errors"),
    MessageFlag(96, "GraftTool debugging"),
    MessageFlag(97, "Quality details"),
    MessageFlag(98, "Color code imported THEX meshes"),
    MessageFlag(99, "Geometry creation"),
    MessageFlag(100, "Skew Control debugging"),
    MessageFlag(101, "Parsing debugging"),
    MessageFlag(102, "CAEntityId debugging"),
    MessageFlag(103, "Print compact interval assignment constraints"),
    MessageFlag(104, "Report interval matching progress"),
    MessageFlag(105, "Mesh Database (MeshContainer) debugging"),
    MessageFlag(106, "Mesh Cleaver debugging"),
    MessageFlag(107, "Midpoint_subdivision debugging"),
    MessageFlag(108, "Simulog tetmesher debugging"),
    MessageFlag(109, "Transition schemes debugging"),
    MessageFlag(110, "Mesh Defined Geometry"),
    MessageFlag(111, "Tri mesher debugging"),
    MessageFlag(112, "Auto Detail Suppression"),
    MessageFlag(113, "Extra Multi-sweep/sweep debugging"),
    MessageFlag(114, "Blend Finder Debugging"),
    MessageFlag(115, "Exporting Feature Debugging Files"),
    MessageFlag(116, "Sizing function tool data information"),
    MessageFlag(117, "Extra Information on Autoscheme Decision Making"),
    MessageFlag(118, "Blend finding optimization file"),
    MessageFlag(119, "Laminate Tool debugging"),
    MessageFlag(120, "Print unassociated node locations on import mesh"),
    MessageFlag(121, "Print verbose infeasible match interval messages"),
    MessageFlag(122, "Mesh-Based Geometry Debug Information"),
    MessageFlag(123, "Collect memory statistics from Tetmesh-GHS3D"),
    MessageFlag(124, "Print verbose debugging information from Tetmesh-GHS3D"),
    MessageFlag(125, "Mesh refinement debugging"),
    MessageFlag(126, "Surface Splicer debugging"),
    MessageFlag(127, "SculptingTool debug flag"),
    MessageFlag(128, "DualToMeshTool debugging"),
    MessageFlag(129, "Virtual Imprint Debugging"),
    MessageFlag(130, "Hexsheet Insertion Debugging"),
    MessageFlag(131, "Mesh Cutting Debugging"),
    MessageFlag(132, "Global Collection Smoothing"),
    MessageFlag(133, "Print verbose import mesh progress"),
    MessageFlag(134, "TriAdvance or TriDelaunay scheme only"),
    MessageFlag(135, "Keep WhiskerWeave data after meshing"),
    MessageFlag(136, "Stair Step mesh for Sculpting"),
    MessageFlag(137, "GJoin"),
    MessageFlag(138, "Parallel CGM timing"),
    MessageFlag(139, "RTree Debugging"),
    MessageFlag(140, "Facet-based intersection debugging"),
    MessageFlag(141, "Settings save/restore"),
    MessageFlag(142, "Decompose Sweep Debugging"),
    MessageFlag(143, "Decomp Sweep Imprint Debugging"),
    MessageFlag(144, "Medial Axis/Chordal Axis Debugging"),
    MessageFlag(145, "Virtual Geometry Facet Operations"),
    MessageFlag(146, "Sector Tool Meshing Scheme"),
    MessageFlag(147, "Skip Vertex Correction in Paver/CleanUp"),
    MessageFlag(148, "Meshing Benchmarks"),
    MessageFlag(149, "MeshCutting Graphical debugging"),
    MessageFlag(150, "MBG to Acis conversion debugging"),
    MessageFlag(151, "KDDTree Debugging"),
    MessageFlag(152, "Boundary Conditions Debugging"),
    MessageFlag(153, "Print Body information in Geometry operations"),
    MessageFlag(154, "Split Surface Debugging"),
    MessageFlag(155, "Meshing Benchmarks Summary"),
    MessageFlag(156, "CAMAL Paver CleanUp debuging"),
    MessageFlag(157, "Skeleton Sizing Function Debugging (timing & counts)"),
    MessageFlag(158, "Write a CAMAL debug file"),
    MessageFlag(159, "Read a CAMAL debug file"),
    MessageFlag(160, "CAMAL debug file format is binary"),
    MessageFlag(161, "Tread Sweep debugging"),
    MessageFlag(162, "Disable tetmesher field points"),
    MessageFlag(163, "Execute tetmesher until two identical solutions"),
    MessageFlag(164, "Sweep with weighted residual projection"),
    MessageFlag(165, "Sweep with hybrid algorithm"),
    MessageFlag(166, "Unconstrained Paving debugging"),
    MessageFlag(167, "Skeleton Sizing Function Debugging (messages)"),
    MessageFlag(168, "Tweak Target Multiple Debugging "),
    MessageFlag(169, "Enable Knupp affine transformation instead of Roca"),
    MessageFlag(170, "Enable legacy (pre-CAMAL) sweeper"),
    MessageFlag(171, "Enable Quad-Coarsening Edge-Swap inside corners"),
    MessageFlag(172, "Enable Quad-Coarsening Edge-Swap outside corners"),
    MessageFlag(173, "Enable Quad-Coarsening Pinch inside corners"),
    MessageFlag(174, "Enable Quad-Coarsening Pinch outside corners"),
    MessageFlag(175, "Disable creation of crashbackup.cub during crash"),
    MessageFlag(176, "Enable UCP database checking"),
    MessageFlag(177, "Enable Unconstrained Plastering Debug Drawing")

      // IMPORTANT!!!
      // If you add a new debug flag, make sure that you change
      // the result of CubitMessage::number_of_debug_flags().
      // In order to use this type of static initialization,
      // we can't use the sizeof operator, so change it manually.
  };
  debugFlag = staticDebugFlag;

  // Check initialization of debugFlag array.
  for (int i=number_of_debug_flags(); i > 0; i--)
  {
    debugFlag[i].setting = CUBIT_FALSE;
    assert(i == debugFlag[i].flagNumber);
    assert(debugFlag[i].description != NULL);
  }
}

CubitMessage::~CubitMessage()
{
  // Close all streams associated with debug flags.
  // If the same stream is being used for debug and logging, it
  // will get closed below.
  for (int i=number_of_debug_flags(); i > 0; i--)
    remove_debug_stream(i);

  // At this time, the only open streams possible are loggingStream and loggingErrorStream.
  if (loggingStream != NULL)
  {
    loggingStream->close();
    delete loggingStream;
    delete loggingFile;
  }
  if (loggingErrorStream != NULL)
  {
    loggingErrorStream->close();
    delete loggingErrorStream;
    delete loggingErrorFile;
  }

  // Set static instance_ to zero to indicated that we are dead.
  instance_ = 0;
}

int CubitMessage::number_of_debug_flags()
{
  return NUM_DEBUG_FLAGS;
//  return sizeof(debugFlag)/sizeof(debugFlag[0])-1;
}

void CubitMessage::internal_error ( const int message_type,
                                    std::ofstream *output_stream,
                                    const char *format, va_list &args )
{
  int print_it = CUBIT_FALSE;

  const char* prefix = "";

  switch (message_type)
  {
    case CUBIT_ERROR:

       print_it = CUBIT_TRUE;
       prefix = "ERROR: ";
       break;
    case CUBIT_WARNING:
      if (warningFlag)
      {
        print_it = CUBIT_TRUE;
        prefix = "WARNING: ";
      }
      break;
    case CUBIT_INFO:
      if (infoFlag)
        print_it = CUBIT_TRUE;
      break;
    case CUBIT_DIAGNOSTIC:
      if (diagnosticFlag)
      {
        print_it = CUBIT_TRUE;
        prefix = "DIAGNOSTIC: ";
      }
      break;
    default:
      if (message_type >= CUBIT_DEBUG_1 && message_type <= number_of_debug_flags()+10)
      {
        if (debugFlag[message_type-10].setting) print_it = CUBIT_TRUE;
        break;
      }
  }

  if (print_it)
  {
    // Avoid a potential buffer overflow for large strings.
    // Simply don't print the entire command to the commandline
    // and give an indication that there is more.  The entire
    // command is journaled correctly. -- KGM
    char msgbuf[2*BUFSIZ];
    int sz = vsnprintf(msgbuf, 2*BUFSIZ, format, args);
    if (sz > 2*BUFSIZ || sz < 0)
    {
      msgbuf[2*BUFSIZ - 5] = '.';
      msgbuf[2*BUFSIZ - 4] = '.';
      msgbuf[2*BUFSIZ - 3] = '.';
      msgbuf[2*BUFSIZ - 2] = '\n';
      msgbuf[2*BUFSIZ - 1] = '\0'; // just to make sure
    }

      // loggingStream is used to journal error, warning, and info messages.
      // debug messages can also be journalled there by setting the
      // output stream for the debug flag to the same file.
    if (loggingStream != NULL && (message_type == CUBIT_ERROR ||
                                  message_type == CUBIT_WARNING ||
                                  message_type == CUBIT_INFO))
    {
      *loggingStream << prefix << msgbuf;
      loggingStream->flush();
    }
      //loggingErrorStream is used to (if the user has requested it)
      // log only ERROR: messages
    if (loggingErrorStream != NULL && message_type == CUBIT_ERROR)
    {
      *loggingErrorStream << prefix << msgbuf;
      loggingErrorStream->flush();
    }

    if (output_stream == NULL)
    {
      if(CubitMessage::mHandler == NULL)
      {
        std::cout << prefix << msgbuf;
        std::cout.flush();
      }
      else
      {
        if(message_type == CUBIT_ERROR)
        {
          CubitMessage::mHandler->print_error_prefix(prefix);
          CubitMessage::mHandler->print_error(msgbuf);
        }
        else
        {
          CubitMessage::mHandler->print_message_prefix(prefix);
          CubitMessage::mHandler->print_message(msgbuf);
        }
      }
    }
    else
    {
      *output_stream << prefix << msgbuf;
      output_stream->flush();
    }
  }
}

int CubitMessage::print_error ( const char *format, ... )
{
  va_list args;
  va_start(args, format);

#ifdef XTERM
  char esc = 0x1B;
  // Turn on reverse video on VT102 (xterm also)
  // (0=normal, 1-bold, 4-underscore, 5-blink, 7-inverse)
  std::cout << esc << '[' << '7' << 'm';
#endif

  internal_error(CUBIT_ERROR, NULL, format, args);

#ifdef XTERM
  std::cout << esc << '[' << '0' << 'm';
  std::cout.flush();
#endif

  va_end(args);
  add_to_error_count();
  return CUBIT_FAILURE;
}

int CubitMessage::print_warning ( const char *format, ... )
{
  va_list args;
  va_start(args, format);

  internal_error(CUBIT_WARNING, NULL, format, args);

  va_end(args);
  add_to_warning_count();
  return CUBIT_FAILURE;
}

int CubitMessage::print_info ( const char *format, ... )
{
  va_list args;
  va_start(args, format);

  internal_error(CUBIT_INFO, NULL, format, args);

  va_end(args);
  return CUBIT_FAILURE;
}

int CubitMessage::is_debug_flag_set( int flag )
{
   if( DEBUG_FLAG( flag ))
   {
      currentDebugFlag = flag;
      return CUBIT_TRUE;
   }
   return CUBIT_FALSE;
}

int CubitMessage::print_debug( const char *format, ... )
{
  va_list args;
  va_start(args, format);

  internal_error(currentDebugFlag+10,
                 debugFlag[currentDebugFlag].outputStream,
                 format, args);

  va_end(args);
  return CUBIT_FAILURE;
}

void CubitMessage::print_diagnostic ( const char *format, ... )
{
  va_list args;
  va_start(args, format);

  internal_error(CUBIT_DIAGNOSTIC, NULL, format, args);

  va_end(args);
}

int CubitMessage::reset_error_count(int value)
{
  int current_value = errorCount;
  if (errorCount != value) {
    errorCount = value;
    PRINT_WARNING("Error count manually changed from %d to %d\n\n",
		  current_value, value);
  }
  return current_value;
}

int CubitMessage::error_count()
{
  return errorCount;
}

void CubitMessage::add_to_error_count()
{
  errorCount++;
}

int CubitMessage::reset_warning_count(int value)
{
  int current_value = warningCount;
  if (warningCount != value) {
    warningCount = value;
    PRINT_INFO("Warning count manually changed from %d to %d\n\n",
		  current_value, value);
  }
  return current_value;
}

int CubitMessage::warning_count()
{
  return warningCount;
}

void CubitMessage::add_to_warning_count()
{
  warningCount++;
}

void CubitMessage::output_debug_information(int from, int to, int step)
{
  if (to == -1)
    to = number_of_debug_flags();

  PRINT_INFO("Debug Flag Settings "
	     "(flag number, setting, output to, description):\n");
   for (int i=from; i <= to; i+=step) {
      debugFlag[i].output();
   }
  PRINT_INFO("\n");
}

void CubitMessage::output_debug_information(CubitString &match)
{
  int count = 0;
  for (int i=1; i <= number_of_debug_flags(); i++) {
    char *tmp = CubitUtil::util_strdup((char*)(debugFlag[i].description));
    if (tmp && strlen(tmp) > 0) {
      CubitString debug_description(tmp);
      debug_description.to_lower();
      if (debug_description.find(match, 0) < debug_description.length()) {
	if (count == 0) {
	  PRINT_INFO("Debug Flag Settings "
		     "(flag number, setting, output to, description):\n");
	}
	debugFlag[i].output();
	count++;
      }
    }
    CubitUtil::util_strdup_free(tmp);
  }
  if (count == 0) {
    PRINT_WARNING("No debug descriptions contain the "
		  "substring '%s'\n", match.c_str());
  }
  PRINT_INFO("\n");
}

void CubitMessage::output_logging_information()
{
  if (loggingStream != NULL)
     PRINT_INFO("logging           = On, log file = '%s'\n", loggingFile->c_str());
  else
     PRINT_INFO("logging           = Off\n");
  if (loggingErrorStream != NULL)
     PRINT_INFO("logging Errors    = On, log file = '%s'\n",loggingErrorFile->c_str());

}

void MessageFlag::output()
{
  CubitMessage::instance()->
    print_info("%2d  %3s  %-16s   %s\n",
               flagNumber, (setting == 1 ? "ON " : "OFF"),
               (filename == NULL ? "terminal" : filename->c_str()),
               description);
}

int CubitMessage::find_file_use(const CubitString &filename)
{
  if (filename == "terminal") {
    // remove_debug_stream has set the outputStream and filename to NULL.
    return -1;
  }

  // See if any of the other debug flags have this file open
  for (int i=number_of_debug_flags(); i > 0; i--) {
    if (debugFlag[i].filename && *(debugFlag[i].filename) == filename) {
      return i;
    }
  }
  if (loggingFile && *(loggingFile) == filename)
    return -2;

  if (loggingErrorFile && *(loggingErrorFile) == filename)
    return -3;

  return 0;
}

int CubitMessage::count_stream_users(const std::ofstream *stream)
{
  int match = 0;
  if (stream != NULL)
  {
    for (int i=number_of_debug_flags(); i > 0; i--)
    {
      if (debugFlag[i].outputStream == stream)
      {
        match++;
      }
    }

    if (loggingStream == stream)
      match++;
    if (loggingErrorStream == stream)
       match++;
  }
  return match;
}

void CubitMessage::set_logging_file_setting(const CubitString &filename, CubitBoolean resume_flag)
{
  // If logging is currently outputting to a file, close it if
  // it is the only thing using that file. (and the filenames don't match)
  if (loggingFile && *loggingFile == filename)
    return;

  if (loggingErrorFile && *loggingErrorFile == filename)
  {
    PRINT_ERROR("Can't set the logging file to be the same as the Error logging file.\n");
    return;
  }

  int users = count_stream_users(loggingStream);
  if (users == 1) { // Just us...
    loggingStream->close();
    delete loggingStream;
    loggingStream = NULL;
    delete loggingFile;
    loggingFile = NULL;
  }

  int match = find_file_use(filename);

  if (match == -1) // Filename is 'terminal'
    return;
  else if (match != 0)
  {
    loggingFile   = debugFlag[match].filename;
    loggingStream = debugFlag[match].outputStream;
  }
  else
  {
    loggingFile   = new CubitString(filename);
    if(resume_flag)
       loggingStream = new std::ofstream(filename.c_str(), std::ios::out | std::ios::app);
    else
       loggingStream = new std::ofstream(filename.c_str());
  }
}

void CubitMessage::set_debug_file_setting(const int index, const CubitString &filename)
{
  // If this flag is currently outputting to a file, close it if
  // this is the only flag using that file.
  remove_debug_stream(index);

  int match = find_file_use(filename);

  if (match == -1) // Filename is 'terminal'
    return;
  if (match == -2 || match == -3) {// Filename is same as loggingFile or loggingErrorFile;
    debugFlag[index].filename = loggingFile;
    debugFlag[index].outputStream = loggingStream;
  }
  else if (match == index)
    return;
  else if (match != 0) {
    debugFlag[index].filename = debugFlag[match].filename;
    debugFlag[index].outputStream = debugFlag[match].outputStream;
  }
  else {
    debugFlag[index].filename = new CubitString(filename);
    debugFlag[index].outputStream = new std::ofstream(filename.c_str());
  }
}

void CubitMessage::remove_debug_stream(const int index)
{
  // NOTE: DO NOT USE PRINT_* CALLS, THIS IS CALLED FROM DESTRUCTOR.

  // Multiple debug flags may be using the same output stream,
  // Go through the list and count who is using this stream,
  // If only one use, close and delete the stream.
  if (debugFlag[index].outputStream == NULL)
    return;

  int match = count_stream_users(debugFlag[index].outputStream);

  if (match == 1) {
    debugFlag[index].outputStream->close();
    delete debugFlag[index].outputStream;
    delete debugFlag[index].filename;
  }
  debugFlag[index].filename = NULL;
  debugFlag[index].outputStream = NULL;
}

CubitBoolean CubitMessage::Interrupt()
{
  return cubit_intr;
}

void CubitMessage::set_message_handler(CubitMessageHandler *handler)
{
  CubitMessage::mHandler = handler;
}

CubitMessageHandler* CubitMessage::get_message_handler()
{
  return CubitMessage::mHandler;
}

MessageFlag::MessageFlag(int flag_number, const char *desc)
    : flagNumber(flag_number), setting(CUBIT_FALSE),
      description(desc), filename(NULL), outputStream(NULL)
{
}

MessageFlag::MessageFlag()
{
  flagNumber   = 0;
  setting      = CUBIT_FALSE;
  description  = NULL;
  filename     = NULL;
  outputStream = NULL;
}


void CubitMessage::set_logging_file_setting(char* filename)
{
  if (loggingFile && *loggingFile == filename)
     return;

  if (CubitUtil::compare(filename,"terminal")) { // Filename is 'terminal'
    if (loggingStream != NULL) {
      loggingStream->close();
      delete loggingStream;
      loggingStream = NULL;
      delete loggingFile;
      loggingFile = NULL;
    }
    return;
  }
  else {
    loggingFile   = new CubitString(filename);
    loggingStream = new std::ofstream(filename, std::ios::out | std::ios::app );
  }
}

void CubitMessage::set_error_logging_file_setting(char* filename)
{
  if (loggingErrorFile && *loggingErrorFile == filename)
     return;

  if (CubitUtil::compare(filename,"terminal")) { // Filename is 'terminal'
    if (loggingErrorStream != NULL) {
      loggingErrorStream->close();
      delete loggingErrorStream;
      loggingErrorStream = NULL;
      delete loggingErrorFile;
      loggingErrorFile = NULL;
    }
    return;
  }
  else {
    loggingErrorFile   = new CubitString(filename);
    loggingErrorStream = new std::ofstream(filename, std::ios::out | std::ios::app );
  }
}

//Initialize all settings in this class
void CubitMessage::initialize_settings()
{

  SettingHandler::instance()->add_setting("Info",
                                          CubitMessage::set_info_flag,
					  CubitMessage::get_info_flag);

  /*SettingHandler::instance()->add_setting("Logging",
                                          CubitMessage::set_logging_file_setting,
					  CubitMessage::get_logging_file_setting);*/

  SettingHandler::instance()->add_setting("Diagnostic",
					 CubitMessage::set_diagnostic_flag,
					 CubitMessage::get_diagnostic_flag);

  SettingHandler::instance()->add_setting("Warning",
					  CubitMessage::set_warning_flag,
					  CubitMessage::get_warning_flag);
}

CubitString CubitMessage::logging_filename() const
{
  CubitString temp_string;
  if(loggingStream != NULL)
    temp_string = loggingFile->c_str();

  return temp_string;
}

CubitString CubitMessage::logging_errors_filename() const
{
  CubitString temp_string;
  if(loggingErrorStream != NULL)
     temp_string = loggingErrorFile->c_str();

  return temp_string;
}



MessageFlag::~MessageFlag()
{
  // It is not safe to delete either the stream or the filename here
  // since multiple instances (debug flags) may refer to the same memory
}

#ifdef STANDALONE
void main() {
CubitMessage::instance()->output_debug_information(1, 10, 2);
CubitMessage::instance()->output_debug_information(12);
CubitMessage::instance()->set_debug_file(5, "Debug_Test.file");
DEBUG_FLAG(5, CUBIT_TRUE);
PRINT_DEBUG_5("This is a test\n");
CubitMessage::instance()->output_debug_information(5,5);
}
#endif
