//- Class: CubitMessage
//- Description: CubitMessage class - used for reporting messages to
//-              the user
//- Owner: Tim Tautges
//- Checked By:
//- Version:

#ifndef CUBITMESSAGE_HPP
#define CUBITMESSAGE_HPP

// *****************
// Class Declaration
// *****************

#include <fstream>

#include <stdarg.h>
#include "CubitDefines.h"
#include "CubitUtilConfigure.h"


// GCC compiler checking of format arguments
// Note:  For class member functions, the first argument is the
//        implicit 'this' pointer, so if the format string is the
//        first explicit argument A should be 2, not 1.
#ifdef __GNUC__
# define PRINTF_FORMAT(A,B) __attribute__ ((format (printf, A, B)))
#else
# define PRINTF_FORMAT(A,B)
#endif

// CUBIT debug flags; see CubitMessage.cc for a definition of each flag

#define CUBIT_ERROR                        0
#define CUBIT_WARNING                      1
#define CUBIT_INFO                         2
#define CUBIT_DIAGNOSTIC                   3
#define CUBIT_DEBUG_1                     11

#define PRINT_DEBUG_1 PRINT_DEBUG(1)
#define PRINT_DEBUG_2 PRINT_DEBUG(2)
#define PRINT_DEBUG_3 PRINT_DEBUG(3)
#define PRINT_DEBUG_4 PRINT_DEBUG(4)
#define PRINT_DEBUG_5 PRINT_DEBUG(5)
#define PRINT_DEBUG_6 PRINT_DEBUG(6)
#define PRINT_DEBUG_7 PRINT_DEBUG(7)
#define PRINT_DEBUG_8 PRINT_DEBUG(8)
#define PRINT_DEBUG_9 PRINT_DEBUG(9)
#define PRINT_DEBUG_10 PRINT_DEBUG(10)
#define PRINT_DEBUG_11 PRINT_DEBUG(11)
#define PRINT_DEBUG_12 PRINT_DEBUG(12)
#define PRINT_DEBUG_13 PRINT_DEBUG(13)
#define PRINT_DEBUG_14 PRINT_DEBUG(14)
#define PRINT_DEBUG_15 PRINT_DEBUG(15)
#define PRINT_DEBUG_16 PRINT_DEBUG(16)
#define PRINT_DEBUG_17 PRINT_DEBUG(17)
#define PRINT_DEBUG_18 PRINT_DEBUG(18)
#define PRINT_DEBUG_19 PRINT_DEBUG(19)
#define PRINT_DEBUG_20 PRINT_DEBUG(20)
#define PRINT_DEBUG_21 PRINT_DEBUG(21)
#define PRINT_DEBUG_22 PRINT_DEBUG(22)
#define PRINT_DEBUG_23 PRINT_DEBUG(23)
#define PRINT_DEBUG_24 PRINT_DEBUG(24)
#define PRINT_DEBUG_25 PRINT_DEBUG(25)
#define PRINT_DEBUG_26 PRINT_DEBUG(26)
#define PRINT_DEBUG_27 PRINT_DEBUG(27)
#define PRINT_DEBUG_28 PRINT_DEBUG(28)
#define PRINT_DEBUG_29 PRINT_DEBUG(29)
#define PRINT_DEBUG_30 PRINT_DEBUG(30)
#define PRINT_DEBUG_31 PRINT_DEBUG(31)
#define PRINT_DEBUG_32 PRINT_DEBUG(32)
#define PRINT_DEBUG_33 PRINT_DEBUG(33)
#define PRINT_DEBUG_34 PRINT_DEBUG(34)
#define PRINT_DEBUG_35 PRINT_DEBUG(35)
#define PRINT_DEBUG_36 PRINT_DEBUG(36)
#define PRINT_DEBUG_37 PRINT_DEBUG(37)
#define PRINT_DEBUG_38 PRINT_DEBUG(38)
#define PRINT_DEBUG_39 PRINT_DEBUG(39)
#define PRINT_DEBUG_40 PRINT_DEBUG(40)
#define PRINT_DEBUG_41 PRINT_DEBUG(41)
#define PRINT_DEBUG_42 PRINT_DEBUG(42)
#define PRINT_DEBUG_43 PRINT_DEBUG(43)
#define PRINT_DEBUG_44 PRINT_DEBUG(44)
#define PRINT_DEBUG_45 PRINT_DEBUG(45)
#define PRINT_DEBUG_46 PRINT_DEBUG(46)
#define PRINT_DEBUG_47 PRINT_DEBUG(47)
#define PRINT_DEBUG_48 PRINT_DEBUG(48)
#define PRINT_DEBUG_49 PRINT_DEBUG(49)
#define PRINT_DEBUG_50 PRINT_DEBUG(50)
#define PRINT_DEBUG_51 PRINT_DEBUG(51)
#define PRINT_DEBUG_52 PRINT_DEBUG(52)
#define PRINT_DEBUG_53 PRINT_DEBUG(53)
#define PRINT_DEBUG_54 PRINT_DEBUG(54)
#define PRINT_DEBUG_55 PRINT_DEBUG(55)
#define PRINT_DEBUG_56 PRINT_DEBUG(56)
#define PRINT_DEBUG_57 PRINT_DEBUG(57)
#define PRINT_DEBUG_58 PRINT_DEBUG(58)
#define PRINT_DEBUG_59 PRINT_DEBUG(59)
#define PRINT_DEBUG_60 PRINT_DEBUG(60)
#define PRINT_DEBUG_61 PRINT_DEBUG(61)
#define PRINT_DEBUG_62 PRINT_DEBUG(62)
#define PRINT_DEBUG_63 PRINT_DEBUG(63)
#define PRINT_DEBUG_64 PRINT_DEBUG(64)
#define PRINT_DEBUG_65 PRINT_DEBUG(65)
#define PRINT_DEBUG_66 PRINT_DEBUG(66)
#define PRINT_DEBUG_67 PRINT_DEBUG(67)
#define PRINT_DEBUG_68 PRINT_DEBUG(68)
#define PRINT_DEBUG_69 PRINT_DEBUG(69)
#define PRINT_DEBUG_70 PRINT_DEBUG(70)
#define PRINT_DEBUG_71 PRINT_DEBUG(71)
#define PRINT_DEBUG_72 PRINT_DEBUG(72)
#define PRINT_DEBUG_73 PRINT_DEBUG(73)
#define PRINT_DEBUG_74 PRINT_DEBUG(74)
#define PRINT_DEBUG_75 PRINT_DEBUG(75)
#define PRINT_DEBUG_76 PRINT_DEBUG(76)
#define PRINT_DEBUG_77 PRINT_DEBUG(77)
#define PRINT_DEBUG_78 PRINT_DEBUG(78)
#define PRINT_DEBUG_79 PRINT_DEBUG(79)
#define PRINT_DEBUG_80 PRINT_DEBUG(80)
#define PRINT_DEBUG_81 PRINT_DEBUG(81)
#define PRINT_DEBUG_82 PRINT_DEBUG(82)
#define PRINT_DEBUG_83 PRINT_DEBUG(83)
#define PRINT_DEBUG_84 PRINT_DEBUG(84)
#define PRINT_DEBUG_85 PRINT_DEBUG(85)
#define PRINT_DEBUG_86 PRINT_DEBUG(86)
#define PRINT_DEBUG_87 PRINT_DEBUG(87)
#define PRINT_DEBUG_88 PRINT_DEBUG(88)
#define PRINT_DEBUG_89 PRINT_DEBUG(89)
#define PRINT_DEBUG_90 PRINT_DEBUG(90)
#define PRINT_DEBUG_91 PRINT_DEBUG(91)
#define PRINT_DEBUG_92 PRINT_DEBUG(92)
#define PRINT_DEBUG_93 PRINT_DEBUG(93)
#define PRINT_DEBUG_94 PRINT_DEBUG(94)
#define PRINT_DEBUG_95 PRINT_DEBUG(95)
#define PRINT_DEBUG_96 PRINT_DEBUG(96)
#define PRINT_DEBUG_97 PRINT_DEBUG(97)
#define PRINT_DEBUG_98 PRINT_DEBUG(98)
#define PRINT_DEBUG_99 PRINT_DEBUG(99)
#define PRINT_DEBUG_100 PRINT_DEBUG(100)
#define PRINT_DEBUG_101 PRINT_DEBUG(101)
#define PRINT_DEBUG_102 PRINT_DEBUG(102)
#define PRINT_DEBUG_103 PRINT_DEBUG(103)
#define PRINT_DEBUG_104 PRINT_DEBUG(104)
#define PRINT_DEBUG_105 PRINT_DEBUG(105)
#define PRINT_DEBUG_106 PRINT_DEBUG(106)
#define PRINT_DEBUG_107 PRINT_DEBUG(107)
#define PRINT_DEBUG_108 PRINT_DEBUG(108)
#define PRINT_DEBUG_109 PRINT_DEBUG(109)
#define PRINT_DEBUG_110 PRINT_DEBUG(110)
#define PRINT_DEBUG_111 PRINT_DEBUG(111)
#define PRINT_DEBUG_112 PRINT_DEBUG(112)
#define PRINT_DEBUG_113 PRINT_DEBUG(113)
#define PRINT_DEBUG_114 PRINT_DEBUG(114)
#define PRINT_DEBUG_115 PRINT_DEBUG(115)
#define PRINT_DEBUG_116 PRINT_DEBUG(116)
#define PRINT_DEBUG_117 PRINT_DEBUG(117)
#define PRINT_DEBUG_118 PRINT_DEBUG(118)
#define PRINT_DEBUG_119 PRINT_DEBUG(119)
#define PRINT_DEBUG_120 PRINT_DEBUG(120)
#define PRINT_DEBUG_121 PRINT_DEBUG(121)
#define PRINT_DEBUG_122 PRINT_DEBUG(122)
#define PRINT_DEBUG_123 PRINT_DEBUG(123)
#define PRINT_DEBUG_124 PRINT_DEBUG(124)
#define PRINT_DEBUG_125 PRINT_DEBUG(125)
#define PRINT_DEBUG_126 PRINT_DEBUG(126)
#define PRINT_DEBUG_127 PRINT_DEBUG(127)
#define PRINT_DEBUG_128 PRINT_DEBUG(128)
#define PRINT_DEBUG_129 PRINT_DEBUG(129)
#define PRINT_DEBUG_130 PRINT_DEBUG(130)
#define PRINT_DEBUG_131 PRINT_DEBUG(131)
#define PRINT_DEBUG_132 PRINT_DEBUG(132)
#define PRINT_DEBUG_133 PRINT_DEBUG(133)
#define PRINT_DEBUG_134 PRINT_DEBUG(134)
#define PRINT_DEBUG_135 PRINT_DEBUG(135)
#define PRINT_DEBUG_136 PRINT_DEBUG(136)
#define PRINT_DEBUG_137 PRINT_DEBUG(137)
#define PRINT_DEBUG_138 PRINT_DEBUG(138)
#define PRINT_DEBUG_139 PRINT_DEBUG(139)
#define PRINT_DEBUG_140 PRINT_DEBUG(140)
#define PRINT_DEBUG_141 PRINT_DEBUG(141)
#define PRINT_DEBUG_142 PRINT_DEBUG(142)
#define PRINT_DEBUG_143 PRINT_DEBUG(143)
#define PRINT_DEBUG_144 PRINT_DEBUG(144)
#define PRINT_DEBUG_145 PRINT_DEBUG(145)
#define PRINT_DEBUG_146 PRINT_DEBUG(146)
#define PRINT_DEBUG_147 PRINT_DEBUG(147)
#define PRINT_DEBUG_148 PRINT_DEBUG(148)
#define PRINT_DEBUG_149 PRINT_DEBUG(149)
#define PRINT_DEBUG_150 PRINT_DEBUG(150)
#define PRINT_DEBUG_151 PRINT_DEBUG(151)
#define PRINT_DEBUG_152 PRINT_DEBUG(152)
#define PRINT_DEBUG_153 PRINT_DEBUG(153)
#define PRINT_DEBUG_154 PRINT_DEBUG(154)
#define PRINT_DEBUG_155 PRINT_DEBUG(155)
#define PRINT_DEBUG_156 PRINT_DEBUG(156)
#define PRINT_DEBUG_157 PRINT_DEBUG(157)
#define PRINT_DEBUG_158 PRINT_DEBUG(158)
#define PRINT_DEBUG_159 PRINT_DEBUG(159)
#define PRINT_DEBUG_160 PRINT_DEBUG(160)
#define PRINT_DEBUG_161 PRINT_DEBUG(161)
#define PRINT_DEBUG_162 PRINT_DEBUG(162)
#define PRINT_DEBUG_163 PRINT_DEBUG(163)
#define PRINT_DEBUG_164 PRINT_DEBUG(164)
#define PRINT_DEBUG_165 PRINT_DEBUG(165)
#define PRINT_DEBUG_166 PRINT_DEBUG(166)
#define PRINT_DEBUG_167 PRINT_DEBUG(167)
#define PRINT_DEBUG_168 PRINT_DEBUG(168)
#define PRINT_DEBUG_169 PRINT_DEBUG(169)
#define PRINT_DEBUG_170 PRINT_DEBUG(170)
#define PRINT_DEBUG_171 PRINT_DEBUG(171)
#define PRINT_DEBUG_172 PRINT_DEBUG(172)
#define PRINT_DEBUG_173 PRINT_DEBUG(173)
#define PRINT_DEBUG_174 PRINT_DEBUG(174)
#define PRINT_DEBUG_175 PRINT_DEBUG(175)
#define PRINT_DEBUG_176 PRINT_DEBUG(176)
#define PRINT_DEBUG_177 PRINT_DEBUG(177)
#define NUM_DEBUG_FLAGS 177

#define PRINT_ERROR CubitMessage::instance()->print_error
#define PRINT_WARNING CubitMessage::instance()->print_warning
#define PRINT_INFO CubitMessage::instance()->print_info
#define DIAGNOSTIC CubitMessage::instance()->print_diagnostic
#define DIAGNOSTIC_FLAG CubitMessage::instance()->diagnostic_flag
#define DEBUG_FLAG CubitMessage::instance()->debug_flag
#define GET_INFO_FLAG CubitMessage::instance()->get_info_flag
#define SET_INFO_FLAG CubitMessage::instance()->set_info_flag
#define SET_WARNING_FLAG CubitMessage::instance()->set_warning_flag
#define GET_WARNING_FLAG CubitMessage::instance()->get_warning_flag
#define DEBUG_FLAG_SET CubitMessage::instance()->is_debug_flag_set
#define PRINT_DEBUG(x) if(!DEBUG_FLAG_SET(x));else CubitMessage::instance()->print_debug

class CubitString;
class CubitMessage;

class CUBIT_UTIL_EXPORT MessageFlag
{
  friend class CubitMessage;
public:
  ~MessageFlag();
private:
  MessageFlag();
  MessageFlag(int flag_number, const char *desc);

  void output();

    // Member variables
  int flagNumber;
  int setting;
  const char *description;
  CubitString *filename;
  std::ofstream *outputStream;
};

class CubitMessageHandler
{
public:
  CubitMessageHandler() {}
  virtual ~CubitMessageHandler() {}

  virtual void print_message_prefix(const char *prefix) = 0;
  virtual void print_message(const char *message) = 0;
  virtual void print_error_prefix(const char *prefix) = 0;
  virtual void print_error(const char *message) = 0;
};

class CUBIT_UTIL_EXPORT CubitMessage
{
protected:

  static CubitMessage* instance_;
  //- static pointer to unique instance of this class

  static CubitMessageHandler* mHandler;
  //- static pointer to the message output handler

  static MessageFlag staticDebugFlag[];

  MessageFlag *debugFlag;
  //- debug flag, used with internal_error

  static int infoFlag;
  //- info flag, used with internal_error

  static int warningFlag;
  //- warning flag, used with internal_error

  static int diagnosticFlag;
  //- diagnostic flag, used with internal_error

  int currentDebugFlag;

  static int errorCount;
  //- static variable to track the errors that occured in the
  //- a session.  Only gets set when PRINT_ERROR is called.

  static int warningCount;
  //- static variable to track the warnings that occured in the
  //- a session.  Only gets set when PRINT_WARNING is called.

  static std::ofstream *loggingStream;
  static CubitString *loggingFile;
  //- Stream pointer for logging of internal_error messages.
  //- If NULL, output goes to terminal only (except for debug)
  //- If Non-NULL, output goes to both terminal and stream.

  static std::ofstream *loggingErrorStream;
  static CubitString *loggingErrorFile;
  //- Stream pointer for logging of only ERROR messages.
  //- If NULL, ERROR output goes to normal places only
  //- If Non-NULL, output goes to both this stream, and all other places.

  void add_to_error_count();
  //- Increments the errorCount variable. Keep private (GDS).

  void add_to_warning_count();
  //- Increments the errorCount variable. Keep private (GDS).
  void set_debug_stream(const int index, std::ofstream *output_stream);
  //- Set the output stream for this debug flag to output_stream.

  void remove_debug_stream(const int index);
  //- Close and delete the stream if only one use.

  int find_file_use(const CubitString &filename);
  int count_stream_users(const std::ofstream *stream);

  CubitMessage ();
    //- Class Constructor. (Not callable by user code. Class is constructed
    //- by the {instance()} member function.

public:

  static CubitMessage* instance();
  //- Controlled access and creation of the sole instance of this class.

  virtual ~CubitMessage();
  //- Class Destructor.


  void set_logging_file_setting(const CubitString &filename, CubitBoolean resume_flag = CUBIT_FALSE);
  void set_debug_file_setting(const int index, const CubitString &filename);

  CubitString logging_filename() const;
  CubitString logging_errors_filename() const;

  int is_debug_flag_set( int flag );
  //- for use with the PRINT_DEBUG macro only

  //static int get_debug_for_setting_handler(const int index)
  //{return staticDebugFlag[index].setting;};
  // static void set_debug_for_setting_handler(const int index, const int value)
  //{staticDebugFlag[index].setting = value;};
  //- for use with SettingHandler.cpp code only

  int  debug_flag(const int index);
  void debug_flag(const int index, const int flag);
  int  number_of_debug_flags();
  //- debug flag, used with internal_error

  static bool get_info_flag();
  static void set_info_flag(bool flag);
  //- info flag, used with internal_error

  static bool get_warning_flag();
  static void set_warning_flag(bool flag);
  //- warning flag, used with internal_error

  static bool get_diagnostic_flag();
  static void set_diagnostic_flag(bool flag);
  //- diagnostic flag, used with internal_error

  virtual void internal_error(const int message_type, std::ofstream *output_stream,
                              const char *format, va_list &args);
  //- write out a debug/info/error/warning message

  int print_error(const char *format, ... ) PRINTF_FORMAT(2,3);
  //- write out an error message

  int print_warning(const char *format, ... ) PRINTF_FORMAT(2,3);
  //- write out a warning message

  int print_info(const char *format, ... ) PRINTF_FORMAT(2,3);
  //- write out an info message

  int print_debug( const char *format, ... ) PRINTF_FORMAT(2,3);
  //- write out a debug message

  void print_diagnostic(const char *format, ... ) PRINTF_FORMAT(2,3);
  //- write out a diagnostic message

  int reset_error_count(int value = 0);
  //- Sets the errorCount variable to 0;
  //- Returns current value of variable.

  int error_count();
  //- Returns the value of the errorCount variable;
  //- This errorCount variable is incremented only if print_error is called;
  //- there is not a public  interface to only set the flag.
  //- My reasoning for that is that there should be
  //- some notification of an error so that the
  //- user can figure out why this function returns TRUE.

  int reset_warning_count(int value = 0);
  //- Sets the warningCount variable to 0;
  //- Returns current value of variable.

  int warning_count();
  //- Returns the value of the warningCount variable;
  //- This warningCount variable is incremented only if print_warning is called;
  //- there is not a public  interface to only set the flag.
  //- My reasoning for that is that there should be
  //- some notification of an error so that the
  //- user can figure out why this function returns TRUE.

  void output_debug_information(int from=1, int to=-1, int step=1);
  void output_debug_information(CubitString &match);
  void output_logging_information();

  static char* get_logging_file_setting();
  static void set_logging_file_setting(char* file);
  static void set_error_logging_file_setting(char* file);

  static void initialize_settings();

  virtual CubitBoolean Interrupt();
    //- passes back value of interrupt flag (see CubitDefines.h for how
    //- this flag is stored)

  static void set_message_handler(CubitMessageHandler *handler);
  static CubitMessageHandler* get_message_handler();

}; // End of Class CubitMessage

inline int
CubitMessage::debug_flag(const int index)
{return debugFlag[index].setting;}

inline void
CubitMessage::debug_flag(const int index, const int flag)
{debugFlag[index].setting = flag;}

inline bool
CubitMessage::get_info_flag()
{return !!infoFlag;}

inline void
CubitMessage::set_info_flag(bool flag)
{infoFlag = flag;}

inline bool
CubitMessage::get_warning_flag()
{return !!warningFlag;}

inline void
CubitMessage::set_warning_flag(bool flag)
{warningFlag = flag;}

inline bool
CubitMessage::get_diagnostic_flag()
{return !!diagnosticFlag;}

inline void
CubitMessage::set_diagnostic_flag(bool flag)
{diagnosticFlag = flag;}

#endif

