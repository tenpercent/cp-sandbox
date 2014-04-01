//- Class: GetLongOpt
//- Description: GetLongOpt manages the definition and parsing of 
//- long options. Long options can be abbreviated as long as there 
//- is no ambiguity. If an option requires a value, the value should
//- be separated from the option either by whitespace or an "=".
//- 
//- Other features: $
//- o	GetLongOpt can be used to parse options given through environments.$
//- o	GetLongOpt provides a usage function to print usage.$
//- o	Flags & options with optional or mandatory values are supported.$
//- o	The option marker ('-' in Unix) can be customized.$
//- o	Parsing of command line returns optind (see getopt(3)).$
//- o	Descriptive error messages.$
//-
//- Author: S Manoharan. Advanced Computer Research Institute. Lyon. France
//- Owner: Greg Sjaardema
//- Checked By:
//- Version $Id: 

#ifndef GETLONGOPT_HPP
#define GETLONGOPT_HPP

#include "CubitDefines.h"

#include <iostream>

#include <string.h>
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT GetLongOpt {
public:
  enum OptType { 
     Toggle, Valueless, OptionalValue, MandatoryValue
   };

  static const char* const TOGGLE_ON;
  static const char* const TOGGLE_OFF;

  GetLongOpt(const char optmark = '-'); //- Constructor
  //- Constructor for GetLongOpt takes an optional argument: the option
  //- marker. If unspecified, this defaults to '-', the standard (?)
  //- Unix option marker. 

  ~GetLongOpt();                        //- Destructor
  
  int parse(int argc, char * const *argv);
  int parse(const char * const str, const char *p);
  //- {GetLongOpt::parse} is overloaded. It can either parse a string of
  //- options (typically given from the environment), or it can parse
  //- the command line args ({argc}, {argv}). In either case a return
  //- value < 1 represents a parse error. Appropriate error messages
  //- are printed when errors are seen. {GetLongOpt::parse}, in its first
  //- form, takes two strings: the first one is the string to be
  //- parsed and the second one is a string to be prefixed to the
  //- parse errors. In ts second form, {GetLongOpt::parse} returns the
  //- the {optind} (see @i{getopt(3)}) if parsing is successful.

  int enroll(const char * const opt, const OptType t,
	     const char * const desc, const char * const val);
  //- Add an option to the list of valid command options.$
  //- {GetLongOpt::enroll} adds option specifications to its internal
  //- database. The first argument is the option sting. The second
  //- is an enum saying if the option is a flag ({GetLongOpt::Valueless}),
  //- if it requires a mandatory value ({GetLongOpt::MandatoryValue}) or
  //- if it takes an optional value ({GetLongOpt::OptionalValue}).
  //- The third argument is a string giving a brief description of
  //- the option. This description will be used by {GetLongOpt::usage}.
  //- GetLongOpt, for usage-printing, uses {$val} to represent values
  //- needed by the options. {<$val>} is a mandatory value and {[$val]}
  //- is an optional value. The final argument to {GetLongOpt::enroll}
  //- is the default string to be returned if the option is not
  //- specified. For flags (options with {Valueless}), use "" (empty
  //- string, or in fact any arbitrary string) for specifying {TRUE}
  //- and 0 (null pointer) to specify {FALSE}.
  
  const char *retrieve(const char * const opt) const;
  //- Retrieve value of option$
  //- The values of the options that are enrolled in the database
  //- can be retrieved using {GetLongOpt::retrieve}. This returns a string
  //- and this string should be converted to whatever type you want.
  //- See @i{atoi}, @i{atof}, @i{atol} etc.
  //- If a "parse" is not done before retrieving all you will get
  //- are the default values you gave while enrolling!
  //- Ambiguities while retrieving (may happen when options are
  //- abbreviated) are resolved by taking the matching option that 
  //- was enrolled last. For example, -{v} will expand to {-verify}.$
  //- If you try to retrieve something you didn't enroll, you will
  //- get a warning message. 

  void options(std::ostream &outfile = std::cout) const;
  //- Print command line options only. Called by usage().

  void usage(std::ostream &outfile = std::cout ) const;
  //- Print usage information to {outfile}

  void usage(const char *str)		{ ustring = str; }
  //- Change header of usage output to {str}
  //- GetLongOpt::usage is overloaded. If passed a string "str", it sets the
  //- internal usage string to "str". Otherwise it simply prints the
  //- command usage. 

private:
  struct Cell {
     const char *option;	//- option name
     OptType type;		//- option type
     const char *description;	//- a description of option
     const char *value;	        //- value of option (string)
     int wasSet;                //- 0 if not set by user, 1 if set
     Cell *next;		//- pointer to the next cell
     
     Cell() { option = description = value = 0; wasSet = 0; next = 0; }
   };

  Cell *table;		//- option table
  const char *ustring;	//- usage message
  const char *pname;	//- program basename
  char optmarker;		//- option marker
  
  int enroll_done;		//- finished enrolling
  Cell *last;		//- last entry in option table 
  
  const char *basename(const char *p) const;
  int setcell(Cell *c, const char *valtoken, const char *nexttoken, const char *name);
};

#endif // GETLONGOPT_HPP

