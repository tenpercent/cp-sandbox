//- Class:          CubitFileUtil
//- Description:    Class with functions to get files from a directory, etc.
//- Owner:          Steve Storm

#ifndef CUBITFILEUTIL_HPP
#define CUBITFILEUTIL_HPP

#include "DLIList.hpp"
#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitUtilConfigure.h"

class CUBIT_UTIL_EXPORT CubitFileUtil
{
public:

  CubitFileUtil();
  //- Constructor
  
   ~CubitFileUtil();
  //- Destructor

   static CubitStatus get_files( const char *path,
                                 char *dir,
                                 DLIList<CubitString*> &filenames,
                                 CubitBoolean include_dirs = CUBIT_FALSE );
   // Gets the files in the given directory that match the given wildcard (case insensitive)
   // specification (sent in as part of path).  Returns directory string.
   // CubitStrings in filenames list do not include path to files, only
   // the file names.  If directories are requested, directories will be
   // included and will have slash (or backslash on NT) at the end.  List
   // is currently returned unsorted.
   // NOTE:  Be sure to free the allocated CubitStrings in filenames when 
   //        you are done.

   static CubitStatus get_current_working_directory( char *wd );
   // Get current working directory (send in char[PATH_MAX])

   static CubitStatus get_full_path_str( const char *part, 
                                         CubitString &full_path );
   // Given a relative path, return the full path to the filename or 
   // directory.

   static void make_path_platform_compatible( char *path );
   // Replaces slashes with backslashes or vice versa in the given path.

   static CubitString get_nice_filename( const char *file_in );
   // Gets a nice looking filename.  All it does is return just the filename
   // without the path if the input filename is in the current working
   // directory (ie., so the filename doesn't contain the full path).

   static void split_path( const char *path, char *dirpart, char *filepart );
   // Split a path string into directory and filename parts.  The function 
   // returns the directory part with a slash (or backslash for NT) on the end.
   // Note the logic is based entirely on looking for the last directory separator 
   // (slash or backslash) in the path - only the string is used (not lstat, for 
   // instance).  Because of this if you know you have a directory coming in put 
   // a slash on the end of your path, otherwise the last directory in your path 
   // will be parsed out as a file.

   static CubitStatus get_default_cubit_file_name( char *filename );
   // Gets the default Cubit filename, when the user doesn't specify one (ie., 
   // at Cubit startup.  Format cubitXX.cub, where XX = integer.  Returned 
   // filename does not include full path (assumed current working directory).

   static int get_next_backup_filenumber( const char *basename,
                                          int &num_backups,
                                          int &first_backup_num );
   // Gets the next backup filenumer (i.e., file.cub.2, if file.cub.1 already
   // exists).  Also returns the number of current backups in the directory,
   // and the first backup number found.  Handles cases like file.cub.3.4, etc.. 
   // (user started from a backup file).
   // Returns -1 if failure.

private:
  static CubitBoolean all_chars_are( char ch, const char *str );
  // Check if all chars from src are ch

  static CubitBoolean is_int_number( const char *str );
  // Check if the string represents an integer number
  
  static CubitBoolean contains_char( char ch, const char *str );
  // Check if the string contains the given character
};


// The FRegExp (Filename Regular Expression) class was adapted 
// from an article in the April 2000 C/C++ Users Journal, by Michal 
// Niklas.  It is used to match filenames against wildcard expressions.  
// It handles '*' and '?'.
class FRegExp
{
public:
  
  //- Heading: Constructors and Destructor
  FRegExp();
  //- Constructor
  
   ~FRegExp();
  //- Destructor

   static int match(const char *re, const char *text);
   // Match text against re (regular expression) (wildcard).
  
#ifdef BOYD15
   static int matchi( const char *reg,  const char *txt);
   // Same except case insensitive
#endif
  
private:

   static int match_here(const char *re, const char *text);
   // Search for re at beginning of text
   static int match_star( const char *reg,  const char *txt);
   // Omit chars until txt[j] = reg[0]

   static int matchi_here( const char *reg,  const char *txt);
   static int matchi_star( const char *reg,  const char *txt);
};

#endif

