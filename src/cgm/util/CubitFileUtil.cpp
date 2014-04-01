//- Class:          CubitFileUtil
//- Description:    Class with functions to get files from a directory, etc.
//- Owner:          Steve Storm

#include "CubitFileUtil.hpp"
#include "CubitString.hpp"

#ifdef NT
  #include <direct.h>
  #include <io.h>
  #include <sys/types.h>
  #include <sys/stat.h>
  #ifndef PATH_MAX
    #define PATH_MAX _MAX_PATH
  #endif
#else
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <dirent.h>
  #include <stdlib.h>
  #include <sys/param.h>
  #include <unistd.h>
  #include <pwd.h>
#endif


#ifdef NT
static const char* DIR_SEP_STR = "\\";
static const char DIR_SEP_CHAR = '\\';
#else
static const char* DIR_SEP_STR = "/";
static const char DIR_SEP_CHAR = '/';
#endif

//#define DEFAULT_NAME  "cubit"
//#define DEFAULT_EXT   "cub"

CubitStatus
CubitFileUtil::get_current_working_directory( char *wd )
{
#ifdef NT
  if( _getcwd( wd, PATH_MAX ) == NULL )
#else
  if( getcwd( wd, PATH_MAX ) == NULL )
#endif
  {
    PRINT_WARNING( "Unable to get new working directory\n" );
    return CUBIT_FAILURE;
  }
  else
  {
    // Add a slash at the end, if not already there
    int wd_len = strlen(wd);
    if( wd[wd_len-1] != DIR_SEP_CHAR )
    {
      wd[wd_len] = DIR_SEP_CHAR;
      wd[wd_len+1] = '\0';
    }

#ifdef NT
    // Make sure format is compatible with full path format
    CubitString full_path_str;
    if( get_full_path_str( wd, full_path_str ) == CUBIT_SUCCESS )
      strcpy( wd, full_path_str.c_str() );
#endif

    return CUBIT_SUCCESS;
  }
}

CubitStatus
CubitFileUtil::get_full_path_str( const char *part, 
                                  CubitString &full_path_str )
{
  char full[PATH_MAX];
  
  char my_part[PATH_MAX];
  strcpy( my_part, part );
  make_path_platform_compatible( my_part );

#ifdef NT
  if( _fullpath( full, my_part, PATH_MAX ) != NULL )
  
    full_path_str = full;

  else
  {
    PRINT_ERROR( "problem getting full path to %s\n", part );
    return CUBIT_FAILURE;
  }
  
#else

  // UNIX
  if( !contains_char( DIR_SEP_CHAR, my_part ) )
  {
    // Add current working directory
    char wd[PATH_MAX];
    get_current_working_directory( wd );
    
    strcat( wd, my_part );
    
    strcpy( my_part, wd );
    
    full_path_str = my_part;
    
    return CUBIT_SUCCESS;
  }
  else
  {
  
    // Okay, realpath will not work of the file does not exist
    // (our Windows function does).  So strip off the path and
    // just check that (for now assumes user will only be using
    // wildcards on file names).
    
    char dpart[PATH_MAX], fpart[PATH_MAX];
    split_path( my_part, dpart, fpart );

    //PRINT_INFO( "Calling realpath on %s\n", part );
    if( realpath( dpart, full ) != NULL )
    {
      //PRINT_INFO( "Got realpath as %s\n", full );
      
      // Append filename back on
      strcat( dpart, fpart );

      full_path_str = dpart;
    }
    else
    {
      PRINT_ERROR( "problem getting full path to %s\n", my_part );
      return CUBIT_FAILURE;
    }
  }
#endif

  return CUBIT_SUCCESS;
}

void 
CubitFileUtil::make_path_platform_compatible( char *path )
{
#ifdef NT
  // Replace '/' with '\\'
  for( char* path_ptr = path; *path_ptr; path_ptr++ )
    if( *path_ptr == '/' )
      *path_ptr = '\\';
#else
     // Replace '\\' with '/'
  for( char* path_ptr = path; *path_ptr; path_ptr++ )
    if( *path_ptr == '\\' )
      *path_ptr = '/';
#endif
}

CubitString
CubitFileUtil::get_nice_filename( const char *file_in )
{
  CubitString ret_str;

  char dpart_in[PATH_MAX], fpart_in[PATH_MAX];
  split_path( file_in, dpart_in, fpart_in );

  char wd[PATH_MAX];
  get_current_working_directory( wd );

  if( strcmp( dpart_in, wd ) == 0 )
    ret_str = fpart_in;
  else
    ret_str = file_in;

  return ret_str;
}

void
CubitFileUtil::split_path( const char *path, char *dirpart, char *filepart )
{
  const char *fpart;
  
  if( (fpart = strrchr(path, DIR_SEP_CHAR)) == NULL )
  {
    // No separator - could be filename or directory.  We assume
    // it's a directory.
    strcpy( filepart, "." );
    strcpy( dirpart, path );
  }
  else
  {
    strcpy( dirpart, path );
    dirpart[fpart - path + 1] = '\0';
    strcpy( filepart, ++fpart );
  }

  // Add slash on end of dirpart if not already there
  int len = strlen( dirpart );
  if( len && dirpart[len-1] != DIR_SEP_CHAR)
    strcat( dirpart, DIR_SEP_STR );

  return;
}

CubitStatus
CubitFileUtil::get_files( const char *path,
                          char *dirname,
                          DLIList<CubitString*> &filenames,
                          CubitBoolean include_dirs )
{
  CubitString *str;

  char my_path[PATH_MAX];
  strcpy( my_path, path );

  // Replace '/' with '\\' or vice versa
  make_path_platform_compatible( my_path );

  // Get full path
  CubitString full_path_str;
  if( get_full_path_str( my_path, full_path_str ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  //PRINT_INFO( "Full path is '%s'\n", full_path_str.c_str() );

  // Get directory and filename parts of pathname
  char dpart[PATH_MAX];
  char fpart[PATH_MAX];
  split_path( full_path_str.c_str(), dpart, fpart );

  strcpy( dirname, dpart );
  
#ifdef NT
  
  // If path points to a directory, add "\*" to the path, to get all the 
  // files in the directory.
  struct _stat file_info;
  if(!_stat(full_path_str.c_str(), &file_info) && (_S_IFDIR & file_info.st_mode))
  {
    // Add the backslash only if needed
    if (full_path_str.get_at(full_path_str.length()-1) == DIR_SEP_CHAR)
      full_path_str += "*";
    else
    {
      full_path_str += DIR_SEP_STR;
      full_path_str += "*";
    }
  }

  // Note, here full_path_str is actually the wildcard expression we are matching
  // against.

  // Now gather up the filenames from the directory.
  _finddata_t c_file;
  long h_file;

  // Note _findfirst accepts wildcard expressions, so we will only get the
  // files we want.
  if( (h_file = _findfirst(full_path_str.c_str(), &c_file )) != -1L )
  {    
    do 
    {
      if( c_file.attrib & _A_SUBDIR )
      {
        // Directory
        if( include_dirs && !all_chars_are('.', c_file.name) ) // Avoid '.' and '..'
        {
          str = new CubitString( c_file.name ); // Without path
          *str += DIR_SEP_STR; // Add slash at end
          filenames.append( str );
        }
      }
      else
      {
        str = new CubitString( c_file.name ); // Without path
        filenames.append( str );
      }
    } while( _findnext( h_file, &c_file ) == 0 );
    
    _findclose( h_file );
  }
  
#else

  // UNIX

  DIR *dir;

  dir = opendir( dirname );
  if( dir == NULL )
  {
    PRINT_ERROR( "%s\n", strerror(0) );
    return CUBIT_FAILURE;
  }

  dirent *direntry;       // Info about what dir contains 
  struct stat stat_buf;   // Info about particular file (or dir) 
  
  while( (direntry = readdir(dir)) != NULL )
  {
    const char* fname = direntry->d_name;
    
    // Append fullpath to fname
    char full_name[PATH_MAX];
    strcpy( full_name, dirname );
    strcat( full_name, fname );
    
    // Get info about file
    if( stat(full_name, &stat_buf) != 0 )
    {
//      PRINT_WARNING( "%s\n", strerror(0) );
      continue;
    }
    if( (stat_buf.st_mode & S_IFMT) == S_IFDIR )
    {
      // Directory
      if( include_dirs && !all_chars_are('.', fname) ) // Avoid '.' and '..'
      {
        // Do wildcard match.  We do this ourselves for Unix (handles
        // '*' and '?').  Another possiblity would be to use either glob()
        // or fnmatch(), but they may not be available on all platforms.
        if( FRegExp::match( fpart, fname ) )
        {
          str = new CubitString( fname ); // Without path
          *str += DIR_SEP_STR; // Add slash at end
          filenames.append( str );
        }
      }
    }
    else
    {
      // Do wildcard match.  We do this ourselves for Unix (handles
      // '*' and '?').  Another possiblity would be to use either glob()
      // or fnmatch(), but they may not be available on all platforms.
      if( FRegExp::match( fpart, fname ) )
      {
        str = new CubitString( fname ); // Without path
        filenames.append( str );
      }
    }
  }
  closedir( dir );

#endif

  return CUBIT_SUCCESS;
}
 
CubitBoolean 
CubitFileUtil::all_chars_are( char ch, const char *str )
{
  while( *str )
    if( *str++ != ch )
      return CUBIT_FALSE;
    return CUBIT_TRUE;
}

CubitBoolean 
CubitFileUtil::is_int_number( const char *str )
{
  while( *str )
    if( !isdigit(*str++) )
      return CUBIT_FALSE;
    return CUBIT_TRUE;
}

CubitBoolean 
CubitFileUtil::contains_char( char ch, const char *str )
{
  while( *str )
    if( *str++ == ch )
      return CUBIT_TRUE;
    return CUBIT_FALSE;
}

CubitStatus 
CubitFileUtil::get_default_cubit_file_name( char *filename )
{
  // Get a list of all files matching 'cubit*.cub', in current directory
  DLIList<CubitString*> filenames;
  char dir_name[PATH_MAX];
  
  if( get_files( "cubit*.cub", dir_name, filenames, CUBIT_FALSE ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( !filenames.size() )
  {
    strcpy( filename, "cubit01.cub" );
    return CUBIT_SUCCESS;
  }

  // Loop through each file, keeping track of max file number
  int i;
  int curr_num, max_num = 0;
  CubitString *fname_str;
  for( i=filenames.size(); i--; )
  {
    fname_str = filenames.get_and_step();

    char fname[PATH_MAX];
    strcpy( fname, fname_str->c_str() );
    delete fname_str;

    // Extract what's between cubit and ending .cub
    fname[strlen(fname)-4] = '\0';

    char num_str[PATH_MAX];
    strcpy( num_str, fname+5 );

    if( !is_int_number( num_str ) )
      continue;

    curr_num = atoi( num_str );

    if( curr_num > max_num )
      max_num = curr_num;
  }

  max_num++;

  // Place the number in the string
  if( max_num<10 )
    sprintf( filename, "cubit0%d.cub", max_num );
  else
    sprintf( filename, "cubit%d.cub", max_num );

  return CUBIT_SUCCESS;
}

int
CubitFileUtil::get_next_backup_filenumber( const char *p_basename,
                                           int &num_backups,
                                           int &first_backup_id )
{
  num_backups = 0;
  first_backup_id = -1;
  DLIList<CubitString*> filenames;
  
  char basename[PATH_MAX];
  strcpy(basename, p_basename);
  // Replace '/' with '\\' or vice versa
  make_path_platform_compatible( basename );

  // Get a list of all files matching 'basename.*'
  char match_name[PATH_MAX];
  strcpy( match_name, basename );
  strcat( match_name, ".*" );

  char dir_name[PATH_MAX];
  if( get_files( match_name, dir_name, filenames, CUBIT_FALSE ) == CUBIT_FAILURE )
    return -1;

  if( !filenames.size() )
    return 1;

  // Loop through each file, keeping track of which are valid backup files
  // and the maximum number.

  char fname[PATH_MAX];
  split_path( basename, dir_name, fname );
  unsigned int base_len = strlen( fname );

  int i;
  int curr_num, max_num = 0, first_num = CUBIT_INT_MAX;
  CubitString *fname_str;
  for( i=filenames.size(); i--; )
  {
    fname_str = filenames.get_and_step();

    char fname[PATH_MAX];
    strcpy( fname, fname_str->c_str() );
    delete fname_str;

    // Continue if nothing after dot.  Also, I noticed the NT will match 
    // file.cub against wildcard file.cub.*
    if( strlen( fname ) <= base_len+1 )
      continue;

    // Extract what's after the dot
    char aft_dot[PATH_MAX];
    strcpy( aft_dot, fname+base_len+1 );

    if( !is_int_number( aft_dot ) )
      continue;

    // Must be a valid integer - extract it

    curr_num = atoi( aft_dot );

    if( curr_num > max_num )
      max_num = curr_num;

    if( curr_num < first_num )
      first_num = curr_num;

    num_backups++;
  }

  if( first_num != CUBIT_INT_MAX )
    first_backup_id = first_num;

  return max_num+1;
}

#if 0
string expand_path(const string& path)
{
  if (path.length() == 0 || path[0] != '~')
    return path;
  
  const char *pfx = NULL;
  string::size_type pos = path.find_first_of(DIR_SEP_CHAR);
  
  if (path.length() == 1 || pos == 1)
  {
    pfx = getenv("HOME");
    if (!pfx)
    {
      // Punt. We're trying to expand ~/, but HOME isn't set
      struct passwd *pw = getpwuid(getuid());
      if (pw)
        pfx = pw->pw_dir;
    }
  }
  else
  {
    string user(path,1,(pos==string::npos) ? string::npos : pos-1);
    struct passwd *pw = getpwnam(user.c_str());
    if (pw)
      pfx = pw->pw_dir;
  }
  
  // if we failed to find an expansion, return the path unchanged.
  
  if (!pfx)
    return path;
  
  string result(pfx);
  
  if (pos == string::npos)
    return result;
  
  if (result.length() == 0 || result[result.length()-1] != DIR_SEP_CHAR)
    result += DIR_SEP_CHAR;
  
  result += path.substr(pos+1);
  
  return result;
}
#endif


FRegExp::FRegExp()
{}

FRegExp::~FRegExp()
{}

int
FRegExp::match(const char *reg, const char *txt)
{
  return match_here(reg, txt);
} 

int
FRegExp::match_here(const char *reg, const char *txt)
{
  if (reg[0] == '\0')
    return *txt == '\0';
  if (reg[0] == '*')
    return match_star(reg + 1, txt);
  if (reg[0] == DIR_SEP_CHAR && (reg[1] == '*' || reg[1] == '?'))
    return *txt == reg[1];
  if (*txt != '\0' && (reg[0] == '?' || tolower(reg[0]) == tolower(*txt)))
    return match_here(reg + 1, txt + 1);
  return 0;
}

int
FRegExp::match_star(const char *reg, const char *txt)
{
  do
  {
    if (match_here(reg, txt))
      return 1;
  } while (*txt != '\0' && *txt++ != *reg);
  return 0;
}

int
FRegExp::matchi_here(const char *reg, const char *txt)
{
  if (reg[0] == '\0')
    return *txt == '\0';
  if (reg[0] == '*')
    return matchi_star(reg + 1, txt);
  if (reg[0] == DIR_SEP_CHAR && (reg[1] == '*' || reg[1] == '?'))
    return tolower(*txt) == tolower(reg[1]);
  if (*txt != '\0' && (reg[0] == '?' || tolower(reg[0]) == tolower(*txt)))
    return matchi_here(reg + 1, txt + 1);
  return 0;
}

int
FRegExp::matchi_star(const char *reg, const char *txt)
{
  do
  {
    if (matchi_here(reg, txt))
      return 1;
  } while (*txt != '\0' && tolower(*txt++) != tolower(*reg));
  return 0;
}
