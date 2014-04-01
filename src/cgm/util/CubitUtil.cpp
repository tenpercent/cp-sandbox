//-------------------------------------------------------------------------
// Filename      : CubitUtil.cc 
//
// Purpose       : This file contains utility functions that can be used
//                 throughout Cubit.
//
// Special Notes : This is a pure virtual class, to prevent instantiation.
//                 All functions are static, called like this:
//                 CubitUtil::function_name();
//
// Creator       : Darryl Melander
//
// Date          : 06/08/98
//
// Owner         : Darryl Melander
//-------------------------------------------------------------------------

#include "CubitUtil.hpp"
#include "CubitString.hpp"
#include "CubitEntity.hpp"
#include "AppUtil.hpp"
#include <ctype.h>
#include <time.h>

#ifdef NT
#include "Windows.h"
#else
#include <unistd.h>
#endif

FILE *CubitUtil::fp = NULL;

void CubitUtil::convert_string_to_lowercase(char *string)
{
   register char *p = string;
   while (*p)
   {
      if (isupper(*p))
         *p = tolower(*p);
      p++;
   }
}

int CubitUtil::strcmp_case_insensitive (const char *s1, const char *s2)
{
   char c1, c2;
   
   do
   {
      c1 = *s1++;
      if(isupper(c1))
         c1 = tolower(c1);
      
      c2 = *s2++;
      if(isupper(c2))
         c2 = tolower(c2);
      
      if(c1 != c2)
         return c1 - c2;
      
   } while(c1 != '\0');
   
   return 0;
}

int CubitUtil::strncmp_case_insensitive (const char *s1, const char *s2,
                                         int n)
{
   char c1, c2;
   
   do
   {
      c1 = *s1++;
      if(isupper(c1))
         c1 = tolower(c1);
      
      c2 = *s2++;
      if(isupper(c2))
         c2 = tolower(c2);
      
      if(c1 != c2)
         return c1 - c2;
      
      n--;
   } while(c1 && n > 0);
   
   return 0;
}

void CubitUtil::list_entity_ids( const char *pre_string, 
                                 const DLIList<CubitEntity*> &entity_list, 
                                 int width, const char *post_string,
                                 int sort, int unique,
                                 int tab, const char *sep_string,
                                 const char *post_string_none )
{
  DLIList <int> id_list( entity_list.size() );
  for ( int i=0; i<entity_list.size(); i++ ) 
    id_list.append( entity_list.next(i)->id() );

  list_entity_ids( pre_string, id_list, width, post_string, sort,
                   unique, tab, sep_string, post_string_none );
}

void CubitUtil::list_entity_ids( const char *pre_string, 
                                 DLIList<int> &id_list,
                                 int width, 
                                 const char *post_string,
                                 int sort, int unique,
                                 int tab_len, const char *sep_string,
                                 const char *post_string_none )
{
  CubitString ret_str;
  process_entity_ids( 1, ret_str, pre_string, id_list, width, post_string,
                      sort, unique, tab_len, sep_string, post_string_none );
}

CubitString CubitUtil::get_entity_ids_str( const char *pre_string, 
                                           DLIList<int> &id_list, 
                                           int width, const char *post_string,
                                           int sort, int unique, int tab_len, 
                                           const char *sep_string,
                                           const char *post_string_none )
{
  CubitString ret_str;

  process_entity_ids( 0, ret_str, pre_string, id_list, width, post_string,
                      sort, unique, tab_len, sep_string, post_string_none );
  return ret_str;
}

void CubitUtil::process_entity_ids( int method,
                                    CubitString &ret_str,
                                    const char *pre_string, 
                                    DLIList<int> &id_list,
                                    int max_len, 
                                    const char *post_string,
                                    int sort, int unique,
                                    int tab_len, const char *sep_string,
                                    const char* post_string_none ) 
{
  // Method: 0 - to a string
  //         1 - to PRINT_INFO
  char temp[200];

  if ( id_list.size() == 0 ) {
    if( method )
      PRINT_INFO("%s%s", pre_string, post_string_none );
    else
    {
      sprintf( temp, "%s%s", pre_string, post_string_none );
      ret_str = temp;
    }
    if( fp )
      fprintf( fp, "%s%s", pre_string, post_string_none );
    return;
  }

  // sort
  if( sort )
  {
    id_list.sort();
    
    // make unique
    if( unique ) 
    {
      int i;
      DLIList <int> id_list_2( id_list );
      id_list_2.reset();
      id_list.clean_out();
      id_list.append( id_list_2.get_and_step() );
      for ( i=id_list_2.size()-1; i--; ) 
      {
        if ( id_list_2.get() != id_list_2.prev() )
          id_list.append( id_list_2.get() );
        id_list_2.step();
      }
    }
  }

  if( max_len < 0 )
    max_len = CUBIT_INT_MAX/2;
    
  // TODO: wrap prestring, if necessary
  if( method )
    PRINT_INFO( "%s", pre_string );
  else
    ret_str = pre_string;
  if( fp )
    fprintf( fp, "%s", pre_string );

  // Keep track of length printed
  int curr_len = strlen(pre_string);
  
  int num = 0;
  int begin = id_list.get();
  int previous = begin;
  int current;
  int comma = 0; // Is comma needed
  int beg_len, prev_len;
  int sep_len = strlen( sep_string );

  // Setup the tab
  char* tab = new char[tab_len+1];
  for( int i=0; i<tab_len; i++ )
     tab[i] = ' ';
  tab[tab_len] = '\0';

  // Loop until all the ids are printed.  Use ranges if possible.
  while( num < id_list.size()+1 )
  {
    current = id_list.get_and_step();
    num++;

    // Handle last entity
    if( num <= id_list.size() )
    {
      if( num==1 ) // Handle 1st time in loop
        continue;
      
      if( current==previous+1 )
      {
        previous = current;
        continue;
      }
    }

    // If we are here, we are no longer tracking a range and
    // need to print the range or a number.
    if( comma )
    {
      if( method )
        PRINT_INFO( sep_string );
      else
        ret_str += sep_string;
      if( fp )
        fprintf( fp, "%s", sep_string );
      curr_len += sep_len;
    }

    if( begin==previous )
    {
      // a single number
      prev_len = int_len(previous);

      if( curr_len+1+prev_len+sep_len > max_len )
      {
        if( method )
        {
          PRINT_INFO( "\n" );
          PRINT_INFO( "%s%d", tab, previous );
        }
        else
        {
          sprintf( temp, "\n%s%d", tab, previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "\n%s%d", tab, previous );
        curr_len = tab_len + prev_len;
      }
      else
      {
        if( comma ) // Don't print space before first item
        {
          if( method )
            PRINT_INFO( " " );
          else
            ret_str += " ";
          if( fp )
            fprintf( fp, " " );
          curr_len++;
        }

        if( method )
          PRINT_INFO( "%d", previous );
        else
        {
          sprintf( temp, "%d", previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "%d", previous );
        curr_len = curr_len + prev_len;
      }
    }
    else if( previous==begin+1 )
    {
      // a range, but only 2 consecutive numbers
      prev_len = int_len(previous);
      beg_len = int_len(begin);

      // Print 1st
      if( curr_len+1+beg_len+sep_len > max_len )
      {
        if( method )
        {
          PRINT_INFO( "\n" );
          PRINT_INFO( "%s%d%s", tab, begin, sep_string );
        }
        else
        {
          sprintf( temp, "\n%s%d%s", tab, begin, sep_string );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "\n%s%d%s", tab, begin, sep_string );
        curr_len = tab_len + beg_len + sep_len;
      }
      else
      {
        if( comma ) // Don't print space before first item
        {
          if( method )
            PRINT_INFO( " " );
          else
            ret_str += " ";
          if( fp )
            fprintf( fp, " " );
          curr_len++;
        }

        if( method )
          PRINT_INFO( "%d%s", begin, sep_string );
        else
        {
          sprintf( temp, "%d%s", begin, sep_string );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "%d%s", begin, sep_string );
        curr_len = curr_len + beg_len + sep_len;
      }

      // Print 2nd
      if( curr_len+1+prev_len+sep_len > max_len )
      {
        if( method )
        {
          PRINT_INFO( "\n" );
          PRINT_INFO( "%s%d", tab, previous );
        }
        else
        {
          sprintf( temp, "\n%s%d", tab, previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "\n%s%d", tab, previous );
        curr_len = tab_len + prev_len;
      }
      else
      {
        if( method )
          PRINT_INFO( " %d", previous );
        else
        {
          sprintf( temp, " %d", previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, " %d", previous );
        curr_len = curr_len + 1+prev_len;
      }
    }
    else
    {
      // a range of 3 or more consecutive numbers
      prev_len = int_len(previous);
      beg_len = int_len(begin);

      if( curr_len+beg_len+prev_len+5+sep_len > max_len )
      {
        if( method )
        {
          PRINT_INFO( "\n" );
          PRINT_INFO( "%s%d to %d", tab, begin, previous );
        }
        else
        {
          sprintf( temp, "\n%s%d to %d", tab, begin, previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "\n%s%d to %d", tab, begin, previous );
        curr_len = tab_len + beg_len+prev_len+4;
      }
      else
      {
        if( comma ) // Don't print space before first item
        {
          if( method )
            PRINT_INFO( " " );
          else
            ret_str += " ";
          if( fp )
            fprintf( fp, " " );
          curr_len++;
        }

        if( method )
          PRINT_INFO( "%d to %d", begin, previous );
        else
        {
          sprintf( temp, "%d to %d", begin, previous );
          ret_str += temp;
        }
        if( fp )
          fprintf( fp, "%d to %d", begin, previous );
        curr_len = curr_len + beg_len+4+prev_len;
      }
    }

    begin = current;
    previous = current;
    comma = 1;

  }

  //TODO: wrap poststring, if required
  if (post_string) {
    
    if( method )
      PRINT_INFO( post_string );
    else
      ret_str += post_string;
    if( fp )
      fprintf( fp, "%s", post_string );
  }
  
  delete [] tab;
}

#define INTABS(n) ((n) >= 0 ? (n) : (-(n)))
int CubitUtil::int_len( int num )
{
  int len = 0; // length of the string to hold the integer number
  unsigned long n; // absolute value of the integer value
  
  // If the number is negative, add 1 for the negative sign
  if (num < 0) len++;
  
  // Loop until the absolute value of the number reaches 0
  n = INTABS(num);
  do {
    // Increment the length and divide the number by 10
    len++;
    n /= 10;
  } while (n);

  return len;
}

namespace
{
  // Unix: Returns the TEMPDIR, TMP, or TEMP directory, whichever is set to a
  // writeable directory.  If neither is set, return /tmp.
  // Windows: Return path returned by GetTempPath() windows function.
  // If it doesn't return a writeable directory, use current directory.
  CubitString get_temp_directory()
  {
#ifdef WIN32

    //get a place to put the temporary file
    const DWORD buf_size = MAX_PATH;
    char temp_path[buf_size];
    DWORD ret_val = GetTempPath(buf_size, temp_path);

    // If the path is not writeable, use the current directory instead
    DWORD atts = GetFileAttributes(temp_path);
#if _MSC_VER > 1200 // after VC6.0
    if (atts == INVALID_FILE_ATTRIBUTES ||      // File doesn't exist
        (atts & FILE_ATTRIBUTE_DIRECTORY) == 0 || // File isn't a directory
        atts & FILE_ATTRIBUTE_READONLY)         // File is read only
#else
    if ((atts & FILE_ATTRIBUTE_DIRECTORY) == 0 || // File isn't a directory
        atts & FILE_ATTRIBUTE_READONLY)         // File is read only
#endif
    {
      if (DEBUG_FLAG(141))
      {
        PRINT_DEBUG_141("\nUsing cwd because ");
#if _MSC_VER > 1200
        if (atts == INVALID_FILE_ATTRIBUTES)
          PRINT_DEBUG_141("directory doesn't exist: %s\n", temp_path);
        else if ((atts & FILE_ATTRIBUTE_DIRECTORY) == 0)
          PRINT_DEBUG_141("file isn't a directory: %s\n", temp_path);
        else if (atts & FILE_ATTRIBUTE_READONLY)         
          PRINT_DEBUG_141("directory is read only: %s\n", temp_path);
#else
        if ((atts & FILE_ATTRIBUTE_DIRECTORY) == 0)
          PRINT_DEBUG_141("file isn't a directory: %s\n", temp_path);
        else if (atts & FILE_ATTRIBUTE_READONLY)         
          PRINT_DEBUG_141("directory is read only: %s\n", temp_path);
#endif
      }
      temp_path[0] = '.';
      temp_path[1] = '\0';
    }
    else
      PRINT_DEBUG_141("\nUsing GetTempPath: %s\n", temp_path);
    return temp_path;

#else

    const char* tmpdir = "/tmp";
    const char* env_tmpdir = getenv("TMPDIR");
    if (!env_tmpdir)
      env_tmpdir = getenv("TMP");
    if (!env_tmpdir)
      env_tmpdir = getenv("TEMP");
    if(env_tmpdir)
      tmpdir = env_tmpdir;

    return tmpdir;

#endif
  }
}

CubitString CubitUtil::get_temporary_filename()
{

  CubitString ret_str;

#ifdef WIN32
  
  //get a place to put the temporary file
  CubitString temp_path = get_temp_directory();

  // make an empty temporary and return the name for it
  char temp_file_name[MAX_PATH];
  if( GetTempFileName(temp_path.c_str(), "CBT", 0, temp_file_name) != 0 )
    ret_str = temp_file_name; 

#else

  CubitString tmpdir = get_temp_directory();
  const char* filepattern = "CBT.XXXXXX";
  char *temp_file_name = new char[tmpdir.length() + strlen(filepattern) + 1];
  sprintf(temp_file_name, "%s/%s", tmpdir.c_str(), filepattern);

  // make an empty file and return the name for it
  int fd = mkstemp(temp_file_name);
  if( fd != -1 )
  {
    ret_str = temp_file_name; 
    // release the open done by mkstemp,
    // temporary file still exists
    close(fd);
  }
  delete [] temp_file_name;

#endif

  return ret_str;
}

void CubitUtil::print_columns( const char* const* array, int count, 
                               const char* indent, size_t step )
{
  assert( step >= sizeof(char*) );
  
  int term_height, term_width, i, j;
  if( !indent ) indent = "";
  
  if( ! AppUtil::instance()->get_terminal_size( term_height, term_width ) )
  {
    for( i = 0; i < count; i++ )
    {
      const char* ptr = *((char**)((long)array + i*step));
      PRINT_INFO("%s%s\n", indent, ptr );
    }
    return;
  }
  
    // find lenth of longest string
  int maxlen = 0;
  for( i = 0; i < count; i++ )
  {
    const char* ptr = *((char**)((long)array + i*step));
    int len = string_length( ptr );
    if( len > maxlen )
      maxlen = len;
  }
  
  char* const line = new char[CUBIT_MAX(maxlen,term_width)+2];
  
    // calculate number of columns of output
  term_width -= string_length(indent);
  int width = maxlen + 1;
  int columns = term_width > width ? term_width / width : 1;
  
    // calculate number of rows of output
  int rows = count / columns;
  if( count % columns )
    rows++;
  
    // calculate the width of one column
  if (columns > 1)
    width = maxlen + (term_width - columns * maxlen) / (columns - 1);
  else
    width = term_width;
  
    // now write output
  for( i = 0; i < rows; i++ )
  {
    int idx;
    const char* str;
    char* ptr = line + sprintf( line, "%s", indent );
    for( j = 0; j < columns - 1; j++ )
    {
      idx = j * rows + i;
      if (idx < count )
        str = *((char**)((long)array + idx*step));
      else
        str = "";
        
      ptr += sprintf( ptr, "%-*s", width, str);
    }
    
    idx = (columns - 1) * rows + i;
    if (idx < count )
      sprintf(ptr, "%s\n", *((char**)((long)array + idx*step)));
    else
      sprintf(ptr, "\n");
    
    PRINT_INFO( "%s", line );
  }
  
  delete [] line;
}

int CubitUtil::string_length( const char* string, int tabsize )
{
  int result = 0;
  for( ; *string ; string++ )
  {
    if( *string == '\t' )
      result += tabsize;
    else if( *string >= ' ' )
      result ++;
  }
  return result;
}


    //does the same thing as strdup... strdup is not supported by some
    // compilers
char* CubitUtil::util_strdup(const char *s1)
{
#ifdef CUBIT_NO_STRDUP
  int len = strlen(s1)+1;
  char* ret_char = (char*) malloc ( (unsigned long) len * sizeof(char));
  strcpy(ret_char, s1);
  return ret_char;
#else
#ifdef NT
  return _strdup(s1);
#else   
  return strdup(s1);
#endif  
#endif  
}

void CubitUtil::cubit_sleep(int duration_in_seconds)
{
#ifdef NT
  ::Sleep(duration_in_seconds*1000);
#else
  sleep(duration_in_seconds);
#endif
}
