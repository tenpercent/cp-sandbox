#ifndef CUBITINPUTFILE_HPP
#define CUBITINPUTFILE_HPP

#include "CubitString.hpp"
#include <stdio.h>
#include "CubitMessage.hpp"

struct CubitInputFile
{
  enum FileType
  {
    FILE_NORMAL=1,
    FILE_FASTQ=2,
    FILE_TEMPORARY=3
  };

  CubitString              filename;
  FILE*                    filePointer;
  int                      lineNumber;
  CubitInputFile::FileType fileType;
  int                      loopCount;
  int                      breakPoint;
  
  CubitInputFile(FILE *file, FileType type=FILE_NORMAL, int loop=1);
  CubitInputFile(const char *fileName,
                 FileType type=FILE_NORMAL,
                 int loop=1,
                 char *default_path=NULL);
  ~CubitInputFile();
  
};

inline CubitInputFile::CubitInputFile(FILE *file,
                                      CubitInputFile::FileType type,
                                      int loop)
  : breakPoint(0)
{
  if (file)
  {
    if (file == stdin)
      filename    = "<stdin>";
    else
      filename    = "<unknown filename>";
    filePointer = file;
  }
  else
  {
    filename = "<Invalid File>";
    filePointer = NULL;
  }
  lineNumber  = 1;
  fileType = type;
  loopCount = --loop;
}

inline CubitInputFile::CubitInputFile(const char *fileName,
                                      CubitInputFile::FileType type,
                                      int loop,
                                      char *includePath)
  : breakPoint(0)
{
  CubitString file_and_path;
  FILE *file = fopen(fileName, "r");

  if (!file && includePath) {
    file_and_path = includePath;
    file_and_path += "/";
    file_and_path += fileName;
    file = fopen(file_and_path.c_str(), "r");
  }

  if (file) {
    if (includePath)
      filename  = file_and_path;
    else
      filename  = fileName;
    filePointer = file;
  }
  else {
    filename = "<Invalid File>";
    filePointer = NULL;
    PRINT_WARNING("Could not open file: %s\n", fileName );
  }
  lineNumber  = 1;
  fileType = type;
  loopCount = --loop;
}

inline CubitInputFile::~CubitInputFile() {
  if (fileType == FILE_TEMPORARY)
    remove(filename.c_str());	/* Delete file if temporary */
}

#endif 

