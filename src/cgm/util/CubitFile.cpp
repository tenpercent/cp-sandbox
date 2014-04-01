/*******************************************************************************
    COPYRIGHT 2002 CATERPILLAR INC.  ALL RIGHTS RESERVED

    This program is the property of Caterpillar Inc., includes Caterpillar's
    confidential and trade secret information, and is maintained as an
    unpublished copyrighted work.  It is not to be copied or used by others
    except under license from Caterpillar.  This program is also protected as an
    unpublished work in accordance with the copyright act of 1976.  In the event
    of either inadvertent or deliberate publication, Caterpillar Inc. intends to
    maintain copyright protection for this work under the relevant copyright
    laws pertaining to published works.  The inclusion of a copyright notice
    hereon is precautionary only, and does not imply publication or disclosure.

 
 Filename      : CubitFile.cpp

 Purpose       : Implements reading and writing data in the Cubit file format.
           
 Special Notes :

 Creator       : Will A. Helden

 Creation Date : 02/15/02

 Owner         : Will A. Helden

*******************************************************************************/

#include "CubitFile.hpp"
#include "CubitFileIOWrapper.hpp"
#include "CubitFileMetaData.hpp"
#include "CubitFileFEModel.hpp"
#include <string.h>

using namespace NCubitFile;

///////////////////////////////////////////////////////////////////////////////
// Construction/Destruction
///////////////////////////////////////////////////////////////////////////////

/* Note: some platforms define both BIG_ENDIAN and LITTLE_ENDIAN
   so don't do if defined(LITTLE_ENDIAN) here.  Compare to BYTE_ORDER instead. */
#if defined(NT) || defined(DA) || defined(CUBIT_LINUX) || (defined(LITTLE_ENDIAN) && (BYTE_ORDER==LITTLE_ENDIAN)) /* should be little endian platforms */
    const UnsignedInt32 CCubitFile::mintNativeEndian = 0;
#else // big endian platforms
    const UnsignedInt32 CCubitFile::mintNativeEndian = 0xFFFFFFFF;
#endif
    const UnsignedInt32 CCubitFile::mintSizeOfContents =
        sizeof(SCubitFileContentsHeader) / sizeof(UnsignedInt32);
    const UnsignedInt32 CCubitFile::mintSizeOfModel =
        sizeof(SCubitFileModelEntry) / sizeof(UnsignedInt32);

CCubitFile::CCubitFile() 
{
    mpReadFile = mpWriteFile = NULL;
    meErrorState = eSuccess;
    mintNextModelID = 1;

    mpaReadModels = mpaWriteModels = NULL;
    mpaReadModelStat = NULL;
    mpMetaData = NULL;
    mstrReadFileName = mstrWriteFileName = mstrBackupFileName = NULL;
    mpReadFEModel = mpWriteFEModel = NULL;

    mModelBuff.mintNumModels = 0;
    mModelBuff.mpaModelData = NULL;
}

CCubitFile::~CCubitFile()
{
    // auto-save?
    // delete incomplete temp files
    FreeAll();

    // Assume if any FE models are still around now (i.e. EndRead/Write hasn't
    // been called) there was an error, so free their memory.
    if(mpReadFEModel)
        delete mpReadFEModel;
    if(mpWriteFEModel)
        delete mpWriteFEModel;

    if(mModelBuff.mpaModelData)
        delete [] mModelBuff.mpaModelData;
}

void CCubitFile::FreeAll()
// Frees all allocated memory, closes open files, and resets class members to
// their constructed state.
{
    // Free all allocated memory.
    if(mstrReadFileName)
        delete [] mstrReadFileName;
    if(mstrWriteFileName)
        delete [] mstrWriteFileName;
    if(mstrBackupFileName)
        delete [] mstrBackupFileName;
    mstrReadFileName = mstrWriteFileName = mstrBackupFileName = NULL;
    if(mpaReadModels)
        delete [] mpaReadModels;
    if(mpaWriteModels)
        delete [] mpaWriteModels;
    mpaReadModels = mpaWriteModels = NULL;
    if(mpaReadModelStat)
        delete [] mpaReadModelStat;
    mpaReadModelStat = NULL;
    if(mpMetaData)
        delete mpMetaData;
    mpMetaData = NULL;

    // Close any open files.
    try {
        if(mpReadFile)
            fclose(mpReadFile);
        if(mpWriteFile)
            fclose(mpWriteFile);
    }
    catch(...)  { }
    mpReadFile = mpWriteFile = NULL;
}

CCubitFile::EErrorCode CCubitFile::Open(const char* xstrReadFileName,
                                        const char* xstrWriteFileName,
                                        const char* xstrBackupFileName)
// Open a file for reading or writing.
// === NOTES ===
// *  The read file is a file to be read only.  The file must exist and be
//    accessable or this function will fail.  
// *  The write file is an output file.  The file must be creatable or this
//    function will fail.  If the file already exists, it will be overwritten.
//    If the write file name is the same as the read file name, a temporary
//    write file name will be generated to be written to until the completion
//    of the write operation at which time the read file will be deleted and
//    the write file renamed to take its place.
// *  The backup file name is used for preventing the loss of old read files.
//    If the read file was to be deleted automatically by a write operation,
//    specification of a backup name, will cause the read file to be renamed
//    instead.
// *  Any of these three filenames may be NULL, but at least one name, read or
//    write, must be specified.
{
    try {
        if(!xstrReadFileName && !xstrWriteFileName)
            throw ePassedNullPointer;

        mpMetaData = new CMetaData;

        // If a file has been designated for reading, copy its name and try to
        // open it.
        if(xstrReadFileName ) {
            UnsignedInt32 lintFileName = strlen(xstrReadFileName);
            mstrReadFileName = new char[lintFileName + 1];
            strcpy(mstrReadFileName, xstrReadFileName);     

            if( !xstrWriteFileName )
            {
              mpReadFile = fopen(mstrReadFileName, "rb");
              if(!mpReadFile)
                  throw eFileReadError;

              // Test to see if the file is a cubit file... it should begin with
              // the string 'CUBE'.
              unsigned char lachrMagic[4];
              if(fread(lachrMagic, 1, 4, mpReadFile) != 4)
                  throw eFileUnrecognizedFormat;
              if((lachrMagic[0] != 'C') || (lachrMagic[1] != 'U') ||
                  (lachrMagic[2] != 'B') || (lachrMagic[3] != 'E'))
                  throw eFileUnrecognizedFormat;
              
              // Load the file's "table of contents".
              CIOWrapper lIO(mpReadFile, 4, 0);
              lIO.BeginReadBlock(4);
              lIO.Read((UnsignedInt32*)&mReadContents, mintSizeOfContents);
              lIO.EndReadBlock();
              
              if(mReadContents.mintHeaderSchema > 0)
                  throw eFileUnrecognizedFormat;
              
              if(mReadContents.mintNumModels) {
                  mpaReadModels =
                      new SCubitFileModelEntry[mReadContents.mintNumModels];
                  mpaReadModelStat = new EModelStat[mReadContents.mintNumModels];
                  if(!mpaReadModels || !mpaReadModelStat)  throw eMemoryError;

                  lIO.BeginReadBlock(mReadContents.mintModelTableOffset);
                  lIO.Read((UnsignedInt32*)mpaReadModels,
                      mReadContents.mintNumModels * mintSizeOfModel);
                  lIO.EndReadBlock();

                  UnsignedInt32 lintMod;
                  for(lintMod = 0; lintMod < mReadContents.mintNumModels; lintMod++) {
                      mpaReadModelStat[lintMod] = eStatNotWritten;
                      if(mpaReadModels[lintMod].mintModelHandle >= mintNextModelID)
                          mintNextModelID = mpaReadModels[lintMod].mintModelHandle + 1;
                  }
              }
              
              // Read the file-scoped (model) meta-data.
              mpMetaData->ReadMetaData(mpReadFile,
                  mReadContents.mintModelMetaDataOffset, 0,
                  mReadContents.mintHeaderSourceEndian);
            }
            else
            {
               mReadContents.mintHeaderSourceEndian = mintNativeEndian;
               mReadContents.mintHeaderSchema = 0;
               mReadContents.mintNumModels = 0;
               mReadContents.mintModelTableOffset = 0;
               mReadContents.mintModelMetaDataOffset = 0;
               mReadContents.mintActiveFEModel = 0;

              // Copy the backup file name, if specified.
              if(xstrBackupFileName) {
                  UnsignedInt32 lintFileName = strlen(xstrBackupFileName);
                  mstrBackupFileName = new char[lintFileName + 1];
                  strcpy(mstrBackupFileName, xstrBackupFileName);
              }
            }
          }
          else {
              // If there wasn't a file to read, initialize the read contents to
              // nothing.
              mReadContents.mintHeaderSourceEndian = mintNativeEndian;
              mReadContents.mintHeaderSchema = 0;
              mReadContents.mintNumModels = 0;
              mReadContents.mintModelTableOffset = 0;
              mReadContents.mintModelMetaDataOffset = 0;
              mReadContents.mintActiveFEModel = 0;
        }
        
        // If a file has been designated for writing
        if(xstrWriteFileName) {
            UnsignedInt32 lintFileName = strlen(xstrWriteFileName);
            mstrWriteFileName = new char[lintFileName + 10];

            // If the read and write file names don't match, then the write
            // file can be opened by name.
            if(!mstrReadFileName || strcmp(mstrReadFileName, xstrWriteFileName)) {
                mintWriteTempFile = 0;
                strcpy(mstrWriteFileName, xstrWriteFileName);       
                mpWriteFile = fopen(mstrWriteFileName, "wb");
            }
            // Otherwise, generate a temporary file name to write to so that
            // the read file's contents will not be destroyed until it is 
            // verifiable that any new writes are successful.
            else {
                mintWriteTempFile = 1;
                UnsignedInt32 lintTempFile = 0;
                while(!mpWriteFile && (lintTempFile < 0xFF)) {
                    sprintf(mstrWriteFileName, "%s~%.2x.tmp",
                        mstrReadFileName, lintTempFile);
                    mpWriteFile = fopen(mstrWriteFileName, "wb");
                    lintTempFile++;
                }
            }
            if(!mpWriteFile)
                throw eFileWriteError;
            
            // Initialize the write file contents - copy the model inventory
            // (but not actual contents yet) from the read file.
            mWriteContents.mintHeaderSourceEndian = mintNativeEndian;
            mWriteContents.mintHeaderSchema = 0;
            mWriteContents.mintNumModels = mReadContents.mintNumModels;
            mWriteContents.mintModelTableOffset = 0;
            mWriteContents.mintModelMetaDataOffset = 0;
            mWriteContents.mintActiveFEModel = mReadContents.mintActiveFEModel;
            mintWriteBuffNumModels = 0;

            if(mWriteContents.mintNumModels) {
                mpaWriteModels =
                    new SCubitFileModelEntry[mWriteContents.mintNumModels];
                if(!mpaWriteModels)  throw eMemoryError;

                UnsignedInt32 lintMod;
                for(lintMod = 0; lintMod < mWriteContents.mintNumModels; lintMod++) {
                    mpaWriteModels[lintMod].mintModelHandle =
                        mpaReadModels[lintMod].mintModelHandle;
                    mpaWriteModels[lintMod].mintModelOffset = 0;
                    mpaWriteModels[lintMod].mintModelLength = 0;
                    mpaWriteModels[lintMod].mintModelType =
                        mpaReadModels[lintMod].mintModelType;
                    mpaWriteModels[lintMod].mintModelOwner =
                        mpaReadModels[lintMod].mintModelOwner;
                }
            }

            // Initialize the write file by writing its identity and an initial
            // table of contents.
            CIOWrapper lIO(mpWriteFile);
            lIO.BeginWriteBlock(0);
            lIO.Write("CUBE", 4);
            lIO.Write((UnsignedInt32*)&mWriteContents, mintSizeOfContents);
            lIO.EndWriteBlock();

            WriteModelTable();
        }
        else {
            // If there wasn't a file to written, initialize the write contents
            // to nothing.
            mWriteContents.mintHeaderSourceEndian = mintNativeEndian;
            mWriteContents.mintHeaderSchema = 0;
            mWriteContents.mintNumModels = 0;
            mWriteContents.mintModelTableOffset = 0;
            mWriteContents.mintModelMetaDataOffset = 0;
            mWriteContents.mintActiveFEModel = 0;
        }
        
        return eSuccess;
    }

    // Handle all open errors by reseting the file object.
    catch(EErrorCode xeErrorCode) {
        FreeAll();
        return xeErrorCode;
    }
    catch(...)  {
        FreeAll();
        return eUnknownError;
    }
}

CCubitFile::EErrorCode CCubitFile::Close()
// Close any open files, completing any pending write operations, and free any
// allocated resources.  See Open() comments for more info.
{
    try {

        if(mpWriteFile) {
            // Write file scoped (model) metadata and update the file's table
            // of contents one last time.
            UnsignedInt32 lintMetaDataLength;
            mpMetaData->WriteMetaData(mpWriteFile,
                mWriteContents.mintModelMetaDataOffset, lintMetaDataLength);
            WriteModelTable();

// throw away erroneously written files?!?

            // If the written file was a temporary file, delete the original
            // file and replace it with the temporary file.
            if(mintWriteTempFile) {             
                // Close any open files.
                if(mpReadFile)
                    fclose(mpReadFile);
                fclose(mpWriteFile);
                mpReadFile = mpWriteFile = NULL;

                // If there was a backup name specified for the old read file
                // rename it instead of deleting it.
                if(mstrBackupFileName)
                    rename(mstrReadFileName, mstrBackupFileName);
                else
                    remove(mstrReadFileName);
                if(rename(mstrWriteFileName, mstrReadFileName))
                    throw eUnknownError;
            }
        }

        // Close any open files and clean up allocated memory.
        FreeAll();
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode) {
        FreeAll();
        return xeErrorCode;
    }
    catch(...)  {
        FreeAll();
        return eUnknownError;
    }
}

/*UnsignedInt32 CCubitFile::GetModelList(UnsignedInt32& xintNumModels,
                                       const SModelData*& xpaModels)
// Returns a list of all models contained in the file.
{
    SCubitFileModelEntry* lpaModels = NULL;
    if(mpaWriteModels) {
        lpaModels = mpaWriteModels;
        xintNumModels = mWriteContents.mintNumModels;
    }
    else if(mpaReadModels) {
        lpaModels = mpaReadModels;
        xintNumModels = mReadContents.mintNumModels;
    }

    if(!lpaModels || !xintNumModels) {
        xpaModels = NULL;
        return eNotFound;
    }

    if(mModelBuff.mintNumModels < xintNumModels) {
        mModelBuff.mintNumModels = xintNumModels;
        if(mModelBuff.mpaModelData)
            delete [] mModelBuff.mpaModelData;
        mModelBuff.mpaModelData = new SModelData[xintNumModels];
    }
    xpaModels = mModelBuff.mpaModelData;
    
    for(UnsignedInt32 lintModel = 0; lintModel < xintNumModels; lintModel++) {
        mModelBuff.mpaModelData[lintModel].mintModelHandle =
            lpaModels[lintModel].mintModelHandle;
        mModelBuff.mpaModelData[lintModel].mintModelType =
            lpaModels[lintModel].mintModelType;
        mModelBuff.mpaModelData[lintModel].mintModelOwner =
            lpaModels[lintModel].mintModelOwner;
    }
    return eSuccess;
}*/

UnsignedInt32 CCubitFile::IsModel(HModel xintModel)
{
    UnsignedInt32 lintIndex;
    if(mWriteContents.mintNumModels)
        return FindModel(xintModel, mpaWriteModels, mWriteContents.mintNumModels, lintIndex);
    else if(mReadContents.mintNumModels)
        return FindModel(xintModel, mpaReadModels, mReadContents.mintNumModels, lintIndex);
    return 0;  // failure
}

UnsignedInt32 CCubitFile::GetReadModelLength(HModel xintModel,
                                             UnsignedInt32& xintLength)
// Looks up the length of the data block in the cubit filr for the requested
// model.  Returns success if the model is found, otherwise an error is returned.
{
    UnsignedInt32 lintIndex;
    if(FindModel(xintModel, mpaReadModels, mReadContents.mintNumModels, lintIndex)) {
        xintLength = mpaReadModels[lintIndex].mintModelLength;
        return eSuccess;
    }
    else {
        xintLength = 0;
        return eNotFound;
    }
}

UnsignedInt32 CCubitFile::CreateModel(EModelType xeModelType, HModel& xintModel)
// Define a new model in the temp file.
{
    if(!mpWriteFile)  return eFileWriteError;

    // Allocate a new bigger model table, copy the existing table (if any)
    // into the new one and then replace the old one.
    SCubitFileModelEntry* lpaModels =
        new SCubitFileModelEntry[mWriteContents.mintNumModels + 1];
    if(!lpaModels)  return eMemoryError;
    if(mpaWriteModels) {
        memcpy(lpaModels, mpaWriteModels,
            sizeof(SCubitFileModelEntry) * mWriteContents.mintNumModels);
        delete [] mpaWriteModels;
    }
    mpaWriteModels = lpaModels;
    
    // Initialize the new model's table entry.
    mpaWriteModels[mWriteContents.mintNumModels].mintModelHandle =
        xintModel;// = mintNextModelID++;
    mpaWriteModels[mWriteContents.mintNumModels].mintModelOffset = 0;
    mpaWriteModels[mWriteContents.mintNumModels].mintModelLength = 0;
    mpaWriteModels[mWriteContents.mintNumModels].mintModelType = xeModelType;
    mpaWriteModels[mWriteContents.mintNumModels].mintModelOwner = 0;
    mpaWriteModels[mWriteContents.mintNumModels].mintModel64bitOSPad = 0;
    mWriteContents.mintNumModels++;
    
    return eSuccess;
}

UnsignedInt32 CCubitFile::DeleteModel(HModel xintModel)
// Flag an model that exists in the read file not to be copied into the write file.
{
    if(!mpWriteFile)  return eFileWriteError;

    // Flag the model in the old file as marked for deletion.
    UnsignedInt32 lintIndex;
    if(!FindModel(xintModel, mpaReadModels, mReadContents.mintNumModels, lintIndex))
        return eNotFound;
    if(mpaReadModelStat[lintIndex] == eStatWritten)
        return eOrderError;
    mpaReadModelStat[lintIndex] = eStatDelete;

    // Remove the model from the write file's table and close up the position in
    // the model table (if applicable).
    if(FindModel(xintModel, mpaWriteModels, mWriteContents.mintNumModels, lintIndex)) {
        mWriteContents.mintNumModels--;
        if(mWriteContents.mintNumModels &&
            (lintIndex < mWriteContents.mintNumModels)) {
            for(UnsignedInt32 lintCopyTo = lintIndex; lintCopyTo <
                mWriteContents.mintNumModels; lintCopyTo++) {
                memcpy(&mpaWriteModels[lintCopyTo],
                    &mpaWriteModels[lintCopyTo + 1],
                    sizeof(SCubitFileModelEntry));
            }
        }
    }

    // Find any models that have the deleted model as their owner and reset
    // the owner to no owner (0).
    if(mWriteContents.mintNumModels) {
        for(UnsignedInt32 lintModel = 0; lintModel <
            mWriteContents.mintNumModels; lintModel++) {
            if(mpaWriteModels[lintModel].mintModelOwner == xintModel)
                mpaWriteModels[lintModel].mintModelOwner = 0;
        }
    }

    // Delete any meta-data that belonged with the deleted model.
    mpMetaData->ClearMetaDataForID(xintModel);

    return eSuccess;
}

UnsignedInt32 CCubitFile::GetModelOwner(HModel xintModel, HModel& xintOwner)
// Return the owning model of the passed model.
{
    // Flag the model in the old file as marked for deletion.
    UnsignedInt32 lintIndex;
    if(!FindModel(xintModel, mpaWriteModels, mWriteContents.mintNumModels, lintIndex))
        return eNotFound;
    xintOwner = mpaWriteModels[lintIndex].mintModelOwner;
    return eSuccess;
}

UnsignedInt32 CCubitFile::SetModelOwner(HModel xintModel, HModel xintOwner)
// Set the owning model for the passed model.
{
    // Flag the model in the old file as marked for deletion.
    UnsignedInt32 lintIndex;
    if(xintOwner &&  // don't allow the model to have a non-existant owner.
        !FindModel(xintOwner, mpaWriteModels, mWriteContents.mintNumModels, lintIndex))
        return eNotFound;
    if(!FindModel(xintModel, mpaWriteModels, mWriteContents.mintNumModels, lintIndex))
        return eNotFound;
    mpaWriteModels[lintIndex].mintModelOwner = xintOwner;
    return eSuccess;
}

UnsignedInt32 CCubitFile::FindModel(HModel xintModel,
                                    SCubitFileModelEntry* xpaModels,
                                    UnsignedInt32 xintNumModels,
                                    UnsignedInt32& xintIndex)
// Determine the index of the passed model handle in the passed model table.
// Private function, return true/false, no exceptions thrown.
{
    if(xintNumModels) {
        for(UnsignedInt32 lintMod = 0; lintMod < xintNumModels; lintMod++) {
            if(xpaModels[lintMod].mintModelHandle == xintModel) {
                xintIndex = lintMod;
                return 1; // success
            }
        }
    }
    return 0; // failure
}

void CCubitFile::WriteModelTable()
// Write (or rewrite) the model table to the write file.
{
    if(!mpWriteFile)  throw eFileWriteError;

    if(!mWriteContents.mintNumModels)  return;  // no models... nothing to do!

    CIOWrapper lIO(mpWriteFile);

    // Write the model table, if the number of models has increased, write the
    // table to a new location in the file, otherwise reuse the previous
    // table's file space.
    if(mWriteContents.mintNumModels > mintWriteBuffNumModels) {
        mintWriteBuffNumModels = mWriteContents.mintNumModels;
        mWriteContents.mintModelTableOffset = lIO.BeginWriteBlock();
    }
    else
        lIO.BeginRewriteBlock(mWriteContents.mintModelTableOffset, 0);
    lIO.Write((UnsignedInt32*)mpaWriteModels,
        mWriteContents.mintNumModels * mintSizeOfModel);
    lIO.EndWriteBlock();

    // Rewrite the contents header to reflect any model table changes.
    lIO.BeginRewriteBlock(4, 0);
    lIO.Write((UnsignedInt32*)&mWriteContents, mintSizeOfContents);
    lIO.EndWriteBlock();
}

void CCubitFile::CopyModel(UnsignedInt32 xintReadOffset,
                           UnsignedInt32& xintWriteOffset,
                           UnsignedInt32 xintLength,
                           FILE* xpReadFile, FILE* xpWriteFile)
// Copy a model from one file to another.  The copied model is byte identical
// to the original.
{
    // Create a 16 kilobyte memory buffer to try to improve copy IO efficiency
    // over a byte-per-byte copy.
    char* lpachrCopyBuffer = new char[0x4000]; //16KB
    if(!lpachrCopyBuffer)
        throw eMemoryError;
    
    // Position the files pointers for the copy.
    CIOWrapper lReadIO(xpReadFile);
    CIOWrapper lWriteIO(xpWriteFile);
    lReadIO.BeginReadBlock(xintReadOffset);
    xintWriteOffset = lWriteIO.BeginWriteBlock();
    
    // Copy the models 16K at a time.
    UnsignedInt32 lintRemaining = xintLength;
    UnsignedInt32 lintCopyBytes = 0x4000;
    while(lintRemaining) {
        if(lintRemaining > 0x4000)  // More than 16K left, copy 16K.
            lintRemaining -= 0x4000;
        else {  // 16K or less left, copy all of what's left.
            lintCopyBytes = lintRemaining;
            lintRemaining = 0;
        }
        lReadIO.Read(lpachrCopyBuffer, lintCopyBytes);
        lWriteIO.Write(lpachrCopyBuffer, lintCopyBytes);
    }

    // Make sure the copy was complete and free the buffer.
    lReadIO.EndReadBlock();
    if(lWriteIO.EndWriteBlock() != xintLength)
        throw eCorruptBlock;
    delete [] lpachrCopyBuffer;
}


///////////////////////////////////////////////////////////////////////////////
// Geometry Model Read/Write
///////////////////////////////////////////////////////////////////////////////

UnsignedInt32 CCubitFile::BeginWriteGeomModel(HModel xintGeomModel,
                                              const char* xstrGeomFile)
// Copy an external geometry file into the cubit file as a geometry model.
{
    if(!mpWriteFile)  return eFileWriteError;

    // Try to open the geometry model file for reading.
    FILE* lpGeomFile = fopen(xstrGeomFile, "rb");
    if(!lpGeomFile)
        return eFileReadError;

    // Determine the size of the geometry file and then copy it into the
    // cubit file.
    EErrorCode leReturn = eSuccess;
    try {
        UnsignedInt32 lintReadIndex, lintWriteIndex;
        // Locate the model's index in the write contents and make sure it is
        // not an FEA model.
        if(!FindModel(xintGeomModel, mpaWriteModels,
            mWriteContents.mintNumModels, lintWriteIndex))
            throw eNotFound;
        if(mpaWriteModels[lintWriteIndex].mintModelType == eFEModel)
            throw eNotFound;
        // Locate the model's index in the read contents, if possible, and
        // make sure the model has not already been written or deleted, and
        // then mark it written.
        if(FindModel(xintGeomModel, mpaReadModels,
            mReadContents.mintNumModels, lintReadIndex)) {
            if(mpaReadModelStat[lintReadIndex] != eStatNotWritten)
                throw eOrderError;
            mpaReadModelStat[lintReadIndex] = eStatWritten;
        }

        // Measure the length of the geometry model file and then copy it into
        // to cubit file.
        if(fseek(lpGeomFile, 0, SEEK_END))
            throw CCubitFile::eFileSeekError;
        long lintGeomLength = ftell(lpGeomFile);
        if(lintGeomLength == -1L)
            throw CCubitFile::eFileTellError;
        CopyModel(0, mpaWriteModels[lintWriteIndex].mintModelOffset,
            lintGeomLength, lpGeomFile, mpWriteFile);
        mpaWriteModels[lintWriteIndex].mintModelLength = lintGeomLength;
    }
    catch(EErrorCode xeErrorCode)  {  leReturn = xeErrorCode;  }
    catch(...)  {  leReturn = eUnknownError;  }

    fclose(lpGeomFile);
    return leReturn;
}

UnsignedInt32 CCubitFile::EndWriteGeomModel()
// The formal end to a geometry model write.
{
    return eSuccess;
}

UnsignedInt32 CCubitFile::BeginReadGeomModel(HModel xintGeomModel,
                                             FILE* xpGeomFile)
// Copy an geometry model out of the cubit file into an external geometry file.
{
    if(!mpReadFile)  return eFileReadError;
    if(!xpGeomFile)  return eFileWriteError;

	UnsignedInt32 lintStartFilePos = ftell(xpGeomFile);

    // Determine the size of the geometry file and then copy it into the
    // cubit file.
    EErrorCode leReturn = eSuccess;
    try {
        UnsignedInt32 lintReadIndex, lintWriteOffset;
        // Locate the model's index in the read contents and make sure it is
        // not an FEA model.
        if(!FindModel(xintGeomModel, mpaReadModels,
            mReadContents.mintNumModels, lintReadIndex))
            throw eNotFound;
        if(mpaReadModels[lintReadIndex].mintModelType == eFEModel)
            throw eNotFound;

        // Copy the geometry model file out of the cubit file.
        CopyModel(mpaReadModels[lintReadIndex].mintModelOffset,
            lintWriteOffset, mpaReadModels[lintReadIndex].mintModelLength,
            mpReadFile, xpGeomFile);
    }
    catch(EErrorCode xeErrorCode)  {  leReturn = xeErrorCode;  }
    catch(...)  {  leReturn = eUnknownError;  }

    fseek(xpGeomFile, lintStartFilePos, SEEK_SET);
    return leReturn;
}

UnsignedInt32 CCubitFile::EndReadGeomModel()
// The formal end to a geometry model read.
{
    return eSuccess;
}


UnsignedInt32 CCubitFile::BeginWriteAssemblyModel(HModel xintModelId,
                                                  FILE*& writeable_file)
{
    // Let's not return a FILE* until we know it's going to be OK to write to.
  writeable_file = NULL;
  
  if(!mpWriteFile)
    return eFileWriteError;
  
  EErrorCode leReturn = eSuccess;
  try
  {
    UnsignedInt32 lintReadIndex, lintWriteIndex;
      // Locate the model's index in the write contents and make sure it is
      // the right kind of model.
    if(!FindModel(xintModelId,
                  mpaWriteModels,
                  mWriteContents.mintNumModels,
                  lintWriteIndex))
      throw eNotFound;
    
    if(mpaWriteModels[lintWriteIndex].mintModelType != eAssemblyModel)
      throw eNotFound;
    
      // Locate the model's index in the read contents, if possible, and
      // make sure the model has not already been written or deleted, and
      // then mark it as being written.
    if(FindModel(xintModelId, mpaReadModels,
                 mReadContents.mintNumModels, lintReadIndex))
    {
      if(mpaReadModelStat[lintReadIndex] != eStatNotWritten)
        throw eOrderError;
      mpaReadModelStat[lintReadIndex] = eStatWriting;
    }
    
      // Move the FILE* to the correct write position so the FILE*
      // can be written to directly.  From CopyFile, it looks like the
      // correct location is the end of the file.
    if (fseek(mpWriteFile, 0, SEEK_END))
      throw eFileSeekError;
    
      // Save our current FILE* position so we can detect the size of
      // data written to the model.
    mpaWriteModels[lintWriteIndex].mintModelOffset = ftell(mpWriteFile);
    if ((int)mpaWriteModels[lintWriteIndex].mintModelOffset == -1L)
      throw eFileTellError;
    
      // set the FILE*
    writeable_file = mpWriteFile;
  }
  catch(EErrorCode xeErrorCode)
  {  leReturn = xeErrorCode;  }
  catch(...)
  {  leReturn = eUnknownError;  }
  
  return leReturn;
}

UnsignedInt32 CCubitFile::EndWriteAssemblyModel(HModel xintModelId)
{
  if (!mpWriteFile)
    return eFileWriteError;
  
    // Find the model's write record
  UnsignedInt32 lintWriteIndex;
  if(!FindModel(xintModelId,
                mpaWriteModels,
                mWriteContents.mintNumModels,
                lintWriteIndex))
    return eNotFound;

    // Make sure we saved the start position
  if (mpaWriteModels[lintWriteIndex].mintModelOffset == 0)
    return eOrderError;
  
    // Get the end of the FILE.
  if (fseek(mpWriteFile, 0, SEEK_END))
    return eFileSeekError;
  long cur_pos = ftell(mpWriteFile);
  if (cur_pos == -1L)
    return eFileTellError;
  
    // See how many bytes that is past our saved position.
  long size_in_bytes = cur_pos - mpaWriteModels[lintWriteIndex].mintModelOffset;
  
    // Locate the model's index in the read contents, if possible, and
    // make sure the model has not already been written or deleted, and
    // then mark it as written.
  UnsignedInt32 lintReadIndex;
  if(FindModel(xintModelId, mpaReadModels,
               mReadContents.mintNumModels, lintReadIndex))
  {
    if(mpaReadModelStat[lintReadIndex] != eStatWriting)
      throw eOrderError;
    mpaReadModelStat[lintReadIndex] = eStatWritten;
  }
  
    // Save the size to the model header.
  mpaWriteModels[lintWriteIndex].mintModelLength = size_in_bytes;
  
  return eSuccess;
}

// Returns a FILE* that can be read from.  The FILE* will NOT return EOF
// at the end of the data...the reader should not go GetReadModelLength()
// bytes beyond the current FILE* position.  You should copy the contents
// to another FILE if this is an issue for the reader.
UnsignedInt32 CCubitFile::BeginReadAssemblyModel(HModel xintModelId,
                                                 FILE*& f)
{
  f = NULL;
  
  if(!mpReadFile)
    return eFileReadError;
  
    // Determine the size of the geometry file and then copy it into the
    // cubit file.
  EErrorCode leReturn = eSuccess;
  try
  {
    UnsignedInt32 lintReadIndex;
      // Locate the model's index in the read contents and make sure it is
      // an assembly model.
    if(!FindModel(xintModelId, mpaReadModels,
                  mReadContents.mintNumModels, lintReadIndex))
      throw eNotFound;
    if(mpaReadModels[lintReadIndex].mintModelType != eAssemblyModel)
      throw eNotFound;
    
      // Set the read pointer to the correct location
    fseek(mpReadFile, mpaReadModels[lintReadIndex].mintModelOffset, SEEK_SET);
    
      // Set the FILE*
    f = mpReadFile;
  }
  catch(EErrorCode xeErrorCode)
  { leReturn = xeErrorCode; }
  catch(...)
  { leReturn = eUnknownError; }
  
  return leReturn;
}

UnsignedInt32 CCubitFile::EndReadAssemblyModel()
{
  return eSuccess;
}


///////////////////////////////////////////////////////////////////////////////
// FEA Model Read/Write
///////////////////////////////////////////////////////////////////////////////

UnsignedInt32 CCubitFile::BeginWriteFEModel(HModel xintFEModel,
                                            UnsignedInt32 xintGeomCount,
                                            UnsignedInt32 xintGroupCount,
                                            UnsignedInt32 xintBlockCount,
                                            UnsignedInt32 xintNodeSetCount,
                                            UnsignedInt32 xintSideSetCount)
// Prepares the file for the writing of a FEA model.  After this call there
// should be no other write oprations on the file not associated with the 
// writing of the FEA model until a matching call is made to EndWriteFEModel.
// RETURNS: 0 on success, other values indicate error code.
{
    try {
        if(!mpWriteFile)  throw eFileWriteError;
        if(mpWriteFEModel)  throw eOrderError;

        // Write the model table first to try to guarentee that it will be
        // at the beginning of the file.
        WriteModelTable();

        // Lookup the requested model in the model table and make sure it is
        // an FEA model.
        if(!FindModel(xintFEModel, mpaWriteModels, mWriteContents.mintNumModels,
            mintFEModelIndex))
            throw eNotFound;
        if(mpaWriteModels[mintFEModelIndex].mintModelType != eFEModel)
            throw eNotFound;

        // Create the object that manages FEA model file transactions and have
        // it begin its model write.
        mpWriteFEModel = new CFEModel();
        if(!mpWriteFEModel)
            throw eMemoryError;
        
        mpaWriteModels[mintFEModelIndex].mintModelOffset =
            mpWriteFEModel->InitWrite(mpWriteFile, xintGeomCount, xintGroupCount,
            xintBlockCount, xintNodeSetCount, xintSideSetCount);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return eUnknownError;  }
}

// Try to write the nodes for the passed geometry index/ID to the writable file
UnsignedInt32 CCubitFile::WriteNodes(UnsignedInt32 xintIndex,
                                     UnsignedInt32 xintGeomID,
                                     UnsignedInt32 xintNodeCount,
                                     UnsignedInt32 *xpaintNodeIDs,
                                     double *xpadblX,
                                     double *xpadblY,
                                     double *xpadblZ)
{
  try
  {
    if(!mpWriteFEModel)
      throw eOrderError;
    mpWriteFEModel->WriteNodes(xintIndex, xintGeomID,
                               xintNodeCount, xpaintNodeIDs,
                               xpadblX, xpadblY, xpadblZ);
  }
  catch(EErrorCode xeErrorCode)  {  return xeErrorCode;  }
  catch(...)  {  return eUnknownError;  }
  return eSuccess;
}

UnsignedInt32 CCubitFile::WriteElems(UnsignedInt32 xintIndex,
                                     UnsignedInt32 xintNumTypes,
                                     SElemData* xpaElemData)
// Try to write the element connectivity for the passed geometry index/ID
// to the writable file.
{
    try {
        if(!mpWriteFEModel)
            throw eOrderError;
        mpWriteFEModel->WriteElems(xintIndex, xintNumTypes, xpaElemData);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return eUnknownError;  }
}

UnsignedInt32 CCubitFile::WriteGroup(UnsignedInt32 xintIndex,
                                     UnsignedInt32 xintGroupID,
                                     UnsignedInt32 xintGroupType,
                                     ConstCharPtr xpachrGroupName,
                                     UnsignedInt32 xintNumTypes,
                                     SGroupData* xpaGroupData)
// Try to write the membership of the passed group index to the writable file.
// RETURNS: 0 on success, other values indicate error code.
{
    try {
        if(!mpWriteFEModel)
            throw eOrderError;
        mpWriteFEModel->WriteGroup(xintIndex, xintGroupID,
            xintGroupType, xpachrGroupName, xintNumTypes, xpaGroupData);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return eUnknownError;  }
}

UnsignedInt32 CCubitFile::WriteBlock(UnsignedInt32 xintIndex,
                                     UnsignedInt32 xintBlockID,
                                     UnsignedInt32 xintBlockType,
                                     UnsignedInt32 xintBlockColor,
                                     UnsignedInt32 xintMixedElemType,
                                     UnsignedInt32 xintDefPyramidType,
                                     UnsignedInt32 xintMaterialID,
                                     UnsignedInt32 xintBlockDimension,
                                     UnsignedInt32 xintNumTypes,
                                     SBlockData* xpaBlockData,
                                     UnsignedInt32 xintAttributeOrder,
                                     double* xpadblAttributes)
// Try to write the membership of the passed group index to the writable file.
// RETURNS: 0 on success, other values indicate error code.
{
    try {
        if(!mpWriteFEModel)
            throw eOrderError;
        mpWriteFEModel->WriteBlock(xintIndex, xintBlockID, xintBlockType,
            xintBlockColor, xintMixedElemType, xintDefPyramidType,
            xintMaterialID, xintBlockDimension, xintNumTypes, xpaBlockData,
            xintAttributeOrder, xpadblAttributes);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return eUnknownError;  }
}

UnsignedInt32 CCubitFile::WriteNodeSet(UnsignedInt32 xintIndex,
                                       UnsignedInt32 xintNodeSetID,
                                       UnsignedInt32 xintColor,
                                       UnsignedInt32 xintPointSymbol,
                                       UnsignedInt32 xintNumTypes,
                                       SNodeSetData* xpaNodeSetData)
// Try to write the membership of the passed group index to the writable file.
// RETURNS: 0 on success, other values indicate error code.
{
    try {
        if(!mpWriteFEModel)
            throw eOrderError;
        mpWriteFEModel->WriteNodeSet(xintIndex, xintNodeSetID, xintColor,
            xintPointSymbol, xintNumTypes, xpaNodeSetData);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return eUnknownError;  }
}

UnsignedInt32 CCubitFile::WriteSideSet_11(UnsignedInt32 xintIndex,
                                       UnsignedInt32 xintSideSetID,
                                       UnsignedInt32 xintColor,
                                       UnsignedInt32 xintUseShells,
                                       UnsignedInt32 xintNumTypes,
                                       SSideSetData_11* xpaSideSetData,
                                       UnsignedInt32 xintNumDistFact,
                                       double* xpadblDistribution)
// Try to write the membership of the passed group index to the writable file.
// RETURNS: 0 on success, other values indicate error code.
{
    try {
        if(!mpWriteFEModel)
            throw eOrderError;
        mpWriteFEModel->WriteSideSet_11(xintIndex, xintSideSetID, xintColor,
            xintUseShells, xintNumTypes, xpaSideSetData, xintNumDistFact,
            xpadblDistribution);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return eUnknownError;  }
}

UnsignedInt32 CCubitFile::EndWriteFEModel()
// Designate the end of writing a FEA model and free any memory allocated for
// the write operations and update the file contents header.
// RETURNS: 0 on success, other values indicate error code.
{
    try {
        if(!mpWriteFEModel)
            throw eOrderError;
        mpaWriteModels[mintFEModelIndex].mintModelLength =
            mpWriteFEModel->EndWrite();
        delete mpWriteFEModel;
        mpWriteFEModel = NULL;

        UnsignedInt32 lintReadIndex;
        if(FindModel(mpaWriteModels[mintFEModelIndex].mintModelHandle,
            mpaReadModels, mReadContents.mintNumModels, lintReadIndex))
            mpaReadModelStat[lintReadIndex] = eStatWritten;

        return eSuccess;
    }
    catch(EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return eUnknownError;  }
}

//-----------------------------------------------------------------------------

UnsignedInt32 CCubitFile::BeginReadFEModel(HModel xintFEModel,
                                           UnsignedInt32& xintGeomCount,
                                           UnsignedInt32& xintGroupCount,
                                           UnsignedInt32& xintBlockCount,
                                           UnsignedInt32& xintNodeSetCount,
                                           UnsignedInt32& xintSideSetCount)
// Prepares the file for the reading of a FEA model.  After this call there
// should be no other read oprations on the file not associated with the 
// reading of the FEA model until a matching call is made to EndReadFEModel.
// RETURNS: 0 on success, other values indicate error code.
{
    try {
        if(!mpReadFile)  throw eFileReadError;
        if(mpReadFEModel)  throw eOrderError;

        UnsignedInt32 lintIndex;
        if(!FindModel(xintFEModel, mpaReadModels, mReadContents.mintNumModels,
            lintIndex))
            throw eNotFound;
        if(mpaReadModels[lintIndex].mintModelType != eFEModel)
            throw eNotFound;

        mpReadFEModel = new CFEModel();
        if(!mpReadFEModel)  throw eMemoryError;
        
        mpReadFEModel->InitRead(mpReadFile, mpaReadModels[lintIndex].mintModelOffset,
            xintGeomCount, xintGroupCount,
            xintBlockCount, xintNodeSetCount, xintSideSetCount);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return eUnknownError;  }
}

UnsignedInt32 CCubitFile::ReadNodes(UnsignedInt32 xintIndex,
									UnsignedInt32& xintGeomID,
                                    UnsignedInt32& xintNodeCount,
                                    UnsignedInt32*& xpaintNodeIDs,
                                    double*& xpadblX,
                                    double*& xpadblY,
                                    double*& xpadblZ)
// Try to read the nodes for the passed geometry index/ID from the
// read-only file.
// RETURNS: 0 on success, other values indicate error code.
{
    UnsignedInt32 lintErrorCode;
    try {
        if(!mpReadFEModel)
            throw eOrderError;
        mpReadFEModel->ReadNodes(xintIndex, xintGeomID,
            xintNodeCount, xpaintNodeIDs, xpadblX, xpadblY, xpadblZ);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode) {
        lintErrorCode = xeErrorCode;
    }
    catch(...) {
        lintErrorCode = eUnknownError;
    }
    xintNodeCount = 0;
    xpaintNodeIDs = NULL;
    xpadblX = xpadblY = xpadblZ = NULL;
    return lintErrorCode;
}

UnsignedInt32 CCubitFile::ReadElems(UnsignedInt32 xintIndex,
									UnsignedInt32& xintGeomID,
                                    UnsignedInt32& xintNumTypes,
                                    SElemData*& xpaElemData)
// Try to read the element connectivity for the passed geometry index/ID
// from the read-only file.
// RETURNS: 0 on success, other values indicate error code.
{
    UnsignedInt32 lintErrorCode;
    try {
        if(!mpReadFEModel)
            throw eOrderError;
        mpReadFEModel->ReadElems(xintIndex, xintGeomID,
			xintNumTypes, xpaElemData);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode) {
        lintErrorCode = xeErrorCode;
    }
    catch(...) {
        lintErrorCode = eUnknownError;
    }
    xintNumTypes = 0;
    xpaElemData = NULL;
    return lintErrorCode;
}

UnsignedInt32 CCubitFile::ReadGroupIdentity(UnsignedInt32 xintIndex,
                                            UnsignedInt32& xintGroupID,
                                            UnsignedInt32& xintGroupType,
                                            ConstCharPtr& xpachrGroupName)
// Try to read the identity (ID, name, type) of the passed group index from the
// read-only file.
// RETURNS: 0 on success, other values indicate error code.
{
    UnsignedInt32 lintErrorCode;
    try {
        if(!mpReadFEModel)
            throw eOrderError;
        mpReadFEModel->ReadGroupIdentity(xintIndex,
            xintGroupID, xintGroupType, xpachrGroupName);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode) {
        lintErrorCode = xeErrorCode;
    }
    catch(...) {
        lintErrorCode = eUnknownError;
    }
    xintGroupID = xintGroupType = 0;
    xpachrGroupName = NULL;
    return lintErrorCode;
}

UnsignedInt32 CCubitFile::ReadGroupMembers(UnsignedInt32 xintIndex,
                                           UnsignedInt32& xintNumTypes,
                                           SGroupData*& xpaGroupData)
// Try to read the membership of the passed group index from the read-only file.
// Note: The 
// RETURNS: 0 on success, other values indicate error code.
{
    UnsignedInt32 lintErrorCode;
    try {
        if(!mpReadFEModel)
            throw eOrderError;
        mpReadFEModel->ReadGroupMembers(xintIndex, xintNumTypes, xpaGroupData);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode) {
        lintErrorCode = xeErrorCode;
    }
    catch(...) {
        lintErrorCode = eUnknownError;
    }
    xintNumTypes = 0;
    xpaGroupData = NULL;
    return lintErrorCode;
}

UnsignedInt32 CCubitFile::ReadBlock(UnsignedInt32 xintIndex,
                                    UnsignedInt32& xintBlockID,
                                    UnsignedInt32& xintBlockType,
                                    UnsignedInt32& xintBlockColor,
                                    UnsignedInt32& xintMixedElemType,
                                    UnsignedInt32& xintDefPyramidType,
                                    UnsignedInt32& xintMaterialID,
                                    UnsignedInt32& xintBlockDimension,
                                    UnsignedInt32& xintNumTypes,
                                    SBlockData*& xpaBlockData,
                                    UnsignedInt32& xintAttributeOrder,
                                    double*& xpadblAttributes)
// Try to read the membership of the passed block index from the read-only file.
// Note: The 
// RETURNS: 0 on success, other values indicate error code.
{
    UnsignedInt32 lintErrorCode;
    try {
        if(!mpReadFEModel)
            throw eOrderError;
        mpReadFEModel->ReadBlock(xintIndex, xintBlockID,
            xintBlockType, xintBlockColor, xintMixedElemType,
            xintDefPyramidType, xintMaterialID, xintBlockDimension,
            xintNumTypes, xpaBlockData, xintAttributeOrder, xpadblAttributes);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode) {
        lintErrorCode = xeErrorCode;
    }
    catch(...) {
        lintErrorCode = eUnknownError;
    }
    xintBlockID = xintBlockType = xintBlockColor = xintMixedElemType =
        xintDefPyramidType = xintMaterialID = xintNumTypes =
        xintAttributeOrder = 0;
    xpaBlockData = NULL;
    xpadblAttributes = NULL;
    return lintErrorCode;
}

UnsignedInt32 CCubitFile::ReadNodeSet(UnsignedInt32 xintIndex,
                                      UnsignedInt32& xintNodeSetID,
                                      UnsignedInt32& xintColor,
                                      UnsignedInt32& xintPointSymbol,
                                      UnsignedInt32& xintNumTypes,
                                      SNodeSetData*& xpaNodeSetData)
// Try to read the membership of the passed node set index from the read-only file.
// Note: The 
// RETURNS: 0 on success, other values indicate error code.
{
    UnsignedInt32 lintErrorCode;
    try {
        if(!mpReadFEModel)
            throw eOrderError;
        mpReadFEModel->ReadNodeSet(xintIndex, xintNodeSetID,
            xintColor, xintPointSymbol, xintNumTypes, xpaNodeSetData);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode) {
        lintErrorCode = xeErrorCode;
    }
    catch(...) {
        lintErrorCode = eUnknownError;
    }
    xintNodeSetID = xintColor = xintPointSymbol = xintNumTypes = 0;
    xpaNodeSetData = NULL;
    return lintErrorCode;
}

UnsignedInt32 CCubitFile::ReadSideSet_10(UnsignedInt32 xintIndex,
                                      UnsignedInt32& xintSideSetID,
                                      UnsignedInt32& xintColor,
                                      UnsignedInt32& xintUseShells,
                                      UnsignedInt32& xintNumTypes,
                                      SSideSetData_10*& xpaSideSetData,
                                      UnsignedInt32& xintNumDistFact,
                                      double*& xpadblDistribution)
// Try to read the membership of the passed side set index from the read-only file.
// Note: The 
// RETURNS: 0 on success, other values indicate error code.
{
    UnsignedInt32 lintErrorCode;
    try {
        if(!mpReadFEModel)
            throw eOrderError;
        mpReadFEModel->ReadSideSet_10(xintIndex, xintSideSetID,
            xintColor, xintUseShells, xintNumTypes, xpaSideSetData,
            xintNumDistFact, xpadblDistribution);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode) {
        lintErrorCode = xeErrorCode;
    }
    catch(...) {
        lintErrorCode = eUnknownError;
    }
    xintSideSetID = xintColor = xintNumTypes = xintNumDistFact = 0;
    xpaSideSetData = NULL;
    xpadblDistribution = NULL;
    return lintErrorCode;
}

UnsignedInt32 CCubitFile::ReadSideSet_11(UnsignedInt32 xintIndex,
                                      UnsignedInt32& xintSideSetID,
                                      UnsignedInt32& xintColor,
                                      UnsignedInt32& xintUseShells,
                                      UnsignedInt32& xintNumTypes,
                                      SSideSetData_11*& xpaSideSetData,
                                      UnsignedInt32& xintNumDistFact,
                                      double*& xpadblDistribution)
// Try to read the membership of the passed side set index from the read-only file.
// Note: The 
// RETURNS: 0 on success, other values indicate error code.
{
    UnsignedInt32 lintErrorCode;
    try {
        if(!mpReadFEModel)
            throw eOrderError;
        mpReadFEModel->ReadSideSet_11(xintIndex, xintSideSetID,
            xintColor, xintUseShells, xintNumTypes, xpaSideSetData,
            xintNumDistFact, xpadblDistribution);
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode) {
        lintErrorCode = xeErrorCode;
    }
    catch(...) {
        lintErrorCode = eUnknownError;
    }
    xintSideSetID = xintColor = xintNumTypes = xintNumDistFact = 0;
    xpaSideSetData = NULL;
    xpadblDistribution = NULL;
    return lintErrorCode;
}

UnsignedInt32 CCubitFile::EndReadFEModel()
// Designate the end of reading a FEA model and free any memory allocated for
// the read operations.
// RETURNS: 0 on success, other values indicate error code.
{
    try {
        if(!mpReadFEModel)
            return eOrderError;
        mpReadFEModel->EndRead();
        delete mpReadFEModel;
        mpReadFEModel = NULL;
        return eSuccess;
    }
    catch(EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return eUnknownError;  }
}


///////////////////////////////////////////////////////////////////////////////
// Meta-data accessors
///////////////////////////////////////////////////////////////////////////////

UnsignedInt32 CCubitFile::GetReadMetaData(EMetaDataOwnerType xeType,
                                          CMetaData*& xpMetaData)
// Return the requested meta-data object for the read file if possible.
{
    xpMetaData = NULL;
    if(xeType == eModelMetaData)
        xpMetaData = mpMetaData;
    else if(mpReadFEModel) {
        switch(xeType) {
        case eGeomMetaData:
            xpMetaData = &mpReadFEModel->GetGeomMetaData();
            break;
        case eNodeMetaData:
            xpMetaData = &mpReadFEModel->GetNodeMetaData();
            break;
        case eElemMetaData:
            xpMetaData = &mpReadFEModel->GetElemMetaData();
            break;
        case eGroupMetaData:
            xpMetaData = &mpReadFEModel->GetGroupMetaData();
            break;
        case eBlockMetaData:
            xpMetaData = &mpReadFEModel->GetBlockMetaData();
            break;
        case eNodeSetMetaData:
            xpMetaData = &mpReadFEModel->GetNodeSetMetaData();
            break;
        case eSideSetMetaData:
            xpMetaData = &mpReadFEModel->GetSideSetMetaData();
            break;
          default:
            break;
        }
    }
    else
        return eOrderError;
    return xpMetaData ? eSuccess : eNotFound;
}

UnsignedInt32 CCubitFile::GetWriteMetaData(EMetaDataOwnerType xeType,
                                           CMetaData*& xpMetaData)
// Return the requested meta-data object for the write file if possible.
{
    xpMetaData = NULL;
    if(xeType == eModelMetaData)
        xpMetaData = mpMetaData;
    else if(mpWriteFEModel) {
        switch(xeType) {
        case eGeomMetaData:
            xpMetaData = &mpWriteFEModel->GetGeomMetaData();
            break;
        case eNodeMetaData:
            xpMetaData = &mpWriteFEModel->GetNodeMetaData();
            break;
        case eElemMetaData:
            xpMetaData = &mpWriteFEModel->GetElemMetaData();
            break;
        case eGroupMetaData:
            xpMetaData = &mpWriteFEModel->GetGroupMetaData();
            break;
        case eBlockMetaData:
            xpMetaData = &mpWriteFEModel->GetBlockMetaData();
            break;
        case eNodeSetMetaData:
            xpMetaData = &mpWriteFEModel->GetNodeSetMetaData();
            break;
        case eSideSetMetaData:
            xpMetaData = &mpWriteFEModel->GetSideSetMetaData();
            break;
          default:
            break;
        }
    }
    else
        return eOrderError;
    return xpMetaData ? eSuccess : eNotFound;
}
