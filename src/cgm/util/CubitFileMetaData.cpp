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

 
 Filename      : CubitFileMetaData.cpp

 Purpose       : Implements meta-data (descriptions) in the Cubit file.
           
 Special Notes :

 Creator       : Will A. Helden

 Creation Date : 02/15/02

 Owner         : Will A. Helden

*******************************************************************************/

#include "CubitFileMetaData.hpp"
#include "CubitFileIOWrapper.hpp"
#include <string.h>

using namespace NCubitFile;

///////////////////////////////////////////////////////////////////////////////
// CMetaData -
///////////////////////////////////////////////////////////////////////////////

CMetaData::CMetaData(UnsignedInt32 xintGrowSizeBy)
{
    mHeader.mintMetaDataSchema = 0;
    mHeader.mintMetaDataCompress = 0;
    mHeader.mintMetaDataCount = 0;

    mintGrowSizeBy = xintGrowSizeBy ? xintGrowSizeBy : 64;
    mintBufferSize = 0;
    mpaEntries = NULL;

    // Initialize the return buffer - allocated structures needed to return
    // the results of meta data queries will be owned and freed by this object.
    mSearchBuff.mintNumOwners = mSearchBuff.mintNumNames =
        mSearchBuff.mintNumMetaData = 0;
    mSearchBuff.mpaintOwners = NULL;
    mSearchBuff.mpapachrNames = NULL;
    mSearchBuff.mintMetaDataType = eUnsignedInt32;
    mSearchBuff.mpaValues = NULL;
    mSearchBuff.mpaintNumValues = NULL;
}

CMetaData::~CMetaData()
{
    // Free the meta data.
    FreeAll();

    // Free the search return buffer.
    if(mSearchBuff.mpaintOwners)
        delete [] mSearchBuff.mpaintOwners;
    if(mSearchBuff.mpapachrNames)
        delete [] mSearchBuff.mpapachrNames;
    if(mSearchBuff.mpaValues) {
        switch(mSearchBuff.mintMetaDataType) {
        case eUnsignedInt32:
            delete [] (UnsignedInt32*)mSearchBuff.mpaValues;
            break;
        case eString:
            delete [] (ConstCharPtr*)mSearchBuff.mpaValues;
            break;
        case eDouble:
            delete [] (double*)mSearchBuff.mpaValues;
            break;
        case eUnsignedInt32Array:
            delete [] (ConstUnsignedInt32Ptr*)mSearchBuff.mpaValues;
            break;
        case eDoubleArray:
            delete [] (ConstDoublePtr*)mSearchBuff.mpaValues;
            break;
        }
    }
    if(mSearchBuff.mpaintNumValues)
        delete [] mSearchBuff.mpaintNumValues;
}

void CMetaData::GrowSize()
// Grows the meta-data storage structure by the predetermined growth size.
{
    // Try to allocate a new storage structure.
    UnsignedInt32 lintNewBufferSize = mintBufferSize + mintGrowSizeBy;
    SCubitFileMetaDataEntry* lpaNewEntries =
        new SCubitFileMetaDataEntry[lintNewBufferSize];
    if(!lpaNewEntries)
        throw CCubitFile::eMemoryError;

    // Copy the existing data into the front new structure.
    if(mpaEntries) {
        memcpy(lpaNewEntries, mpaEntries,
            mintBufferSize * sizeof(SCubitFileMetaDataEntry));
    }
    // Initialize the remaining part of the new structure.
    for(UnsignedInt32 lintEntry = mintBufferSize; lintEntry <
        lintNewBufferSize; lintEntry++) {
        lpaNewEntries[lintEntry].mintMetaDataOwner = 0;
        lpaNewEntries[lintEntry].mintMetaDataType = eUnsignedInt32;
        lpaNewEntries[lintEntry].mpachrName = NULL;
        lpaNewEntries[lintEntry].muValue.mint = 0;
        lpaNewEntries[lintEntry].mintNumValues = 0;
    }
    // Replace the old structure with the new one and free the old one.
    if(mpaEntries)
        delete [] mpaEntries;
    mintBufferSize = lintNewBufferSize;
    mpaEntries = lpaNewEntries;
}

UnsignedInt32 CMetaData::AddEntry(UnsignedInt32 xintOwner,
                                  ConstCharPtr xpachrName)
// Creates a valueless meta-data entry for the passed id and name.
// Private function, returns the index of the new entry, throws exceptions.
{
    if(mHeader.mintMetaDataCount == mintBufferSize)
        GrowSize();
    mpaEntries[mHeader.mintMetaDataCount].mintMetaDataOwner = xintOwner;
    UnsignedInt32 lintNameLength = strlen(xpachrName);
    mpaEntries[mHeader.mintMetaDataCount].mpachrName =
        new char[lintNameLength + 1];
    if(!mpaEntries[mHeader.mintMetaDataCount].mpachrName)
        throw CCubitFile::eMemoryError;
    strcpy(mpaEntries[mHeader.mintMetaDataCount].mpachrName, xpachrName);
    return mHeader.mintMetaDataCount++;
}

UnsignedInt32 CMetaData::Find(UnsignedInt32 xintOwner,
                              ConstCharPtr xpachrName,
                              UnsignedInt32& xintIndex)
// Search for the index of a piece of data with the declared owner and name.
// Private function, returns true/false, no exceptions thrown.
{
    if(!mHeader.mintMetaDataCount)  return 0;  // false, not found
    for(UnsignedInt32 lintEntry = 0; lintEntry <
        mHeader.mintMetaDataCount; lintEntry++) {
        if(mpaEntries[lintEntry].mintMetaDataOwner == xintOwner) {
            if(!strcmp(xpachrName, mpaEntries[lintEntry].mpachrName)) {
                xintIndex = lintEntry;
                return 1;  // true, found
            }
        }
    }
    return 0;  // false, not found
}


///////////////////////////////////////////////////////////////////////////////
// Get/Set Value Methods
///////////////////////////////////////////////////////////////////////////////

UnsignedInt32 CMetaData::GetValue(UnsignedInt32 xintOwner,
                                  ConstCharPtr xpachrName,
                                  UnsignedInt32& xintValue)
// Retrieve integer valued meta-data by id and name.
// Public function, returns CCubitFile::EErrorCode, no exceptions thrown.
{
    UnsignedInt32 lintEntry;
    if(Find(xintOwner, xpachrName, lintEntry)) {
        if(mpaEntries[lintEntry].mintMetaDataType == eUnsignedInt32) {
            xintValue = mpaEntries[lintEntry].muValue.mint;
            return CCubitFile::eSuccess;
        }
    }
    xintValue = 0;
    return CCubitFile::eNotFound;
}

UnsignedInt32 CMetaData::GetValue(UnsignedInt32 xintOwner,
                                  ConstCharPtr xpachrName,
                                  ConstCharPtr& xpachrValue)
// Retrieve string valued meta-data by id and name.
// Public function, returns CCubitFile::EErrorCode, no exceptions thrown.
{
    UnsignedInt32 lintEntry;
    if(Find(xintOwner, xpachrName, lintEntry)) {
        if(mpaEntries[lintEntry].mintMetaDataType == eString) {
            xpachrValue = mpaEntries[lintEntry].muValue.mpachr;
            return CCubitFile::eSuccess;
        }
    }
    xpachrValue = NULL;
    return CCubitFile::eNotFound;
}

UnsignedInt32 CMetaData::GetValue(UnsignedInt32 xintOwner,
                                  ConstCharPtr xpachrName,
                                  double& xdblValue)
// Retrieve double valued meta-data by id and name.
// Public function, returns CCubitFile::EErrorCode, no exceptions thrown.
{
    UnsignedInt32 lintEntry;
    if(Find(xintOwner, xpachrName, lintEntry)) {
        if(mpaEntries[lintEntry].mintMetaDataType == eDouble) {
            xdblValue = mpaEntries[lintEntry].muValue.mdbl;
            return CCubitFile::eSuccess;
        }
    }
    xdblValue = 0.0;
    return CCubitFile::eNotFound;
}

UnsignedInt32 CMetaData::GetValue(UnsignedInt32 xintOwner,
                                  ConstCharPtr xpachrName,
                                  const UnsignedInt32*& xpaintValue,
                                  UnsignedInt32& xintNumValues)
// Retrieve integer array valued meta-data by id and name.
// Public function, returns CCubitFile::EErrorCode, no exceptions thrown.
{
    UnsignedInt32 lintEntry;
    if(Find(xintOwner, xpachrName, lintEntry)) {
        if(mpaEntries[lintEntry].mintMetaDataType == eUnsignedInt32Array) {
            xpaintValue = mpaEntries[lintEntry].muValue.mpaint;
            xintNumValues = mpaEntries[lintEntry].mintNumValues;
            return CCubitFile::eSuccess;
        }
    }
    xpaintValue = NULL;
    xintNumValues = 0;
    return CCubitFile::eNotFound;
}

UnsignedInt32 CMetaData::GetValue(UnsignedInt32 xintOwner,
                                  ConstCharPtr xpachrName,
                                  const double*& xpadblValue,
                                  UnsignedInt32& xintNumValues)
// Retrieve double array valued meta-data by id and name.
// Public function, returns CCubitFile::EErrorCode, no exceptions thrown.
{
    UnsignedInt32 lintEntry;
    if(Find(xintOwner, xpachrName, lintEntry)) {
        if(mpaEntries[lintEntry].mintMetaDataType == eDoubleArray) {
            xpadblValue = mpaEntries[lintEntry].muValue.mpadbl;
            xintNumValues = mpaEntries[lintEntry].mintNumValues;
            return CCubitFile::eSuccess;
        }
    }
    xpadblValue = NULL;
    xintNumValues = 0;
    return CCubitFile::eNotFound;
}

//-----------------------------------------------------------------------------

UnsignedInt32 CMetaData::SetValue(UnsignedInt32 xintOwner,
                                  ConstCharPtr xpachrName,
                                  UnsignedInt32 xintValue)
// Assign integer valued meta-data to an id and name.
// Public function, returns CCubitFile::EErrorCode, no exceptions thrown.
{
    if(!xpachrName)
        return CCubitFile::ePassedNullPointer;

    try {
        UnsignedInt32 lintEntry;
        if(Find(xintOwner, xpachrName, lintEntry))
            FreeValueAt(lintEntry);
        else
            lintEntry = AddEntry(xintOwner, xpachrName);
        mpaEntries[lintEntry].mintMetaDataType = eUnsignedInt32;
        mpaEntries[lintEntry].muValue.mint = xintValue;
        return CCubitFile::eSuccess;
    }
    catch(CCubitFile::EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return CCubitFile::eUnknownError;  }
}

UnsignedInt32 CMetaData::SetValue(UnsignedInt32 xintOwner,
                                  ConstCharPtr xpachrName,
                                  ConstCharPtr xpachrValue)
// Assign string valued meta-data to an id and name.
// Public function, returns CCubitFile::EErrorCode, no exceptions thrown.
{
    if(!xpachrName)
        throw CCubitFile::ePassedNullPointer;

    try {
        UnsignedInt32 lintEntry, lintLength;
        if(Find(xintOwner, xpachrName, lintEntry))
            FreeValueAt(lintEntry);
        else
            lintEntry = AddEntry(xintOwner, xpachrName);

        mpaEntries[lintEntry].mintMetaDataType = eString;
        lintLength = xpachrValue ? strlen(xpachrValue) : 0;
        if(lintLength) {
            mpaEntries[lintEntry].muValue.mpachr = new char[lintLength + 1];
            if(!mpaEntries[lintEntry].muValue.mpachr)
                throw CCubitFile::eMemoryError;
            strcpy(mpaEntries[lintEntry].muValue.mpachr, xpachrValue);
        }
        else
            mpaEntries[lintEntry].muValue.mpachr = NULL;

        return CCubitFile::eSuccess;
    }
    catch(CCubitFile::EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return CCubitFile::eUnknownError;  }
}

UnsignedInt32 CMetaData::SetValue(UnsignedInt32 xintOwner,
                                  ConstCharPtr xpachrName, double xdblValue)
// Assign double valued meta-data to an id and name.
// Public function, returns CCubitFile::EErrorCode, no exceptions thrown.
{
    if(!xpachrName)
        throw CCubitFile::ePassedNullPointer;

    try {
        UnsignedInt32 lintEntry;
        if(Find(xintOwner, xpachrName, lintEntry))
            FreeValueAt(lintEntry);
        else
            lintEntry = AddEntry(xintOwner, xpachrName);
        mpaEntries[lintEntry].mintMetaDataType = eDouble;
        mpaEntries[lintEntry].muValue.mdbl = xdblValue;
        return CCubitFile::eSuccess;
    }
    catch(CCubitFile::EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return CCubitFile::eUnknownError;  }
}

UnsignedInt32 CMetaData::SetValue(UnsignedInt32 xintOwner,
                                  ConstCharPtr xpachrName,
                                  const UnsignedInt32* xpaintValue,
                                  UnsignedInt32 xintNumValues)
// Assign integer array valued meta-data to an id and name.
// Public function, returns CCubitFile::EErrorCode, no exceptions thrown.
{
    if(!xpachrName)
        throw CCubitFile::ePassedNullPointer;

    try {
        UnsignedInt32 lintEntry;
        if(Find(xintOwner, xpachrName, lintEntry))
            FreeValueAt(lintEntry);
        else
            lintEntry = AddEntry(xintOwner, xpachrName);

        mpaEntries[lintEntry].mintMetaDataType = eUnsignedInt32Array;
        if(xpaintValue && xintNumValues) {
            mpaEntries[lintEntry].muValue.mpaint =
                new UnsignedInt32[xintNumValues];
            if(!mpaEntries[lintEntry].muValue.mpaint)
                throw CCubitFile::eMemoryError;
            memcpy(mpaEntries[lintEntry].muValue.mpaint, xpaintValue,
                xintNumValues * sizeof(UnsignedInt32));
            mpaEntries[lintEntry].mintNumValues = xintNumValues;
        }
        else {
            mpaEntries[lintEntry].muValue.mpaint = NULL;
            mpaEntries[lintEntry].mintNumValues = 0;
        }

        return CCubitFile::eSuccess;
    }
    catch(CCubitFile::EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return CCubitFile::eUnknownError;  }
}

UnsignedInt32 CMetaData::SetValue(UnsignedInt32 xintOwner,
                                  ConstCharPtr xpachrName,
                                  const double* xpadblValue,
                                  UnsignedInt32 xintNumValues)
// Assign double array valued meta-data to an id and name.
// Public function, returns CCubitFile::EErrorCode, no exceptions thrown.
{
    if(!xpachrName)
        throw CCubitFile::ePassedNullPointer;

    try {
        UnsignedInt32 lintEntry;
        if(Find(xintOwner, xpachrName, lintEntry))
            FreeValueAt(lintEntry);
        else
            lintEntry = AddEntry(xintOwner, xpachrName);

        mpaEntries[lintEntry].mintMetaDataType = eDoubleArray;
        if(xpadblValue && xintNumValues) {
            mpaEntries[lintEntry].muValue.mpadbl =
                new double[xintNumValues];
            if(!mpaEntries[lintEntry].muValue.mpadbl)
                throw CCubitFile::eMemoryError;
            memcpy(mpaEntries[lintEntry].muValue.mpadbl, xpadblValue,
                xintNumValues * sizeof(double));
            mpaEntries[lintEntry].mintNumValues = xintNumValues;
        }
        else {
            mpaEntries[lintEntry].muValue.mpadbl = NULL;
            mpaEntries[lintEntry].mintNumValues = 0;
        }

        return CCubitFile::eSuccess;
    }
    catch(CCubitFile::EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return CCubitFile::eUnknownError;  }
}



///////////////////////////////////////////////////////////////////////////////
// Search Methods
///////////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

UnsignedInt32 CMetaData::GetMetaDataAll(UnsignedInt32& xintNumFound,
                                        const ConstCharPtr*& xpapachrName,
                                        const UnsignedInt32*& xpaintOwner,
                                        const UnsignedInt32*& xpaintValue)
// Search the meta-data for all entries of type integer.
{
    xintNumFound = 0;
    xpapachrName = NULL;
    xpaintOwner = NULL;
    xpaintValue = NULL;
    try {  RebuildSearchBuff(eUnsignedInt32);  }
    catch(CCubitFile::EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return CCubitFile::eUnknownError;  }

    if((xintNumFound = SearchAll(eUnsignedInt32)) > 0) {
        xpapachrName = mSearchBuff.mpapachrNames;
        xpaintOwner = mSearchBuff.mpaintOwners;
        xpaintValue = (UnsignedInt32*)mSearchBuff.mpaValues;
    }
    return CCubitFile::eSuccess;
}

UnsignedInt32 CMetaData::GetMetaDataAll(UnsignedInt32& xintNumFound,
                                        const ConstCharPtr*& xpapachrName,
                                        const UnsignedInt32*& xpaintOwner,
                                        const ConstCharPtr*& xpapachrValue)
// Search the meta-data for all entries of type string.
{
    xintNumFound = 0;
    xpapachrName = NULL;
    xpaintOwner = NULL;
    xpapachrValue = NULL;
    try {  RebuildSearchBuff(eString);  }
    catch(CCubitFile::EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return CCubitFile::eUnknownError;  }

    if((xintNumFound = SearchAll(eString)) > 0) {
        xpapachrName = mSearchBuff.mpapachrNames;
        xpaintOwner = mSearchBuff.mpaintOwners;
        xpapachrValue = (ConstCharPtr*)mSearchBuff.mpaValues;
    }
    return CCubitFile::eSuccess;
}

UnsignedInt32 CMetaData::GetMetaDataAll(UnsignedInt32& xintNumFound,
                                        const ConstCharPtr*& xpapachrName,
                                        const UnsignedInt32*& xpaintOwner,
                                        const double*& xpadblValue)
// Search the meta-data for all entries of type double.
{
    xintNumFound = 0;
    xpapachrName = NULL;
    xpaintOwner = NULL;
    xpadblValue = NULL;
    try {  RebuildSearchBuff(eDouble);  }
    catch(CCubitFile::EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return CCubitFile::eUnknownError;  }

    if((xintNumFound = SearchAll(eDouble)) > 0) {
        xpapachrName = mSearchBuff.mpapachrNames;
        xpaintOwner = mSearchBuff.mpaintOwners;
        xpadblValue = (double*)mSearchBuff.mpaValues;
    }
    return CCubitFile::eSuccess;
}

UnsignedInt32 CMetaData::GetMetaDataAll(UnsignedInt32& xintNumFound,
                                        const ConstCharPtr*& xpapachrName,
                                        const UnsignedInt32*& xpaintOwner,
                                        const ConstUnsignedInt32Ptr*& xpapaintValue,
                                        const UnsignedInt32*& xpaintNumValues)
// Search the meta-data for all entries of type integer array.
{
    xintNumFound = 0;
    xpapachrName = NULL;
    xpaintOwner = NULL;
    xpapaintValue = NULL;
    xpaintNumValues = NULL;
    try {  RebuildSearchBuff(eUnsignedInt32Array);  }
    catch(CCubitFile::EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return CCubitFile::eUnknownError;  }

    if((xintNumFound = SearchAll(eUnsignedInt32Array)) > 0) {
        xpapachrName = mSearchBuff.mpapachrNames;
        xpaintOwner = mSearchBuff.mpaintOwners;
        xpapaintValue = (ConstUnsignedInt32Ptr*)mSearchBuff.mpaValues;
        xpaintNumValues = mSearchBuff.mpaintNumValues;
    }
    return CCubitFile::eSuccess;
}

UnsignedInt32 CMetaData::GetMetaDataAll(UnsignedInt32& xintNumFound,
                                        const ConstCharPtr*& xpapachrName,
                                        const UnsignedInt32*& xpaintOwner,
                                        const ConstDoublePtr*& xpapadblValue,
                                        const UnsignedInt32*& xpaintNumValues)
// Search the meta-data for all entries of type double array.
{
    xintNumFound = 0;
    xpapachrName = NULL;
    xpaintOwner = NULL;
    xpapadblValue = NULL;
    xpaintNumValues = NULL;
    try {  RebuildSearchBuff(eDoubleArray);  }
    catch(CCubitFile::EErrorCode xeErrorCode)  {  return xeErrorCode;  }
    catch(...)  {  return CCubitFile::eUnknownError;  }

    if((xintNumFound = SearchAll(eDoubleArray)) > 0) {
        xpapachrName = mSearchBuff.mpapachrNames;
        xpaintOwner = mSearchBuff.mpaintOwners;
        xpapadblValue = (ConstDoublePtr*)mSearchBuff.mpaValues;
        xpaintNumValues = mSearchBuff.mpaintNumValues;
    }
    return CCubitFile::eSuccess;
}


//-----------------------------------------------------------------------------

void CMetaData::RebuildSearchBuff(UnsignedInt32 xintSearchType)
// Grow the search return buffer of the class to match the current storage
// buffer size and fit it for the passed return type.
{
    if(!mintBufferSize)  return;

    if(mSearchBuff.mintNumOwners < mintBufferSize) {
        if(mSearchBuff.mpaintOwners)
            delete [] mSearchBuff.mpaintOwners;
        mSearchBuff.mpaintOwners = new UnsignedInt32[mintBufferSize];
        mSearchBuff.mintNumOwners = mintBufferSize;
        if(!mSearchBuff.mpaintOwners) {
            mSearchBuff.mintNumOwners = 0;
            throw CCubitFile::eMemoryError;
        }
    }
    if(mSearchBuff.mintNumNames < mintBufferSize) {
        if(mSearchBuff.mpapachrNames)
            delete [] mSearchBuff.mpapachrNames;
        mSearchBuff.mpapachrNames = new ConstCharPtr[mintBufferSize];
        mSearchBuff.mintNumNames = mintBufferSize;
        if(!mSearchBuff.mpapachrNames) {
            mSearchBuff.mintNumNames = 0;
            throw CCubitFile::eMemoryError;
        }
    }
    if((mSearchBuff.mintNumMetaData < mintBufferSize) ||
        (mSearchBuff.mintMetaDataType != xintSearchType)) {
        if(mSearchBuff.mpaValues) {
            switch(mSearchBuff.mintMetaDataType) {
            case eUnsignedInt32:
                delete [] (UnsignedInt32*)mSearchBuff.mpaValues;  break;
            case eString:
                delete [] (ConstCharPtr*)mSearchBuff.mpaValues;  break;
            case eDouble:
                delete [] (double*)mSearchBuff.mpaValues;  break;
            case eUnsignedInt32Array:
                delete [] (ConstUnsignedInt32Ptr*)mSearchBuff.mpaValues;  break;
            case eDoubleArray:
                delete [] (ConstDoublePtr*)mSearchBuff.mpaValues;  break;
            }
        }
        if(mSearchBuff.mpaintNumValues)
            delete [] mSearchBuff.mpaintNumValues;
        mSearchBuff.mintMetaDataType = xintSearchType;
        mSearchBuff.mintNumMetaData = mintBufferSize;

        switch(mSearchBuff.mintMetaDataType) {
        case eUnsignedInt32:
            mSearchBuff.mpaValues = new UnsignedInt32[mintBufferSize];
            break;
        case eString:
            mSearchBuff.mpaValues = new ConstCharPtr[mintBufferSize];
            break;
        case eDouble:
            mSearchBuff.mpaValues = new double[mintBufferSize];
            break;
        case eUnsignedInt32Array:
            mSearchBuff.mpaValues = new ConstUnsignedInt32Ptr[mintBufferSize];
            break;
        case eDoubleArray:
            mSearchBuff.mpaValues = new ConstDoublePtr[mintBufferSize];
            break;
        }
        if(!mSearchBuff.mpaValues) {
            mSearchBuff.mintNumMetaData = 0;
            throw CCubitFile::eMemoryError;
        }

        switch(mSearchBuff.mintMetaDataType) {
        case eUnsignedInt32Array:
        case eDoubleArray:
            mSearchBuff.mpaintNumValues = new UnsignedInt32[mintBufferSize];
            if(!mSearchBuff.mpaintNumValues) {
                mSearchBuff.mintMetaDataType = 0;
                switch(mSearchBuff.mintMetaDataType) {
                case eUnsignedInt32:
                    delete [] (UnsignedInt32*)mSearchBuff.mpaValues;  break;
                case eString:
                    delete [] (ConstCharPtr*)mSearchBuff.mpaValues;  break;
                case eDouble:
                    delete [] (double*)mSearchBuff.mpaValues;  break;
                case eUnsignedInt32Array:
                    delete [] (ConstUnsignedInt32Ptr*)mSearchBuff.mpaValues;  break;
                case eDoubleArray:
                    delete [] (ConstDoublePtr*)mSearchBuff.mpaValues;  break;
                }
                mSearchBuff.mpaValues = NULL;
                throw CCubitFile::eMemoryError;
            }
            break;
        default:
            mSearchBuff.mpaintNumValues = NULL;
        }
    }
}

UnsignedInt32 CMetaData::SearchByID(UnsignedInt32 xintSearchType,
                                    UnsignedInt32 xintOwner)
// Traverses the meta-data entries to find entries with an owner id matching
// the passed one and copies the values and names of those entries to the
// search buffer.
// Private function, returns the number of entries meeting the search criteria,
// no exceptions thrown.
{
    if(!mHeader.mintMetaDataCount)  return 0;  // none found

    UnsignedInt32 lintNumFound = 0;
    for(UnsignedInt32 lintEntry = 0; lintEntry <
        mHeader.mintMetaDataCount; lintEntry++) {
        if(mpaEntries[lintEntry].mintMetaDataOwner == xintOwner) {
            if(mpaEntries[lintEntry].mintMetaDataType == xintSearchType) {
                mSearchBuff.mpapachrNames[lintNumFound] =
                    mpaEntries[lintEntry].mpachrName;
                switch(xintSearchType) {
                case eUnsignedInt32:
                    ((UnsignedInt32*)mSearchBuff.mpaValues)[lintNumFound] =
                        mpaEntries[lintEntry].muValue.mint;
                    break;
                case eString:
                    ((ConstCharPtr*)mSearchBuff.mpaValues)[lintNumFound] =
                        mpaEntries[lintEntry].muValue.mpachr;
                    break;
                case eDouble:
                    ((double*)mSearchBuff.mpaValues)[lintNumFound] =
                        mpaEntries[lintEntry].muValue.mdbl;
                    break;
                case eUnsignedInt32Array:
                    ((ConstUnsignedInt32Ptr*)mSearchBuff.mpaValues)[lintNumFound] =
                        mpaEntries[lintEntry].muValue.mpaint;
                    mSearchBuff.mpaintNumValues[lintNumFound] =
                        mpaEntries[lintEntry].mintNumValues;
                    break;
                case eDoubleArray:
                    ((ConstDoublePtr*)mSearchBuff.mpaValues)[lintNumFound] =
                        mpaEntries[lintEntry].muValue.mpadbl;
                    mSearchBuff.mpaintNumValues[lintNumFound] =
                        mpaEntries[lintEntry].mintNumValues;
                    break;
                }
                lintNumFound++;
            }
        }
    }
    return lintNumFound;  // return number found
}

UnsignedInt32 CMetaData::SearchByName(UnsignedInt32 xintSearchType,
                                      const char* xpachrName)
// Traverses the meta-data entries to find entries with an name matching
// the passed one and copies the owner id's and names of those entries to the
// search buffer.
// Private function, returns the number of entries meeting the search criteria,
// no exceptions thrown.
{
    if(!mHeader.mintMetaDataCount)  return 0;  // none found

    UnsignedInt32 lintNumFound = 0;
    for(UnsignedInt32 lintEntry = 0; lintEntry <
        mHeader.mintMetaDataCount; lintEntry++) {
        if(!strcmp(xpachrName, mpaEntries[lintEntry].mpachrName)) {
            if(mpaEntries[lintEntry].mintMetaDataType == xintSearchType) {
                mSearchBuff.mpaintOwners[lintNumFound] =
                    mpaEntries[lintEntry].mintMetaDataOwner;
                switch(xintSearchType) {
                case eUnsignedInt32:
                    ((UnsignedInt32*)mSearchBuff.mpaValues)[lintNumFound] =
                        mpaEntries[lintEntry].muValue.mint;
                    break;
                case eString:
                    ((ConstCharPtr*)mSearchBuff.mpaValues)[lintNumFound] =
                        mpaEntries[lintEntry].muValue.mpachr;
                    break;
                case eDouble:
                    ((double*)mSearchBuff.mpaValues)[lintNumFound] =
                        mpaEntries[lintEntry].muValue.mdbl;
                    break;
                case eUnsignedInt32Array:
                    ((ConstUnsignedInt32Ptr*)mSearchBuff.mpaValues)[lintNumFound] =
                        mpaEntries[lintEntry].muValue.mpaint;
                    mSearchBuff.mpaintNumValues[lintNumFound] =
                        mpaEntries[lintEntry].mintNumValues;
                    break;
                case eDoubleArray:
                    ((ConstDoublePtr*)mSearchBuff.mpaValues)[lintNumFound] =
                        mpaEntries[lintEntry].muValue.mpadbl;
                    mSearchBuff.mpaintNumValues[lintNumFound] =
                        mpaEntries[lintEntry].mintNumValues;
                    break;
                }
                lintNumFound++;
            }
        }
    }
    return lintNumFound;  // return number found
}

UnsignedInt32 CMetaData::SearchAll(UnsignedInt32 xintSearchType)
// Traverses the meta-data entries to find all entries of the passed type.
// Private function, returns the number of entries meeting the search criteria,
// no exceptions thrown.
{
    if(!mHeader.mintMetaDataCount)  return 0;  // none found

    UnsignedInt32 lintNumFound = 0;
    for(UnsignedInt32 lintEntry = 0; lintEntry <
        mHeader.mintMetaDataCount; lintEntry++) {
        if(mpaEntries[lintEntry].mintMetaDataType == xintSearchType) {
            mSearchBuff.mpapachrNames[lintNumFound] =
                mpaEntries[lintEntry].mpachrName;
            mSearchBuff.mpaintOwners[lintNumFound] =
                mpaEntries[lintEntry].mintMetaDataOwner;
            switch(xintSearchType) {
            case eUnsignedInt32:
                ((UnsignedInt32*)mSearchBuff.mpaValues)[lintNumFound] =
                    mpaEntries[lintEntry].muValue.mint;
                break;
            case eString:
                ((ConstCharPtr*)mSearchBuff.mpaValues)[lintNumFound] =
                    mpaEntries[lintEntry].muValue.mpachr;
                break;
            case eDouble:
                ((double*)mSearchBuff.mpaValues)[lintNumFound] =
                    mpaEntries[lintEntry].muValue.mdbl;
                break;
            case eUnsignedInt32Array:
                ((ConstUnsignedInt32Ptr*)mSearchBuff.mpaValues)[lintNumFound] =
                    mpaEntries[lintEntry].muValue.mpaint;
                mSearchBuff.mpaintNumValues[lintNumFound] =
                    mpaEntries[lintEntry].mintNumValues;
                break;
            case eDoubleArray:
                ((ConstDoublePtr*)mSearchBuff.mpaValues)[lintNumFound] =
                    mpaEntries[lintEntry].muValue.mpadbl;
                mSearchBuff.mpaintNumValues[lintNumFound] =
                    mpaEntries[lintEntry].mintNumValues;
                break;
            }
            lintNumFound++;
        }
    }
    return lintNumFound;  // return number found
}



///////////////////////////////////////////////////////////////////////////////
// Deletion functions
///////////////////////////////////////////////////////////////////////////////

UnsignedInt32 CMetaData::ClearMetaDataForID(UnsignedInt32 xintOwner)
// Delete all meta-data entries with the same owner.
{
    if(mHeader.mintMetaDataCount) {
        for(UnsignedInt32 lintEntry = 0; lintEntry <
            mHeader.mintMetaDataCount; lintEntry++) {
            if(mpaEntries[lintEntry].mintMetaDataOwner == xintOwner)
                FreeAt(lintEntry);
        }
    }
    return CCubitFile::eSuccess;
}

void CMetaData::FreeValueAt(UnsignedInt32 xintEntry)
// Frees any allocated memory associated with the value of the meta-data stored
// at the passed index.  Private function, no exceptions thrown.
{
    switch(mpaEntries[xintEntry].mintMetaDataType) {
    case eString:
        if(mpaEntries[xintEntry].muValue.mpachr)
            delete [] mpaEntries[xintEntry].muValue.mpachr;
        mpaEntries[xintEntry].muValue.mpachr = NULL;
        break;
    case eUnsignedInt32Array:
        if(mpaEntries[xintEntry].muValue.mpaint)
            delete [] mpaEntries[xintEntry].muValue.mpaint;
        mpaEntries[xintEntry].muValue.mpaint = NULL;
        break;
    case eDoubleArray:
        if(mpaEntries[xintEntry].muValue.mpadbl)
            delete [] mpaEntries[xintEntry].muValue.mpadbl;
        mpaEntries[xintEntry].muValue.mpadbl = NULL;
        break;
    }
}

void CMetaData::FreeAt(UnsignedInt32 xintEntry)
// Frees the meta-data at the passed index and moves all data at higher indices
// down to close the gap.  Private function, no exceptions thrown.
{
    // Free the allocated memory for an entry.
    delete [] mpaEntries[xintEntry].mpachrName;
    FreeValueAt(xintEntry);
    mHeader.mintMetaDataCount--;

    // Move all entries in the meta-data table ahead of the deleted one down
    // to fill in the gap.
    if(xintEntry < mHeader.mintMetaDataCount) {
        for(UnsignedInt32 lintCopyTo = xintEntry; lintCopyTo <
            mHeader.mintMetaDataCount; lintCopyTo++) {
            memcpy(&mpaEntries[lintCopyTo], &mpaEntries[lintCopyTo + 1],
                sizeof(SCubitFileMetaDataEntry));
        }
    }
}

void CMetaData::FreeAll()
// Frees all memory used to store meta-data data and resets the object back
// to its constructed state.  (Has no effect on the return buffers.)
// Private function, no exceptions thrown.
{
    if(mHeader.mintMetaDataCount) {
        for(UnsignedInt32 lintEntry = 0; lintEntry < mintBufferSize; lintEntry++) {
            if(mpaEntries[lintEntry].mpachrName)
                delete [] mpaEntries[lintEntry].mpachrName;
            FreeValueAt(lintEntry);
        }
        delete [] mpaEntries;
        mHeader.mintMetaDataCount = 0;
        mintBufferSize = 0;
        mpaEntries = NULL;
    }
    mHeader.mintMetaDataSchema = 0;
    mHeader.mintMetaDataCompress = 0;
}


///////////////////////////////////////////////////////////////////////////////
// Read / Write File Methods
///////////////////////////////////////////////////////////////////////////////

void CMetaData::WriteMetaData(FILE* xpFile,
                              UnsignedInt32& xintWroteAtOffset,
                              UnsignedInt32& xintLength,
                              UnsignedInt32 xintOffsetFrom)
// Write the meta-data stored in this object to the passed file.  The file
// position, offset from an optional passed value, and length are returned.
{
    if(!xpFile)  throw CCubitFile::eFileWriteError;

    CIOWrapper* lpIO = new CIOWrapper(xpFile);
    xintWroteAtOffset = lpIO->BeginWriteBlock(xintOffsetFrom);

    // Write a description header.
    lpIO->Write(&mHeader.mintMetaDataSchema, 1);
    lpIO->Write(&mHeader.mintMetaDataCompress, 1);
    lpIO->Write(&mHeader.mintMetaDataCount, 1);

    // Write the actual meta-data entries.
    UnsignedInt32 lintEntry, lintNumEntries = mHeader.mintMetaDataCount;
    if(lintNumEntries) {
        for(lintEntry = 0; lintEntry < lintNumEntries; lintEntry++) {
            lpIO->Write(&mpaEntries[lintEntry].mintMetaDataOwner, 1);
            lpIO->Write(&mpaEntries[lintEntry].mintMetaDataType, 1);
            lpIO->Write(mpaEntries[lintEntry].mpachrName);
            switch(mpaEntries[lintEntry].mintMetaDataType) {
            case eUnsignedInt32:
                lpIO->Write(&mpaEntries[lintEntry].muValue.mint, 1);
                break;
            case eString:
                lpIO->Write(mpaEntries[lintEntry].muValue.mpachr);
                break;
            case eDouble:
                lpIO->Write(&mpaEntries[lintEntry].muValue.mdbl, 1);
                break;
            case eUnsignedInt32Array:
                lpIO->Write(&mpaEntries[lintEntry].mintNumValues, 1);
                if(mpaEntries[lintEntry].mintNumValues) {
                    lpIO->Write(mpaEntries[lintEntry].muValue.mpaint,
                        mpaEntries[lintEntry].mintNumValues);
                }
                break;
            case eDoubleArray:
                lpIO->Write(&mpaEntries[lintEntry].mintNumValues, 1);
                if(mpaEntries[lintEntry].mintNumValues) {
                    lpIO->Write(mpaEntries[lintEntry].muValue.mpadbl,
                        mpaEntries[lintEntry].mintNumValues);
                }
                break;
            }
        }
    }
    xintLength = lpIO->EndWriteBlock();
    delete lpIO;
}

void CMetaData::ReadMetaData(FILE* xpFile,
                             UnsignedInt32 xintAbsoluteOffset,
                             UnsignedInt32 xintRelativeOffset,
                             UnsignedInt32 xintSourceEndian)
// Reset the meta-data object and then fill it with data read from a passed
// file at the passed location.  The data will be in the passed endian format
// which is not stored with the meta-data itself as the meta-data is never
// fully independent.
{
    // Delete all stored data
    FreeAll();

    if(!xpFile)  throw CCubitFile::eFileReadError;

    CIOWrapper* lpIO = new CIOWrapper(xpFile, xintSourceEndian);
    lpIO->BeginReadBlock(xintAbsoluteOffset, xintRelativeOffset);

    // Read the description header.
    lpIO->Read(&mHeader.mintMetaDataSchema, 1);
    lpIO->Read(&mHeader.mintMetaDataCompress, 1);
    lpIO->Read(&mHeader.mintMetaDataCount, 1);

    UnsignedInt32 lintEntry, lintNumEntries = mHeader.mintMetaDataCount;
    if(lintNumEntries) {
        // Grow the storage capacity of the meta-data object to match that of
        // the stored data.
        UnsignedInt32 lintGrowSizeBy = mintGrowSizeBy;
        mintGrowSizeBy = lintNumEntries;
        GrowSize();
        mintGrowSizeBy = lintGrowSizeBy;

        // Read the stored meta-data entries.
        for(lintEntry = 0; lintEntry < lintNumEntries; lintEntry++) {
            lpIO->Read(&mpaEntries[lintEntry].mintMetaDataOwner, 1);
            lpIO->Read(&mpaEntries[lintEntry].mintMetaDataType, 1);
            mpaEntries[lintEntry].mpachrName = lpIO->Read();
            switch(mpaEntries[lintEntry].mintMetaDataType) {
            case eUnsignedInt32:
                lpIO->Read(&mpaEntries[lintEntry].muValue.mint, 1);
                break;
            case eString:
                mpaEntries[lintEntry].muValue.mpachr = lpIO->Read();
                break;
            case eDouble:
                lpIO->Read(&mpaEntries[lintEntry].muValue.mdbl, 1);
                break;
            case eUnsignedInt32Array:
                lpIO->Read(&mpaEntries[lintEntry].mintNumValues, 1);
                if(mpaEntries[lintEntry].mintNumValues) {
                    mpaEntries[lintEntry].muValue.mpaint =
                        new UnsignedInt32[mpaEntries[lintEntry].mintNumValues];
                    if(!mpaEntries[lintEntry].muValue.mpaint)
                        throw CCubitFile::eMemoryError;
                    lpIO->Read(mpaEntries[lintEntry].muValue.mpaint,
                        mpaEntries[lintEntry].mintNumValues);
                }
                break;
            case eDoubleArray:
                lpIO->Read(&mpaEntries[lintEntry].mintNumValues, 1);
                if(mpaEntries[lintEntry].mintNumValues) {
                    mpaEntries[lintEntry].muValue.mpadbl =
                        new double[mpaEntries[lintEntry].mintNumValues];
                    if(!mpaEntries[lintEntry].muValue.mpadbl)
                        throw CCubitFile::eMemoryError;
                    lpIO->Read(mpaEntries[lintEntry].muValue.mpadbl,
                        mpaEntries[lintEntry].mintNumValues);
                }
                break;
            }
        }
    }
    delete lpIO;
}

