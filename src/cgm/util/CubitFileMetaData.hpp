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

 
 Filename      : CubitFileMetaData.hpp

 Purpose       : Declares a meta-data (descriptions) management interface for
                 the Cubit (*.cub) format file.
           
 Special Notes :

 Creator       : Will A. Helden

 Creation Date : 02/15/02

 Owner         : Will A. Helden

*******************************************************************************/

#ifndef CubitFileMetaData_HPP
#define CubitFileMetaData_HPP

#include "CubitFile.hpp"
#include "CubitUtilConfigure.h"

namespace NCubitFile {

typedef const UnsignedInt32* ConstUnsignedInt32Ptr;
typedef const double* ConstDoublePtr;

class CUBIT_UTIL_EXPORT CMetaData {
public:
    CMetaData(UnsignedInt32 xintGrowSizeBy = 64);
    virtual ~CMetaData();
    
    UnsignedInt32 GetValue(UnsignedInt32 xintOwner,
        ConstCharPtr xpachrName, UnsignedInt32& xintValue);
    UnsignedInt32 GetValue(UnsignedInt32 xintOwner,
        ConstCharPtr xpachrName, ConstCharPtr& xpachrValue);
    UnsignedInt32 GetValue(UnsignedInt32 xintOwner,
        ConstCharPtr xpachrName, double& xdblValue);
    UnsignedInt32 GetValue(UnsignedInt32 xintOwner,
        ConstCharPtr xpachrName, const UnsignedInt32*& xpaintValue,
        UnsignedInt32& xintNumValues);
    UnsignedInt32 GetValue(UnsignedInt32 xintOwner,
        ConstCharPtr xpachrName, const double*& xpadblValue,
        UnsignedInt32& xintNumValues);

    UnsignedInt32 SetValue(UnsignedInt32 xintOwner,
        ConstCharPtr xpachrName, UnsignedInt32 xintValue);
    UnsignedInt32 SetValue(UnsignedInt32 xintOwner,
        ConstCharPtr xpachrName, ConstCharPtr xpachrValue);
    UnsignedInt32 SetValue(UnsignedInt32 xintOwner,
        ConstCharPtr xpachrName, double xdblValue);
    UnsignedInt32 SetValue(UnsignedInt32 xintOwner,
        ConstCharPtr xpachrName, const UnsignedInt32* xpaintValue,
        UnsignedInt32 xintNumValues);
    UnsignedInt32 SetValue(UnsignedInt32 xintOwner,
        ConstCharPtr xpachrName, const double* xpadblValue,
        UnsignedInt32 xintNumValues);
    
#ifdef BOYD15
    UnsignedInt32 GetMetaDataByID(
        UnsignedInt32 xintOwner, UnsignedInt32& xintNumFound,
        const ConstCharPtr*& xpapachrName, const UnsignedInt32*& xpaintValue);
    UnsignedInt32 GetMetaDataByID(
        UnsignedInt32 xintOwner, UnsignedInt32& xintNumFound,
        const ConstCharPtr*& xpapachrName, const ConstCharPtr*& xpapachrValue);
    UnsignedInt32 GetMetaDataByID(
        UnsignedInt32 xintOwner, UnsignedInt32& xintNumFound,
        const ConstCharPtr*& xpapachrName, const double*& xpadblValue);
    UnsignedInt32 GetMetaDataByID(
        UnsignedInt32 xintOwner, UnsignedInt32& xintNumFound,
        const ConstCharPtr*& xpapachrName, const ConstUnsignedInt32Ptr*& xpapaintValue,
        const UnsignedInt32*& xpaintNumValues);
    UnsignedInt32 GetMetaDataByID(
        UnsignedInt32 xintOwner, UnsignedInt32& xintNumFound,
        const ConstCharPtr*& xpapachrName, const ConstDoublePtr*& xpapadblValue,
        const UnsignedInt32*& xpaintNumValues);

    UnsignedInt32 GetMetaDataByName(
        const char* xpachrName, UnsignedInt32& xintNumFound,
        const UnsignedInt32*& xpaintOwner, const UnsignedInt32*& xpaintValue);
    UnsignedInt32 GetMetaDataByName(
        const char* xpachrName, UnsignedInt32& xintNumFound,
        const UnsignedInt32*& xpaintOwner, const ConstCharPtr*& xpapachrValue);
    UnsignedInt32 GetMetaDataByName(
        const char* xpachrName, UnsignedInt32& xintNumFound,
        const UnsignedInt32*& xpaintOwner, const double*& xpadblValue);
    UnsignedInt32 GetMetaDataByName(
        const char* xpachrName, UnsignedInt32& xintNumFound,
        const UnsignedInt32*& xpaintOwner, const ConstUnsignedInt32Ptr*& xpapaintValue,
        const UnsignedInt32*& xpaintNumValues);
    UnsignedInt32 GetMetaDataByName(
        const char* xpachrName, UnsignedInt32& xintNumFound,
        const UnsignedInt32*& xpaintOwner, const ConstDoublePtr*& xpapadblValue,
        const UnsignedInt32*& xpaintNumValues);

    UnsignedInt32 ClearValue(UnsignedInt32 xintOwner, ConstCharPtr xpachrName);
#endif
    UnsignedInt32 ClearMetaDataForID(UnsignedInt32 xintOwner);

    UnsignedInt32 GetMetaDataAll(UnsignedInt32& xintNumFound,
        const ConstCharPtr*& xpapachrName, const UnsignedInt32*& xpaintOwner,
        const UnsignedInt32*& xpaintValue);
    UnsignedInt32 GetMetaDataAll(UnsignedInt32& xintNumFound,
        const ConstCharPtr*& xpapachrName, const UnsignedInt32*& xpaintOwner,
        const ConstCharPtr*& xpapachrValue);
    UnsignedInt32 GetMetaDataAll(UnsignedInt32& xintNumFound,
        const ConstCharPtr*& xpapachrName, const UnsignedInt32*& xpaintOwner,
        const double*& xpadblValue);
    UnsignedInt32 GetMetaDataAll(UnsignedInt32& xintNumFound,
        const ConstCharPtr*& xpapachrName, const UnsignedInt32*& xpaintOwner,
        const ConstUnsignedInt32Ptr*& xpapaintValue,
        const UnsignedInt32*& xpaintNumValues);
    UnsignedInt32 GetMetaDataAll(UnsignedInt32& xintNumFound,
        const ConstCharPtr*& xpapachrName, const UnsignedInt32*& xpaintOwner,
        const ConstDoublePtr*& xpapadblValue,
        const UnsignedInt32*& xpaintNumValues);
    
    void WriteMetaData(FILE* xpFile,
        UnsignedInt32& xintWroteAtOffset,
        UnsignedInt32& xintLength,
        UnsignedInt32 xintOffsetFrom = 0);
    void ReadMetaData(FILE* xpFile,
        UnsignedInt32 xintAbsoluteOffset,
        UnsignedInt32 xintRelativeOffset,
        UnsignedInt32 xintSourceEndian);

private:
    enum EMetaDataType { eUnsignedInt32, eString, eDouble,
        eUnsignedInt32Array, eDoubleArray };

    void GrowSize();
    UnsignedInt32 Find(UnsignedInt32 xintOwner,
        ConstCharPtr xpachrName, UnsignedInt32& xintIndex);
    void FreeAll();
    void FreeAt(UnsignedInt32 xintEntry);
    void FreeValueAt(UnsignedInt32 xintEntry);
    UnsignedInt32 AddEntry(UnsignedInt32 xintOwner, ConstCharPtr xpachrName);
    void RebuildSearchBuff(UnsignedInt32 xintSearchType);
    UnsignedInt32 SearchByID(UnsignedInt32 xintSearchType,
        UnsignedInt32 xintOwner);
    UnsignedInt32 SearchByName(UnsignedInt32 xintSearchType,
        const char* xpachrName);
    UnsignedInt32 SearchAll(UnsignedInt32 xintSearchType);
    
    struct SCubitFileMetaDataHeader {
        UnsignedInt32 mintMetaDataSchema;
        UnsignedInt32 mintMetaDataCompress;
        UnsignedInt32 mintMetaDataCount;
    } mHeader;
    struct SCubitFileMetaDataEntry {
        UnsignedInt32 mintMetaDataOwner;
        UnsignedInt32 mintMetaDataType;
        char* mpachrName;
        union {
            UnsignedInt32 mint;
            char* mpachr;
            double mdbl;
            UnsignedInt32* mpaint;
            double* mpadbl;
        } muValue;
        UnsignedInt32 mintNumValues;
    } *mpaEntries;

    UnsignedInt32 mintBufferSize;
    UnsignedInt32 mintGrowSizeBy;

    struct {
        UnsignedInt32 mintNumOwners;
        UnsignedInt32* mpaintOwners;
        UnsignedInt32 mintNumNames;
        ConstCharPtr* mpapachrNames;
        UnsignedInt32 mintMetaDataType;
        UnsignedInt32 mintNumMetaData;
        void* mpaValues;
        UnsignedInt32* mpaintNumValues;
    } mSearchBuff;
};

} // namespace NCubitFile

#endif 

