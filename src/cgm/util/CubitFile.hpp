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

 
 Filename      : CubitFile.hpp

 Purpose       : Declares an interface for reading and writing data in the
                 Cubit (*.cub) file format.
           
 Special Notes :

 Creator       : Will A. Helden

 Creation Date : 02/15/02

 Owner         : Will A. Helden

*******************************************************************************/

#ifndef CubitFile_HPP
#define CubitFile_HPP

#include <stdio.h>
#include <memory.h>
#include "CubitUtilConfigure.h"

namespace NCubitFile {

class CMetaData;
class CFEModel;

typedef unsigned int UnsignedInt32;  // This must always be a 32 bit integer!!!
typedef const char* ConstCharPtr;

struct SModelData {
    UnsignedInt32 mintModelHandle;
    UnsignedInt32 mintModelType;
    UnsignedInt32 mintModelOwner;
};
struct SElemData {
    UnsignedInt32 mintElemType;
    UnsignedInt32 mintElemOrder;
    UnsignedInt32 mintElemCount;
    UnsignedInt32* mpaElemIDs;
    UnsignedInt32* mpaElemConnect;
};
struct SGroupData {
    UnsignedInt32 mintMemberType;
    UnsignedInt32 mintMemberCount;
    UnsignedInt32* mpaMemberIDs;
};
typedef SGroupData SBlockData;
typedef SGroupData SNodeSetData;

// first edition of sideset data
struct SSideSetData_10 {
    UnsignedInt32 mintMemberType;
    UnsignedInt32 mintMemberCount;
    UnsignedInt32 mintMemberSenseSize;
    UnsignedInt32* mpaintMemberIDs;
    char* mpachrMemberSense;
    UnsignedInt32* mpaintMemberSense;
    char* mpachrMemberSideNumber;
};
// next edition of sideset data holds more information in a more flexible way
struct SSideSetData_11 {
    UnsignedInt32 mintMemberCount;  // number of members {type,id} ...
    UnsignedInt32* mpaintMemberTypes; // types of each member
    UnsignedInt32* mpaintMemberIDs; // ids of the members
    char* mpachrMemberSenses;        // shell sense (forward, reverse, both)
    UnsignedInt32* mpaintMemberWRTEntities; // wrt entities, format is  ( {num_wrt} {type,id} {type,id} ...
                                            //                            {num_wrt} {type,id} ...
};

class CUBIT_UTIL_EXPORT CCubitFile  
{
public:
    typedef UnsignedInt32 HModel;
    enum EErrorCode { eSuccess = 0,
        eFileReadError, eFileSeekError, eFileTellError, eFileWriteError,
        eFileUnrecognizedFormat,
        eDuplicateWrite, eCorruptBlock,
        eNotFound, eOrderError,
        eMemoryError, ePassedNullPointer,
        eUnknownError
    };
    enum EMetaDataOwnerType {
        eModelMetaData,
        eGeomMetaData, eNodeMetaData, eElemMetaData, eGroupMetaData,
        eBlockMetaData, eNodeSetMetaData, eSideSetMetaData
    };
    enum EModelType {
      eFEModel, eACISText, eACISBinary, eAssemblyModel
    };
    enum ESideSetSenseSize {
        eSideSetSenseNone, eSideSetSenseByte, eSideSetSenseInt32
    };

    CCubitFile();
    virtual ~CCubitFile();

    EErrorCode Open(const char* xstrReadFileName, const char* xstrWriteFileName,
        const char* xstrBackupFileName = NULL);
    EErrorCode Close();

    UnsignedInt32 GetModelList(UnsignedInt32& xintNumModels,
        const SModelData*& xpaModels);
    UnsignedInt32 IsModel(HModel xintModel);
    UnsignedInt32 GetReadModelLength(HModel xintModel,
                                     UnsignedInt32& xintLength);
    UnsignedInt32 CreateModel(EModelType xeModelType, HModel& xintModel);
    UnsignedInt32 DeleteModel(HModel xintModel);
    UnsignedInt32 GetModelOwner(HModel xintModel, HModel& xintOwner);
    UnsignedInt32 SetModelOwner(HModel xintModel, HModel xintOwner);

    // Returns a FILE* that can be directly written to.  Be sure to
    // call EndWriteAssemblyModel() to correctly accomodate the written data.
  UnsignedInt32 BeginWriteAssemblyModel(HModel xintModelId,
                                        FILE*& writeable_file);
  UnsignedInt32 EndWriteAssemblyModel(HModel xintModelId);
    // Returns a FILE* that can be read from.  The FILE* will NOT return EOF
    // at the end of the data...the reader should not go GetReadModelLength()
    // bytes beyond the current FILE* position.  You should copy the contents
    // to another FILE if this is an issue for the reader.
  UnsignedInt32 BeginReadAssemblyModel(HModel xintGeomModel,
                                       FILE*& readable_file);
  UnsignedInt32 EndReadAssemblyModel();
  
  UnsignedInt32 BeginWriteGeomModel(HModel xintGeomModel,
                                    const char* xstrGeomFile);
  UnsignedInt32 EndWriteGeomModel();
  UnsignedInt32 BeginReadGeomModel(HModel xintGeomModel,
                                   FILE* xpGeomFile);
  UnsignedInt32 EndReadGeomModel();
  
  UnsignedInt32 BeginWriteFEModel(HModel xintFEModel,
                                  UnsignedInt32 xintGeomCount,
                                  UnsignedInt32 xintGroupCount,
                                  UnsignedInt32 xintBlockCount,
                                  UnsignedInt32 xintNodeSetCount,
                                  UnsignedInt32 xintSideSetCount);
  UnsignedInt32 WriteNodes(UnsignedInt32 xintIndex, UnsignedInt32 xintGeomID,
                           UnsignedInt32 xintNodeCount, UnsignedInt32 *xpaintNodeIDs,
                           double *xpadblX, double *xpadblY, double *xpadblZ);
  UnsignedInt32 WriteElems(UnsignedInt32 xintIndex,
                           UnsignedInt32 xintNumTypes, SElemData* xpaElemData);
  UnsignedInt32 WriteGroup(UnsignedInt32 xintIndex,
                           UnsignedInt32 xintGroupID, UnsignedInt32 xintGroupType,
                           ConstCharPtr xpachrGroupName,
                           UnsignedInt32 xintNumTypes, SGroupData* xpaGroupData);
  UnsignedInt32 WriteBlock(UnsignedInt32 xintIndex,
                           UnsignedInt32 xintBlockID, UnsignedInt32 xintBlockType,
                           UnsignedInt32 xintBlockColor, UnsignedInt32 xintMixedElemType,
                           UnsignedInt32 xintDefPyramidType, UnsignedInt32 xintMaterialID,
                           UnsignedInt32 xintBlockDimension,
                           UnsignedInt32 xintNumTypes, SBlockData* xpaBlockData,
                           UnsignedInt32 xintAttributeOrder, double* xpadblAttributes);
  UnsignedInt32 WriteNodeSet(UnsignedInt32 xintIndex,
                             UnsignedInt32 xintNodeSetID, UnsignedInt32 xintColor,
                             UnsignedInt32 xintPointSymbol,
                             UnsignedInt32 xintNumTypes, SNodeSetData* xpaNodeSetData);
  UnsignedInt32 WriteSideSet_11(UnsignedInt32 xintIndex,
                             UnsignedInt32 xintSideSetID, UnsignedInt32 xintColor,
                             UnsignedInt32 xintUseShells,
                             UnsignedInt32 xintNumTypes, SSideSetData_11* xpaSideSetData,
                             UnsignedInt32 xintNumDistFact, double* xpadblDistribution);
  UnsignedInt32 EndWriteFEModel();
  
  UnsignedInt32 BeginReadFEModel(HModel xintFEModel,
                                 UnsignedInt32& xintGeomCount,
                                 UnsignedInt32& xintGroupCount,
                                 UnsignedInt32& xintBlockCount,
                                 UnsignedInt32& xintNodeSetCount,
                                 UnsignedInt32& xintSideSetCount);
  UnsignedInt32 ReadNodes(UnsignedInt32 xintIndex, UnsignedInt32& xintGeomID,
                          UnsignedInt32& xintNodeCount, UnsignedInt32*& xpaintNodeIDs,
                          double*& xpadblX, double*& xpadblY, double*& xpadblZ);
  UnsignedInt32 ReadElems(UnsignedInt32 xintIndex, UnsignedInt32& xintGeomID,
                          UnsignedInt32& xintNumTypes, SElemData*& xpaElemData);
  UnsignedInt32 ReadGroupIdentity(UnsignedInt32 xintIndex,
                                  UnsignedInt32& xintGroupID,
                                  UnsignedInt32& xintGroupType,
                                  ConstCharPtr& xpachrGroupName);
  UnsignedInt32 ReadGroupMembers(UnsignedInt32 xintIndex,
                                 UnsignedInt32& xintNumTypes, SGroupData*& xpaGroupData);
    UnsignedInt32 ReadBlock(UnsignedInt32 xintIndex,
        UnsignedInt32& xintBlockID, UnsignedInt32& xintBlockType,
        UnsignedInt32& xintBlockColor, UnsignedInt32& xintMixedElemType,
        UnsignedInt32& xintDefPyramidType, UnsignedInt32& xintMaterialID,
        UnsignedInt32& xintBlockDimension,
        UnsignedInt32& xintNumTypes, SBlockData*& xpaBlockData,
        UnsignedInt32& xintAttributeOrder, double*& xpadblAttributes);
    UnsignedInt32 ReadNodeSet(UnsignedInt32 xintIndex,
        UnsignedInt32& xintNodeSetID, UnsignedInt32& xintColor,
        UnsignedInt32& xintPointSymbol,
        UnsignedInt32& xintNumTypes, SNodeSetData*& xpaNodeSetData);

    // read old sideset format
    UnsignedInt32 ReadSideSet_10(UnsignedInt32 xintIndex,
        UnsignedInt32& xintSideSetID, UnsignedInt32& xintColor,
        UnsignedInt32& xintUseShells,
        UnsignedInt32& xintNumTypes, SSideSetData_10*& xpaSideSetData,
        UnsignedInt32& xintNumDistFact, double*& xpadblDistribution);
    
    // read new sideset format
    UnsignedInt32 ReadSideSet_11(UnsignedInt32 xintIndex,
        UnsignedInt32& xintSideSetID, UnsignedInt32& xintColor,
        UnsignedInt32& xintUseShells,
        UnsignedInt32& xintNumTypes, SSideSetData_11*& xpaSideSetData,
        UnsignedInt32& xintNumDistFact, double*& xpadblDistribution);
    
    UnsignedInt32 EndReadFEModel();

    UnsignedInt32 GetReadMetaData(EMetaDataOwnerType xeType,
        CMetaData*& xpMetaData);
    UnsignedInt32 GetWriteMetaData(EMetaDataOwnerType xeType,
        CMetaData*& xpMetaData);

EErrorCode GetError() const; // ???

//  enum EEndian { eLittleEndian = 0, eBigEndian = 0xFFFFFFFF };
    static const UnsignedInt32 mintNativeEndian;
    enum ECompression { eNoCompression };

private:
    char* mstrReadFileName;
    char* mstrWriteFileName;
    char* mstrBackupFileName;
    UnsignedInt32 mintWriteTempFile;
    FILE* mpReadFile;
    FILE* mpWriteFile;
EErrorCode meErrorState;
    UnsignedInt32 mintWriteBuffNumModels;
    UnsignedInt32 mintNextModelID;

    struct SCubitFileContentsHeader {
        UnsignedInt32 mintHeaderSourceEndian;
        UnsignedInt32 mintHeaderSchema;
        UnsignedInt32 mintNumModels;
        UnsignedInt32 mintModelTableOffset;
        UnsignedInt32 mintModelMetaDataOffset;
        UnsignedInt32 mintActiveFEModel;
    } mReadContents, mWriteContents;
    static const UnsignedInt32 mintSizeOfContents;
    struct SCubitFileModelEntry {
        UnsignedInt32 mintModelHandle;
        UnsignedInt32 mintModelOffset;
        UnsignedInt32 mintModelLength;
        UnsignedInt32 mintModelType;
        UnsignedInt32 mintModelOwner;
        UnsignedInt32 mintModel64bitOSPad;  // make struct 64 bit word aligned
    } *mpaReadModels, *mpaWriteModels;
    static const UnsignedInt32 mintSizeOfModel;
  enum EModelStat { eStatNotWritten, eStatWritten, eStatDelete, eStatWriting };
    EModelStat* mpaReadModelStat;

    void WriteModelTable();
    UnsignedInt32 FindModel(HModel xintModel, SCubitFileModelEntry* xpaModels,
        UnsignedInt32 xintNumModels, UnsignedInt32& xintIndex);
    void CopyModel(UnsignedInt32 xintReadOffset, UnsignedInt32& xintWriteOffset,
        UnsignedInt32 xintLength, FILE* xpReadFile, FILE* xpWriteFile);
    void FreeAll();

    UnsignedInt32 mintFEModelIndex;
    CFEModel* mpReadFEModel;
    CFEModel* mpWriteFEModel;
    CMetaData* mpMetaData;
    struct {
        UnsignedInt32 mintNumModels;
        SModelData* mpaModelData;
    } mModelBuff;
};

} // namespace NCubitFile

#endif

