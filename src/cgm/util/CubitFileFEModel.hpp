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

 
 Filename      : CubitFileFEModel.hpp

 Purpose       : Defines an interface for the reading and writing functionality
                 for a FE model section of a Cubit (*.cub) format file.
           
 Special Notes :

 Creator       : Will A. Helden

 Creation Date : 02/15/02

 Owner         : Will A. Helden

*******************************************************************************/

#include "CubitFileMetaData.hpp"
#include "CubitUtilConfigure.h"

#ifndef CubitFileFEModel_HPP
#define CubitFileFEModel_HPP

namespace NCubitFile {

class CUBIT_UTIL_EXPORT CFEModel {
public:
    CFEModel();
    virtual ~CFEModel();
    
    UnsignedInt32 InitWrite(FILE* xpFile,
        UnsignedInt32 xintGeomCount, UnsignedInt32 xintGroupCount,
        UnsignedInt32 xintBlockCount, UnsignedInt32 xintNodeSetCount,
        UnsignedInt32 xintSideSetCount);
    void WriteNodes(UnsignedInt32 xintIndex, UnsignedInt32 xintGeomID,
        UnsignedInt32 xintNodeCount, UnsignedInt32 *xpaintNodeIDs,
        double *xpadblX, double *xpadblY, double *xpadblZ);
    void WriteElems(UnsignedInt32 xintIndex,
        UnsignedInt32 xintNumTypes, SElemData* xpaElemData);
    void WriteGroup(UnsignedInt32 xintIndex, UnsignedInt32 xintGroupID,
        UnsignedInt32 xintGroupType, const char* xpachrGroupName,
        UnsignedInt32 xintNumTypes, SGroupData* xpaGroupData);
    void WriteBlock(UnsignedInt32 xintIndex,
        UnsignedInt32 xintBlockID, UnsignedInt32 xintBlockType,
        UnsignedInt32 xintBlockColor, UnsignedInt32 xintMixedElemType,
        UnsignedInt32 xintDefPyramidType, UnsignedInt32 xintMaterialID,
        UnsignedInt32 xintBlockDimension,
        UnsignedInt32 xintNumTypes, SBlockData* xpaBlockData,
        UnsignedInt32 xintAttributeOrder, double* xpadblAttributes);
    void WriteNodeSet(UnsignedInt32 xintIndex,
        UnsignedInt32 xintNodeSetID, UnsignedInt32 xintColor,
        UnsignedInt32 xintPointSymbol,
        UnsignedInt32 xintNumTypes, SNodeSetData* xpaNodeSetData);
    void WriteSideSet_11(UnsignedInt32 xintIndex,
        UnsignedInt32 xintSideSetID, UnsignedInt32 xintColor,
        UnsignedInt32 xintUseShells,
        UnsignedInt32 xintNumTypes, SSideSetData_11* xpaSideSetData,
        UnsignedInt32 xintNumDistFact, double* xpadblDistribution);
    UnsignedInt32 EndWrite();
    
    void InitRead(FILE* xpFile, UnsignedInt32 xintAbsoluteOffset,
        UnsignedInt32& xintGeomCount, UnsignedInt32& xintGroupCount,
        UnsignedInt32& xintBlockCount, UnsignedInt32& xintNodeSetCount,
        UnsignedInt32& xintSideSetCount); 
    void ReadNodes(UnsignedInt32 xintIndex, UnsignedInt32& xintGeomID,
        UnsignedInt32& xintNodeCount, UnsignedInt32*& xpaintNodeIDs,
        double*& xpadblX, double*& xpadblY, double*& xpadblZ);
    void ReadElems(UnsignedInt32 xintIndex, UnsignedInt32& xintGeomID,
        UnsignedInt32& xintNumTypes, SElemData*& xpaElemData);
    void ReadGroupIdentity(UnsignedInt32 xintIndex, UnsignedInt32& xintGroupID,
        UnsignedInt32& xintGroupType, const char*& xpachrGroupName);
    void ReadGroupMembers(UnsignedInt32 xintIndex,
        UnsignedInt32& xintNumTypes, SGroupData*& xpaGroupData);
    void ReadBlock(UnsignedInt32 xintIndex,
        UnsignedInt32& xintBlockID, UnsignedInt32& xintBlockType,
        UnsignedInt32& xintBlockColor, UnsignedInt32& xintMixedElemType,
        UnsignedInt32& xintDefPyramidType, UnsignedInt32& xintMaterialID,
        UnsignedInt32& xintBlockDimension,
        UnsignedInt32& xintNumTypes, SBlockData*& xpaBlockData,
        UnsignedInt32& xintAttributeOrder, double*& xpadblAttributes);
    void ReadNodeSet(UnsignedInt32 xintIndex,
        UnsignedInt32& xintNodeSetID, UnsignedInt32& xintColor,
        UnsignedInt32& xintPointSymbol,
        UnsignedInt32& xintNumTypes, SNodeSetData*& xpaNodeSetData);
    void ReadSideSet_10(UnsignedInt32 xintIndex,
        UnsignedInt32& xintSideSetID, UnsignedInt32& xintColor,
        UnsignedInt32& xintUseShells,
        UnsignedInt32& xintNumTypes, SSideSetData_10*& xpaSideSetData,
        UnsignedInt32& xintNumDistFact, double*& xpadblDistribution);
    void ReadSideSet_11(UnsignedInt32 xintIndex,
        UnsignedInt32& xintSideSetID, UnsignedInt32& xintColor,
        UnsignedInt32& xintUseShells,
        UnsignedInt32& xintNumTypes, SSideSetData_11*& xpaSideSetData,
        UnsignedInt32& xintNumDistFact, double*& xpadblDistribution);
    void EndRead();
    
    CMetaData& GetGeomMetaData();
    CMetaData& GetNodeMetaData();
    CMetaData& GetElemMetaData();
    CMetaData& GetGroupMetaData();
    CMetaData& GetBlockMetaData();
    CMetaData& GetNodeSetMetaData();
    CMetaData& GetSideSetMetaData();
    
private:
    FILE* mpReadFile;
    FILE* mpWriteFile;
    UnsignedInt32 mintFEModelOffset;
    CMetaData mGeomMetaData;
    CMetaData mNodeMetaData;
    CMetaData mElemMetaData;
    CMetaData mGroupMetaData;
    CMetaData mBlockMetaData;
    CMetaData mNodeSetMetaData;
    CMetaData mSideSetMetaData;

    // Data storage structures:
    // CAUTION: These structures must be 64 bit word aligned!!!
    struct SCubitFileFEModelHeader {
        UnsignedInt32 mintFEModelEndian;
        UnsignedInt32 mintFEModelSchema;
        UnsignedInt32 mintFEModelCompress;
        UnsignedInt32 mintFEModelLength;
        UnsignedInt32 mintGeometryCount;
        UnsignedInt32 mintGeomTableOffset;
        UnsignedInt32 mintGeomMetaDataOffset;
        UnsignedInt32 mintNodeMetaDataOffset;
        UnsignedInt32 mintElemMetaDataOffset;
        UnsignedInt32 mintGroupCount;
        UnsignedInt32 mintGroupTableOffset;
        UnsignedInt32 mintGroupMetaDataOffset;
        UnsignedInt32 mintBlockCount;
        UnsignedInt32 mintBlockTableOffset;
        UnsignedInt32 mintBlockMetaDataOffset;
        UnsignedInt32 mintNodeSetCount;
        UnsignedInt32 mintNodeSetTableOffset;
        UnsignedInt32 mintNodeSetMetaDataOffset;
        UnsignedInt32 mintSideSetCount;
        UnsignedInt32 mintSideSetTableOffset;
        UnsignedInt32 mintSideSetMetaDataOffset;
        UnsignedInt32 mintPadFor64bitOS;
    } mFEModel;
    static const UnsignedInt32 mintSizeOfFEModelHeader;
    struct SCubitFileGeomEntry {
        UnsignedInt32 mintNodeCount;
        UnsignedInt32 mintNodeOffset;
        UnsignedInt32 mintElemCount;
        UnsignedInt32 mintElemOffset;
        UnsignedInt32 mintElemTypeCount;
        UnsignedInt32 mintElemLength;
		UnsignedInt32 mintGeomID;
        UnsignedInt32 mintPadFor64bitOS;
    } *mpaGeoms;
    static const UnsignedInt32 mintSizeOfGeomEntry;
    struct SCubitFileGroupEntry {
        UnsignedInt32 mintGroupID;
        UnsignedInt32 mintGroupType;
        UnsignedInt32 mintMemberCount;
        UnsignedInt32 mintMemberOffset;
        UnsignedInt32 mintMemberTypeCount;
        UnsignedInt32 mintGroupLength;
    } *mpaGroups;
    static const UnsignedInt32 mintSizeOfGroupEntry;
    struct SCubitFileBlockEntry {
        UnsignedInt32 mintBlockID;
        UnsignedInt32 mintBlockElementType;
        UnsignedInt32 mintMemberCount;
        UnsignedInt32 mintMemberOffset;
        UnsignedInt32 mintMemberTypeCount;
        UnsignedInt32 mintAttributeOrder;
        UnsignedInt32 mintBlockColor;
        UnsignedInt32 mintBlockMixedElemType;
        UnsignedInt32 mintBlockDefPyramidType;
        UnsignedInt32 mintBlockMaterial;
        UnsignedInt32 mintBlockLength;
        UnsignedInt32 mintBlockDimension;
    } *mpaBlocks;
    static const UnsignedInt32 mintSizeOfBlockEntry;
    struct SCubitFileNodeSetEntry {
        UnsignedInt32 mintNodeSetID;
        UnsignedInt32 mintMemberCount;
        UnsignedInt32 mintMemberOffset;
        UnsignedInt32 mintMemberTypeCount;
        UnsignedInt32 mintNodeSetPointSym;
        UnsignedInt32 mintNodeSetColor;
        UnsignedInt32 mintNodeSetLength;
        UnsignedInt32 mintPadFor64bitOS;
    } *mpaNodeSets;
    static const UnsignedInt32 mintSizeOfNodeSetEntry;
    struct SCubitFileSideSetEntry {
        UnsignedInt32 mintSideSetID;
        UnsignedInt32 mintMemberCount;
        UnsignedInt32 mintMemberOffset;
        UnsignedInt32 mintMemberTypeCount;
        UnsignedInt32 mintNumDistFact;
        UnsignedInt32 mintSideSetColor;
        UnsignedInt32 mintUseShells;
        UnsignedInt32 mintSideSetLength;
    } *mpaSideSets;
    static const UnsignedInt32 mintSizeOfSideSetEntry;
    
    // Buffers for read operations.
    struct SNodeReturnBuffer {
        UnsignedInt32 mintNumNodes;
        UnsignedInt32* mpaNodeIDs;
        double* mpadblX;
        double* mpadblY;
        double* mpadblZ;
    } mNodeBuff;
    struct SElemReturnBuffer {
        UnsignedInt32 mintNumTypes;
        UnsignedInt32 mintNumElems;
        UnsignedInt32 mintNumConnect;
        SElemData* mpaElemData;
        UnsignedInt32* mpaElemIDs;
        UnsignedInt32* mpaElemConnect;
    } mElemBuff;
    struct SGroupReturnBuffer {
        UnsignedInt32 mintNumTypes;
        UnsignedInt32 mintNumMembers;
        SGroupData* mpaGroupData;
        UnsignedInt32* mpaintMemberIDs;
    } mGroupBuff;
    struct SBlockReturnBuffer {
        UnsignedInt32 mintNumTypes;
        UnsignedInt32 mintNumMembers;
        UnsignedInt32 mintAttributeOrder;
        SBlockData* mpaBlockData;
        UnsignedInt32* mpaintMemberIDs;
        double* mpadblAttributes;
    } mBlockBuff;
    struct SNodeSetReturnBuffer {
        UnsignedInt32 mintNumTypes;
        UnsignedInt32 mintNumMembers;
        SNodeSetData* mpaNodeSetData;
        UnsignedInt32* mpaintMemberIDs;
    } mNodeSetBuff;

    // these sideset buffers are big arrays and the mpaSideSetData can reference parts of these buffers
    // I guess the reason for that is memory ownership is a bit more clear when passing data around
    struct SSideSetReturnBuffer_10 {
        UnsignedInt32 mintNumTypes;
        UnsignedInt32 mintNumMembersIDs;
        UnsignedInt32 mintNumMembersSense8;
        UnsignedInt32 mintNumMembersSense32;
        UnsignedInt32 mintNumMembersSideNum;
        UnsignedInt32 mintNumDistFact;
        SSideSetData_10* mpaSideSetData;
        UnsignedInt32* mpaintMemberIDs;
        char* mpachrMemberSense;
        UnsignedInt32* mpaintMemberSense;
        char* mpachrMemberSideNum;
        double* mpadblDistribution;
    } mSideSetBuff_10;
    struct SSideSetReturnBuffer_11 {
        UnsignedInt32 mintNumTypes;
        UnsignedInt32 mintNumMembersIDs;
        UnsignedInt32 mintNumMembersSense;
        UnsignedInt32 mintNumMembersSideNum;
        UnsignedInt32 mintNumDistFact;
        SSideSetData_11* mpaSideSetData;
        UnsignedInt32* mpaintMemberTypes;
        UnsignedInt32* mpaintMemberIDs;
        char* mpachrMemberSense;
        UnsignedInt32 mintNumWRTEntities;
        UnsignedInt32* mpaintMemberWRTEntities;
        double* mpadblDistribution;
    } mSideSetBuff_11;
};

template <class TBuffer>
TBuffer* AdjustBuffer(UnsignedInt32 xintRequiredSize,
                      UnsignedInt32& xintActualSize, TBuffer*& xpaBuffer)
{
    if(!xintRequiredSize)  return NULL;  // Nothing requested, return nothing.

    if(xintActualSize < xintRequiredSize) {
        if(xpaBuffer)
            delete [] xpaBuffer;
        xintActualSize = xintRequiredSize;
        xpaBuffer = new TBuffer[xintActualSize];
        if(!xpaBuffer) {
            xintActualSize = 0;
            throw CCubitFile::eMemoryError;
        }
    }
    return xpaBuffer;
}

} // namespace NCubitFile

#endif

