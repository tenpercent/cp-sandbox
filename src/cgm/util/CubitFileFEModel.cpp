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

 
 Filename      : CubitFileFEModel.cpp

 Purpose       : Implements the reading and writing functionality for a FE
                 model section of a Cubit (*.cub) format file.
           
 Special Notes :

 Creator       : Will A. Helden

 Creation Date : 02/15/02

 Owner         : Will A. Helden

*******************************************************************************/

#include "CubitFileFEModel.hpp"
#include "CubitFileIOWrapper.hpp"

using namespace NCubitFile;

///////////////////////////////////////////////////////////////////////////////
// CFEModel
///////////////////////////////////////////////////////////////////////////////

// The number of 32 bit words contained int each of the stored structures.
const UnsignedInt32 CFEModel::mintSizeOfFEModelHeader =
    sizeof(SCubitFileFEModelHeader) / sizeof(UnsignedInt32);
const UnsignedInt32 CFEModel::mintSizeOfGeomEntry =
    sizeof(SCubitFileGeomEntry) / sizeof(UnsignedInt32);
const UnsignedInt32 CFEModel::mintSizeOfGroupEntry =
    sizeof(SCubitFileGroupEntry) / sizeof(UnsignedInt32);
const UnsignedInt32 CFEModel::mintSizeOfBlockEntry =
    sizeof(SCubitFileBlockEntry) / sizeof(UnsignedInt32);
const UnsignedInt32 CFEModel::mintSizeOfNodeSetEntry =
    sizeof(SCubitFileNodeSetEntry) / sizeof(UnsignedInt32);
const UnsignedInt32 CFEModel::mintSizeOfSideSetEntry =
    sizeof(SCubitFileSideSetEntry) / sizeof(UnsignedInt32);

CFEModel::CFEModel()
{
    mpReadFile = mpWriteFile = NULL;
    mpaGeoms = NULL;
    mpaGroups = NULL;
    mpaBlocks = NULL;
    mpaNodeSets = NULL;
    mpaSideSets = NULL;

    mintFEModelOffset = 0;
    memset(&mFEModel, 0, sizeof(SCubitFileFEModelHeader));
    mFEModel.mintFEModelEndian = CCubitFile::mintNativeEndian;

    // Initialize the return buffers for file data - all memory allocated by
    // this object should also be freed up by it.
    memset(&mNodeBuff, 0, sizeof(SNodeReturnBuffer));
    memset(&mElemBuff, 0, sizeof(SElemReturnBuffer));
    memset(&mGroupBuff, 0, sizeof(SGroupReturnBuffer));
    memset(&mBlockBuff, 0, sizeof(SBlockReturnBuffer));
    memset(&mNodeSetBuff, 0, sizeof(SNodeSetReturnBuffer));
    memset(&mSideSetBuff_10, 0, sizeof(SSideSetReturnBuffer_10));
    memset(&mSideSetBuff_11, 0, sizeof(SSideSetReturnBuffer_11));
}

CFEModel::~CFEModel()
{
    if(mpaGeoms)
        delete [] mpaGeoms;
    if(mNodeBuff.mpaNodeIDs)
        delete [] mNodeBuff.mpaNodeIDs;
    if(mNodeBuff.mpadblX)
        delete [] mNodeBuff.mpadblX;
    if(mNodeBuff.mpadblY)
        delete [] mNodeBuff.mpadblY;
    if(mNodeBuff.mpadblZ)
        delete [] mNodeBuff.mpadblZ;
    if(mElemBuff.mpaElemData)
        delete [] mElemBuff.mpaElemData;
    if(mElemBuff.mpaElemIDs)
        delete [] mElemBuff.mpaElemIDs;
    if(mElemBuff.mpaElemConnect)
        delete [] mElemBuff.mpaElemConnect;
    if(mGroupBuff.mpaGroupData)
        delete [] mGroupBuff.mpaGroupData;
    if(mGroupBuff.mpaintMemberIDs)
        delete [] mGroupBuff.mpaintMemberIDs;
    if(mBlockBuff.mpaBlockData)
        delete [] mBlockBuff.mpaBlockData;
    if(mBlockBuff.mpaintMemberIDs)
        delete [] mBlockBuff.mpaintMemberIDs;
    if(mBlockBuff.mpadblAttributes)
        delete [] mBlockBuff.mpadblAttributes;
    if(mNodeSetBuff.mpaNodeSetData)
        delete [] mNodeSetBuff.mpaNodeSetData;
    if(mNodeSetBuff.mpaintMemberIDs)
        delete [] mNodeSetBuff.mpaintMemberIDs;
    
    if(mSideSetBuff_10.mpaSideSetData)
        delete [] mSideSetBuff_10.mpaSideSetData;
    if(mSideSetBuff_10.mpaintMemberIDs)
        delete [] mSideSetBuff_10.mpaintMemberIDs; 
    if(mSideSetBuff_10.mpachrMemberSense)
        delete [] mSideSetBuff_10.mpachrMemberSense;
    if(mSideSetBuff_10.mpaintMemberSense)
        delete [] mSideSetBuff_10.mpaintMemberSense;
    if(mSideSetBuff_10.mpachrMemberSideNum)
        delete [] mSideSetBuff_10.mpachrMemberSideNum;
    if(mSideSetBuff_10.mpadblDistribution)
        delete [] mSideSetBuff_10.mpadblDistribution;
    
    if(mSideSetBuff_11.mpaSideSetData)
        delete [] mSideSetBuff_11.mpaSideSetData;
    if(mSideSetBuff_11.mpaintMemberTypes)
        delete [] mSideSetBuff_11.mpaintMemberTypes; 
    if(mSideSetBuff_11.mpaintMemberIDs)
        delete [] mSideSetBuff_11.mpaintMemberIDs; 
    if(mSideSetBuff_11.mpachrMemberSense)
        delete [] mSideSetBuff_11.mpachrMemberSense;
    if(mSideSetBuff_11.mpaintMemberWRTEntities)
        delete [] mSideSetBuff_11.mpaintMemberWRTEntities;
    if(mSideSetBuff_11.mpadblDistribution)
        delete [] mSideSetBuff_11.mpadblDistribution;
}


///////////////////////////////////////////////////////////////////////////////
// Write Functions
///////////////////////////////////////////////////////////////////////////////

UnsignedInt32 CFEModel::InitWrite(FILE* xpFile,
                                  UnsignedInt32 xintGeomCount,
                                  UnsignedInt32 xintGroupCount,
                                  UnsignedInt32 xintBlockCount,
                                  UnsignedInt32 xintNodeSetCount,
                                  UnsignedInt32 xintSideSetCount)
{
    if(mpWriteFile)  throw CCubitFile::eFileWriteError;

    mpWriteFile = xpFile;
    CIOWrapper lIO(mpWriteFile);

    // Write out the FE model header to reserve its position and size in the
    // file.
    mintFEModelOffset = lIO.BeginWriteBlock();
    lIO.Write((UnsignedInt32*)&mFEModel, mintSizeOfFEModelHeader);
    mFEModel.mintFEModelLength = lIO.EndWriteBlock();

    mFEModel.mintGeometryCount = xintGeomCount;
    if(xintGeomCount) {
        // Create a geometry entity array for storing ownership statistics and
        // initially blank it.
        mpaGeoms = new SCubitFileGeomEntry[xintGeomCount];
        memset(mpaGeoms, 0, sizeof(SCubitFileGeomEntry) * xintGeomCount);

        // Write the blank geometry table to the file to reserve its position
        // and size in the file.
        mFEModel.mintGeomTableOffset = lIO.BeginWriteBlock(mintFEModelOffset);
        lIO.Write((UnsignedInt32*)mpaGeoms,
            xintGeomCount * mintSizeOfGeomEntry);
        mFEModel.mintFEModelLength += lIO.EndWriteBlock();
    }

    mFEModel.mintGroupCount = xintGroupCount;
    if(xintGroupCount) {
        // Create a group array for storing ownership statistics and
        // initially blank it.
        mpaGroups = new SCubitFileGroupEntry[xintGroupCount];
        memset(mpaGroups, 0, sizeof(SCubitFileGroupEntry) * xintGroupCount);

        // Write the blank group table to the file to reserve its position
        // and size in the file.
        mFEModel.mintGroupTableOffset = lIO.BeginWriteBlock(mintFEModelOffset);
        lIO.Write((UnsignedInt32*)mpaGroups,
            xintGroupCount * mintSizeOfGroupEntry);
        mFEModel.mintFEModelLength += lIO.EndWriteBlock();
    }

    mFEModel.mintBlockCount = xintBlockCount;
    if(xintBlockCount) {
        // Create a block array for storing ownership statistics and
        // initially blank it.
        mpaBlocks = new SCubitFileBlockEntry[xintBlockCount];
        memset(mpaBlocks, 0, sizeof(SCubitFileBlockEntry) * xintBlockCount);

        // Write the blank block table to the file to reserve its position
        // and size in the file.
        mFEModel.mintBlockTableOffset = lIO.BeginWriteBlock(mintFEModelOffset);
        lIO.Write((UnsignedInt32*)mpaBlocks,
            xintBlockCount * mintSizeOfBlockEntry);
        mFEModel.mintFEModelLength += lIO.EndWriteBlock();
    }

    mFEModel.mintNodeSetCount = xintNodeSetCount;
    if(xintNodeSetCount) {
        // Create a node set array for storing ownership statistics and
        // initially blank it.
        mpaNodeSets = new SCubitFileNodeSetEntry[xintNodeSetCount];
        memset(mpaNodeSets, 0, sizeof(SCubitFileNodeSetEntry) * xintNodeSetCount);

        // Write the blank geometry table to the file to reserve its position
        // and size in the file.
        mFEModel.mintNodeSetTableOffset = lIO.BeginWriteBlock(mintFEModelOffset);
        lIO.Write((UnsignedInt32*)mpaNodeSets,
            xintNodeSetCount * mintSizeOfNodeSetEntry);
        mFEModel.mintFEModelLength += lIO.EndWriteBlock();
    }

    mFEModel.mintSideSetCount = xintSideSetCount;
    if(xintSideSetCount) {
        // Create a SideSet array for storing ownership statistics and
        // initially blank it.
        mpaSideSets = new SCubitFileSideSetEntry[xintSideSetCount];
        memset(mpaSideSets, 0, sizeof(SCubitFileSideSetEntry) * xintSideSetCount);

        // Write the blank geometry table to the file to reserve its position
        // and size in the file.
        mFEModel.mintSideSetTableOffset = lIO.BeginWriteBlock(mintFEModelOffset);
        lIO.Write((UnsignedInt32*)mpaSideSets,
            xintSideSetCount * mintSizeOfSideSetEntry);
        mFEModel.mintFEModelLength += lIO.EndWriteBlock();
    }

    return mintFEModelOffset;
}

void CFEModel::WriteNodes(UnsignedInt32 xintIndex,
						  UnsignedInt32 xintGeomID,
                          UnsignedInt32 xintNodeCount,
                          UnsignedInt32 *xpaintNodeIDs,
                          double *xpadblX, double *xpadblY, double *xpadblZ)
{
    if(!mpWriteFile)
        throw CCubitFile::eFileWriteError;
    if(xintIndex >= mFEModel.mintGeometryCount)
        throw CCubitFile::eNotFound;
    if(!mpaGeoms)
        throw CCubitFile::eOrderError;
    if(mpaGeoms[xintIndex].mintNodeCount)
        throw CCubitFile::eDuplicateWrite;

	mpaGeoms[xintIndex].mintGeomID = xintGeomID;

    if(xintNodeCount) {
        if(!xpaintNodeIDs || !xpadblX || !xpadblY || !xpadblZ)
            throw CCubitFile::ePassedNullPointer;

        CIOWrapper* lpIO = new CIOWrapper(mpWriteFile);

        mpaGeoms[xintIndex].mintNodeCount = xintNodeCount;
        mpaGeoms[xintIndex].mintNodeOffset =
            lpIO->BeginWriteBlock(mintFEModelOffset);

        lpIO->Write(xpaintNodeIDs, xintNodeCount);
        lpIO->Write(xpadblX, xintNodeCount);
        lpIO->Write(xpadblY, xintNodeCount);
        lpIO->Write(xpadblZ, xintNodeCount);

        mFEModel.mintFEModelLength += lpIO->EndWriteBlock();

        delete lpIO;
    }
}

void CFEModel::WriteElems(UnsignedInt32 xintIndex,
                          UnsignedInt32 xintNumTypes, SElemData* xpaElemData)
{
    if(!mpWriteFile)
        throw CCubitFile::eFileWriteError;
    if(xintIndex >= mFEModel.mintGeometryCount)
        throw CCubitFile::eNotFound;
    if(!mpaGeoms)
        throw CCubitFile::eOrderError;
    if(mpaGeoms[xintIndex].mintElemTypeCount)
        throw CCubitFile::eDuplicateWrite;

    if(xintNumTypes) {
        if(!xpaElemData)
            throw CCubitFile::ePassedNullPointer;
        UnsignedInt32 lintElem;
        for(lintElem = 0; lintElem < xintNumTypes; lintElem++) {
            if(xpaElemData[lintElem].mintElemCount) {
                if(!xpaElemData[lintElem].mpaElemIDs)
                    throw CCubitFile::ePassedNullPointer;
                if(xpaElemData[lintElem].mintElemOrder)
                    if(!xpaElemData[lintElem].mpaElemConnect)
                        throw CCubitFile::ePassedNullPointer;
            }
        }

        CIOWrapper* lpIO = new CIOWrapper(mpWriteFile);
        mpaGeoms[xintIndex].mintElemTypeCount = xintNumTypes;
        mpaGeoms[xintIndex].mintElemOffset =
            lpIO->BeginWriteBlock(mintFEModelOffset);

        for(lintElem = 0; lintElem < xintNumTypes; lintElem++) {
            if(!xpaElemData[lintElem].mintElemCount) {
                mpaGeoms[xintIndex].mintElemTypeCount--;
                continue;
            }
            mpaGeoms[xintIndex].mintElemCount +=
                xpaElemData[lintElem].mintElemCount;
        
            lpIO->Write(&xpaElemData[lintElem].mintElemType, 1);
            lpIO->Write(&xpaElemData[lintElem].mintElemOrder, 1);
            lpIO->Write(&xpaElemData[lintElem].mintElemCount, 1);
            lpIO->Write(xpaElemData[lintElem].mpaElemIDs,
                xpaElemData[lintElem].mintElemCount);
            if(xpaElemData[lintElem].mintElemOrder) {
                lpIO->Write(xpaElemData[lintElem].mpaElemConnect,
                    xpaElemData[lintElem].mintElemCount *
                    xpaElemData[lintElem].mintElemOrder);
            }
        }

        mpaGeoms[xintIndex].mintElemLength = lpIO->EndWriteBlock();
        mFEModel.mintFEModelLength += mpaGeoms[xintIndex].mintElemLength;
        delete lpIO;
    }
}

void CFEModel::WriteGroup(UnsignedInt32 xintIndex,
                          UnsignedInt32 xintGroupID,
                          UnsignedInt32 xintGroupType,
                          const char* xpachrGroupName,
                          UnsignedInt32 xintNumTypes,
                          SGroupData* xpaGroupData)
{
    if(!mpWriteFile)
        throw CCubitFile::eFileWriteError;
    if(xintIndex >= mFEModel.mintGroupCount)
        throw CCubitFile::eNotFound;
    if(!mpaGroups)
        throw CCubitFile::eOrderError;
    if(mpaGroups[xintIndex].mintMemberTypeCount)
        throw CCubitFile::eDuplicateWrite;

    if(xintNumTypes) {  
        // An empty group is valid, but an incompletely defined one is not.
        if(!xpaGroupData)
            throw CCubitFile::ePassedNullPointer;
        for(UnsignedInt32 lintGroup = 0; lintGroup < xintNumTypes; lintGroup++) {
            if(xpaGroupData[lintGroup].mintMemberCount &&
                !xpaGroupData[lintGroup].mpaMemberIDs)
                throw CCubitFile::ePassedNullPointer;
        }
    }

    mpaGroups[xintIndex].mintGroupID = xintGroupID;
    mpaGroups[xintIndex].mintGroupType = xintGroupType;
    mpaGroups[xintIndex].mintMemberTypeCount = xintNumTypes;
    if(mGroupMetaData.SetValue(xintGroupID, "NAME", xpachrGroupName) !=
        CCubitFile::eSuccess)
        throw CCubitFile::eFileWriteError;
    
    if(xintNumTypes) {
        CIOWrapper* lpIO = new CIOWrapper(mpWriteFile);
        mpaGroups[xintIndex].mintMemberOffset =
            lpIO->BeginWriteBlock(mintFEModelOffset);

        for(UnsignedInt32 lintGroup = 0; lintGroup < xintNumTypes; lintGroup++) {
            if(!xpaGroupData[lintGroup].mintMemberCount) {
                mpaGroups[xintIndex].mintMemberTypeCount--;
                continue;
            }
            mpaGroups[xintIndex].mintMemberCount +=
                xpaGroupData[lintGroup].mintMemberCount;
            
            lpIO->Write(&xpaGroupData[lintGroup].mintMemberType, 1);
            lpIO->Write(&xpaGroupData[lintGroup].mintMemberCount, 1);
            lpIO->Write(xpaGroupData[lintGroup].mpaMemberIDs,
                xpaGroupData[lintGroup].mintMemberCount);
        }
        
        mpaGroups[xintIndex].mintGroupLength = lpIO->EndWriteBlock();
        mFEModel.mintFEModelLength += mpaGroups[xintIndex].mintGroupLength;
        delete lpIO;
    }
    else {
        // An empty group does not have a data block in the file.
        mpaGroups[xintIndex].mintMemberOffset = 0;
        mpaGroups[xintIndex].mintGroupLength = 0;
    }
}

void CFEModel::WriteBlock(UnsignedInt32 xintIndex,
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
{
    if(!mpWriteFile)
        throw CCubitFile::eFileWriteError;
    if(xintIndex >= mFEModel.mintBlockCount)
        throw CCubitFile::eNotFound;
    if(!mpaBlocks)
        throw CCubitFile::eOrderError;
    if(mpaBlocks[xintIndex].mintMemberTypeCount)
        throw CCubitFile::eDuplicateWrite;

    UnsignedInt32 lintBlock;
    if(xintNumTypes) {  
        // An empty block is valid, but an incompletely defined one is not.
        if(!xpaBlockData)
            throw CCubitFile::ePassedNullPointer;
        for(lintBlock = 0; lintBlock < xintNumTypes; lintBlock++) {
            if(xpaBlockData[lintBlock].mintMemberCount &&
                !xpaBlockData[lintBlock].mpaMemberIDs)
                throw CCubitFile::ePassedNullPointer;
        }
        if(xintAttributeOrder) {
            if(!xpadblAttributes)
                throw CCubitFile::ePassedNullPointer;
        }
    }
    else if(xintAttributeOrder)
        throw CCubitFile::ePassedNullPointer;

    mpaBlocks[xintIndex].mintBlockID = xintBlockID;
    mpaBlocks[xintIndex].mintBlockElementType = xintBlockType;
    mpaBlocks[xintIndex].mintMemberTypeCount = xintNumTypes;
    mpaBlocks[xintIndex].mintBlockColor = xintBlockColor;
    mpaBlocks[xintIndex].mintBlockMixedElemType = xintMixedElemType;
    mpaBlocks[xintIndex].mintBlockDefPyramidType = xintDefPyramidType;
    mpaBlocks[xintIndex].mintBlockMaterial = xintMaterialID;
    mpaBlocks[xintIndex].mintBlockDimension = xintBlockDimension;
    mpaBlocks[xintIndex].mintAttributeOrder = xintAttributeOrder;
    
    if(xintNumTypes) {
        CIOWrapper* lpIO = new CIOWrapper(mpWriteFile);
        mpaBlocks[xintIndex].mintMemberOffset =
            lpIO->BeginWriteBlock(mintFEModelOffset);

        for(lintBlock = 0; lintBlock < xintNumTypes; lintBlock++) {
            if(!xpaBlockData[lintBlock].mintMemberCount) {
                mpaBlocks[xintIndex].mintMemberTypeCount--;
                continue;
            }
            mpaBlocks[xintIndex].mintMemberCount +=
                xpaBlockData[lintBlock].mintMemberCount;
            
            lpIO->Write(&xpaBlockData[lintBlock].mintMemberType, 1);
            lpIO->Write(&xpaBlockData[lintBlock].mintMemberCount, 1);
            lpIO->Write(xpaBlockData[lintBlock].mpaMemberIDs,
                xpaBlockData[lintBlock].mintMemberCount);
        }
        if(xintAttributeOrder)
            lpIO->Write(xpadblAttributes, xintAttributeOrder);
        
        mpaBlocks[xintIndex].mintBlockLength = lpIO->EndWriteBlock();
        mFEModel.mintFEModelLength += mpaBlocks[xintIndex].mintBlockLength;
        delete lpIO;
    }
    else {
        // An empty block does not have a data block in the file.
        mpaBlocks[xintIndex].mintMemberOffset = 0;
        mpaBlocks[xintIndex].mintBlockLength = 0;
    }
}

void CFEModel::WriteNodeSet(UnsignedInt32 xintIndex,
                            UnsignedInt32 xintNodeSetID,
                            UnsignedInt32 xintColor,
                            UnsignedInt32 xintPointSymbol,
                            UnsignedInt32 xintNumTypes,
                            SNodeSetData* xpaNodeSetData)
{
    if(!mpWriteFile)
        throw CCubitFile::eFileWriteError;
    if(xintIndex >= mFEModel.mintNodeSetCount)
        throw CCubitFile::eNotFound;
    if(!mpaNodeSets)
        throw CCubitFile::eOrderError;
    if(mpaNodeSets[xintIndex].mintMemberTypeCount)
        throw CCubitFile::eDuplicateWrite;

    UnsignedInt32 lintNodeSet;
    if(xintNumTypes) {  
        // An empty node set is valid, but an incompletely defined one is not.
        if(!xpaNodeSetData)
            throw CCubitFile::ePassedNullPointer;
        for(lintNodeSet = 0; lintNodeSet < xintNumTypes; lintNodeSet++) {
            if(xpaNodeSetData[lintNodeSet].mintMemberCount &&
                !xpaNodeSetData[lintNodeSet].mpaMemberIDs)
                throw CCubitFile::ePassedNullPointer;
        }
    }

    mpaNodeSets[xintIndex].mintNodeSetID = xintNodeSetID;
    mpaNodeSets[xintIndex].mintMemberTypeCount = xintNumTypes;
    mpaNodeSets[xintIndex].mintNodeSetPointSym = xintPointSymbol;
    mpaNodeSets[xintIndex].mintNodeSetColor = xintColor;
    
    if(xintNumTypes) {
        CIOWrapper* lpIO = new CIOWrapper(mpWriteFile);
        mpaNodeSets[xintIndex].mintMemberOffset =
            lpIO->BeginWriteBlock(mintFEModelOffset);

        for(lintNodeSet = 0; lintNodeSet < xintNumTypes; lintNodeSet++) {
            if(!xpaNodeSetData[lintNodeSet].mintMemberCount) {
                mpaNodeSets[xintIndex].mintMemberTypeCount--;
                continue;
            }
            mpaNodeSets[xintIndex].mintMemberCount +=
                xpaNodeSetData[lintNodeSet].mintMemberCount;
            
            lpIO->Write(&xpaNodeSetData[lintNodeSet].mintMemberType, 1);
            lpIO->Write(&xpaNodeSetData[lintNodeSet].mintMemberCount, 1);
            lpIO->Write(xpaNodeSetData[lintNodeSet].mpaMemberIDs,
                xpaNodeSetData[lintNodeSet].mintMemberCount);
        }
        
        mpaNodeSets[xintIndex].mintNodeSetLength = lpIO->EndWriteBlock();
        mFEModel.mintFEModelLength += mpaNodeSets[xintIndex].mintNodeSetLength;
        delete lpIO;
    }
    else {
        // An empty node set does not have a data block in the file.
        mpaNodeSets[xintIndex].mintMemberOffset = 0;
        mpaNodeSets[xintIndex].mintNodeSetLength = 0;
    }
}

void CFEModel::WriteSideSet_11(UnsignedInt32 xintIndex,
                            UnsignedInt32 xintSideSetID,
                            UnsignedInt32 xintColor,
                            UnsignedInt32 xintUseShells,
                            UnsignedInt32 xintNumTypes,
                            SSideSetData_11* xpaSideSetData,
                            UnsignedInt32 xintNumDistFact,
                            double* xpadblDistribution)
{
    if(!mpWriteFile)
        throw CCubitFile::eFileWriteError;
    if(xintIndex >= mFEModel.mintSideSetCount)
        throw CCubitFile::eNotFound;
    if(!mpaSideSets)
        throw CCubitFile::eOrderError;
    if(mpaSideSets[xintIndex].mintMemberTypeCount)
        throw CCubitFile::eDuplicateWrite;

    UnsignedInt32 lintSideSet;
    if(xintNumTypes) {  
        // An empty side set is valid, but an incompletely defined one is not.
        if(!xpaSideSetData)
            throw CCubitFile::ePassedNullPointer;
        for(lintSideSet = 0; lintSideSet < xintNumTypes; lintSideSet++) {
            if(xpaSideSetData[lintSideSet].mintMemberCount) {
                if(!xpaSideSetData[lintSideSet].mpaintMemberTypes)
                    throw CCubitFile::ePassedNullPointer;
                if(!xpaSideSetData[lintSideSet].mpaintMemberIDs)
                    throw CCubitFile::ePassedNullPointer;
                if(!xpaSideSetData[lintSideSet].mpachrMemberSenses)
                    throw CCubitFile::ePassedNullPointer;
                if(!xpaSideSetData[lintSideSet].mpaintMemberWRTEntities)
                    throw CCubitFile::ePassedNullPointer;
            }
        }
        if(xintNumDistFact) {
            if(!xpadblDistribution)
                throw CCubitFile::ePassedNullPointer;
        }
    }
    else if(xintNumDistFact)
        throw CCubitFile::ePassedNullPointer;

    mpaSideSets[xintIndex].mintSideSetID = xintSideSetID;
    mpaSideSets[xintIndex].mintMemberTypeCount = xintNumTypes;
    mpaSideSets[xintIndex].mintSideSetColor = xintColor;
    mpaSideSets[xintIndex].mintUseShells = xintUseShells;
    mpaSideSets[xintIndex].mintNumDistFact = xintNumDistFact;

    if(xintNumTypes) {
        CIOWrapper* lpIO = new CIOWrapper(mpWriteFile);
        mpaSideSets[xintIndex].mintMemberOffset =
            lpIO->BeginWriteBlock(mintFEModelOffset);

        for(lintSideSet = 0; lintSideSet < xintNumTypes; lintSideSet++) {
            if(!xpaSideSetData[lintSideSet].mintMemberCount) {
                mpaSideSets[xintIndex].mintMemberTypeCount--;
                continue;
            }
            mpaSideSets[xintIndex].mintMemberCount +=
                xpaSideSetData[lintSideSet].mintMemberCount;
            
            lpIO->Write(&xpaSideSetData[lintSideSet].mintMemberCount, 1);
            lpIO->Write(xpaSideSetData[lintSideSet].mpaintMemberTypes, 
                xpaSideSetData[lintSideSet].mintMemberCount);
            lpIO->Write(xpaSideSetData[lintSideSet].mpaintMemberIDs,
                xpaSideSetData[lintSideSet].mintMemberCount);
            lpIO->Write(xpaSideSetData[lintSideSet].mpachrMemberSenses,
                xpaSideSetData[lintSideSet].mintMemberCount, 1);
            UnsignedInt32* end = xpaSideSetData[lintSideSet].mpaintMemberWRTEntities;
            int i;
            for(i=0; i<(int)xpaSideSetData[lintSideSet].mintMemberCount; i++)
            {
              UnsignedInt32 num_wrt = *end;
              end = end + 1 + (num_wrt * 2);
            }
            UnsignedInt32 wrt_size = end - xpaSideSetData[lintSideSet].mpaintMemberWRTEntities;
            lpIO->Write(&wrt_size, 1);
            lpIO->Write(xpaSideSetData[lintSideSet].mpaintMemberWRTEntities, wrt_size);
        }
        if(xintNumDistFact)
            lpIO->Write(xpadblDistribution, xintNumDistFact);
        
        mpaSideSets[xintIndex].mintSideSetLength = lpIO->EndWriteBlock();
        mFEModel.mintFEModelLength += mpaSideSets[xintIndex].mintSideSetLength;
        delete lpIO;
    }
    else {
        // An empty side set does not have a data block in the file.
        mpaSideSets[xintIndex].mintMemberOffset = 0;
        mpaSideSets[xintIndex].mintSideSetLength = 0;
        mpaSideSets[xintIndex].mintNumDistFact = 0;
    }
}

UnsignedInt32 CFEModel::EndWrite()
{
    CIOWrapper lIO(mpWriteFile);

    if(mFEModel.mintGeometryCount) {
        lIO.BeginRewriteBlock(mintFEModelOffset, mFEModel.mintGeomTableOffset);
        lIO.Write((UnsignedInt32*)mpaGeoms,
            mFEModel.mintGeometryCount * mintSizeOfGeomEntry);
        lIO.EndWriteBlock();
    }

    if(mFEModel.mintGroupCount) {
        lIO.BeginRewriteBlock(mintFEModelOffset, mFEModel.mintGroupTableOffset);
        lIO.Write((UnsignedInt32*)mpaGroups,
            mFEModel.mintGroupCount * mintSizeOfGroupEntry);
        lIO.EndWriteBlock();
    }

    if(mFEModel.mintBlockCount) {
        lIO.BeginRewriteBlock(mintFEModelOffset, mFEModel.mintBlockTableOffset);
        lIO.Write((UnsignedInt32*)mpaBlocks,
            mFEModel.mintBlockCount * mintSizeOfBlockEntry);
        lIO.EndWriteBlock();
    }

    if(mFEModel.mintNodeSetCount) {
        lIO.BeginRewriteBlock(mintFEModelOffset, mFEModel.mintNodeSetTableOffset);
        lIO.Write((UnsignedInt32*)mpaNodeSets,
            mFEModel.mintNodeSetCount * mintSizeOfNodeSetEntry);
        lIO.EndWriteBlock();
    }

    if(mFEModel.mintSideSetCount) {
        lIO.BeginRewriteBlock(mintFEModelOffset, mFEModel.mintSideSetTableOffset);
        lIO.Write((UnsignedInt32*)mpaSideSets,
            mFEModel.mintSideSetCount * mintSizeOfSideSetEntry);
        lIO.EndWriteBlock();
    }

    UnsignedInt32 lintMetaDataLength;
    mGeomMetaData.WriteMetaData(mpWriteFile, mFEModel.mintGeomMetaDataOffset,
        lintMetaDataLength, mintFEModelOffset);
    mFEModel.mintFEModelLength += lintMetaDataLength;
    mNodeMetaData.WriteMetaData(mpWriteFile, mFEModel.mintNodeMetaDataOffset,
        lintMetaDataLength, mintFEModelOffset);
    mFEModel.mintFEModelLength += lintMetaDataLength;
    mElemMetaData.WriteMetaData(mpWriteFile, mFEModel.mintElemMetaDataOffset,
        lintMetaDataLength, mintFEModelOffset);
    mFEModel.mintFEModelLength += lintMetaDataLength;
    mGroupMetaData.WriteMetaData(mpWriteFile, mFEModel.mintGroupMetaDataOffset,
        lintMetaDataLength, mintFEModelOffset);
    mFEModel.mintFEModelLength += lintMetaDataLength;
    mBlockMetaData.WriteMetaData(mpWriteFile, mFEModel.mintBlockMetaDataOffset,
        lintMetaDataLength, mintFEModelOffset);
    mFEModel.mintFEModelLength += lintMetaDataLength;
    mNodeSetMetaData.WriteMetaData(mpWriteFile, mFEModel.mintNodeSetMetaDataOffset,
        lintMetaDataLength, mintFEModelOffset);
    mFEModel.mintFEModelLength += lintMetaDataLength;
    mSideSetMetaData.WriteMetaData(mpWriteFile, mFEModel.mintSideSetMetaDataOffset,
        lintMetaDataLength, mintFEModelOffset);
    mFEModel.mintFEModelLength += lintMetaDataLength;

    lIO.BeginRewriteBlock(mintFEModelOffset, 0);
    lIO.Write((UnsignedInt32*)&mFEModel, mintSizeOfFEModelHeader);
    lIO.EndWriteBlock();

    mpWriteFile = NULL;
    return mFEModel.mintFEModelLength;
}


///////////////////////////////////////////////////////////////////////////////
// Read Functions
///////////////////////////////////////////////////////////////////////////////

void CFEModel::InitRead(FILE* xpFile,
                        UnsignedInt32 xintAbsoluteOffset,
                        UnsignedInt32& xintGeomCount,
                        UnsignedInt32& xintGroupCount,
                        UnsignedInt32& xintBlockCount,
                        UnsignedInt32& xintNodeSetCount,
                        UnsignedInt32& xintSideSetCount)
{
    if(mpReadFile)
        throw CCubitFile::eOrderError;

    mpReadFile = xpFile;
    mintFEModelOffset = xintAbsoluteOffset;
    CIOWrapper lIO(mpReadFile, xintAbsoluteOffset, 0);

    lIO.BeginReadBlock(mintFEModelOffset, 0);
    lIO.Read((UnsignedInt32*)&mFEModel, mintSizeOfFEModelHeader);
    lIO.EndReadBlock();
    xintGeomCount = mFEModel.mintGeometryCount;
    xintGroupCount = mFEModel.mintGroupCount;
    xintBlockCount = mFEModel.mintBlockCount;
    xintNodeSetCount = mFEModel.mintNodeSetCount;
    xintSideSetCount = mFEModel.mintSideSetCount;

    // Restore the geometry definition header table if there is one.
    if(mFEModel.mintGeometryCount) {
        mpaGeoms = new SCubitFileGeomEntry[mFEModel.mintGeometryCount];

        lIO.BeginReadBlock(mintFEModelOffset, mFEModel.mintGeomTableOffset);
        lIO.Read((UnsignedInt32*)mpaGeoms,
            mFEModel.mintGeometryCount * mintSizeOfGeomEntry);
        lIO.EndReadBlock();
    }

    // Restore the group definition header table if there is one.
    if(mFEModel.mintGroupCount) {
        mpaGroups = new SCubitFileGroupEntry[mFEModel.mintGroupCount];

        lIO.BeginReadBlock(mintFEModelOffset, mFEModel.mintGroupTableOffset);
        lIO.Read((UnsignedInt32*)mpaGroups,
            mFEModel.mintGroupCount * mintSizeOfGroupEntry);
        lIO.EndReadBlock();
    }

    // Restore the block definition header table if there is one.
    if(mFEModel.mintBlockCount) {
        mpaBlocks = new SCubitFileBlockEntry[mFEModel.mintBlockCount];

        lIO.BeginReadBlock(mintFEModelOffset, mFEModel.mintBlockTableOffset);
        lIO.Read((UnsignedInt32*)mpaBlocks,
            mFEModel.mintBlockCount * mintSizeOfBlockEntry);
        lIO.EndReadBlock();
    }

    // Restore the node set definition header table if there is one.
    if(mFEModel.mintNodeSetCount) {
        mpaNodeSets = new SCubitFileNodeSetEntry[mFEModel.mintNodeSetCount];

        lIO.BeginReadBlock(mintFEModelOffset, mFEModel.mintNodeSetTableOffset);
        lIO.Read((UnsignedInt32*)mpaNodeSets,
            mFEModel.mintNodeSetCount * mintSizeOfNodeSetEntry);
        lIO.EndReadBlock();
    }

    // Restore the side set definition header table if there is one.
    if(mFEModel.mintSideSetCount) {
        mpaSideSets = new SCubitFileSideSetEntry[mFEModel.mintSideSetCount];

        lIO.BeginReadBlock(mintFEModelOffset, mFEModel.mintSideSetTableOffset);
        lIO.Read((UnsignedInt32*)mpaSideSets,
            mFEModel.mintSideSetCount * mintSizeOfSideSetEntry);
        lIO.EndReadBlock();
    }

    mGeomMetaData.ReadMetaData(mpReadFile, mintFEModelOffset,
        mFEModel.mintGeomMetaDataOffset, mFEModel.mintFEModelEndian);
    mNodeMetaData.ReadMetaData(mpReadFile, mintFEModelOffset,
        mFEModel.mintNodeMetaDataOffset, mFEModel.mintFEModelEndian);
    mElemMetaData.ReadMetaData(mpReadFile, mintFEModelOffset,
        mFEModel.mintElemMetaDataOffset, mFEModel.mintFEModelEndian);
    mGroupMetaData.ReadMetaData(mpReadFile, mintFEModelOffset,
        mFEModel.mintGroupMetaDataOffset, mFEModel.mintFEModelEndian);
    mBlockMetaData.ReadMetaData(mpReadFile, mintFEModelOffset,
        mFEModel.mintBlockMetaDataOffset, mFEModel.mintFEModelEndian);
    mNodeSetMetaData.ReadMetaData(mpReadFile, mintFEModelOffset,
        mFEModel.mintNodeSetMetaDataOffset, mFEModel.mintFEModelEndian);
    mSideSetMetaData.ReadMetaData(mpReadFile, mintFEModelOffset,
        mFEModel.mintSideSetMetaDataOffset, mFEModel.mintFEModelEndian);
}

void CFEModel::ReadNodes(UnsignedInt32 xintIndex,
						 UnsignedInt32& xintGeomID,
                         UnsignedInt32& xintNodeCount,
                         UnsignedInt32*& xpaintNodeIDs,
                         double*& xpadblX,
                         double*& xpadblY,
                         double*& xpadblZ)
{
    if(!mpReadFile)
        throw CCubitFile::eFileReadError;
    if(xintIndex >= mFEModel.mintGeometryCount)
        throw CCubitFile::eNotFound;
    if(!mpaGeoms)
        throw CCubitFile::eOrderError;

	xintGeomID = mpaGeoms[xintIndex].mintGeomID;
    xintNodeCount = mpaGeoms[xintIndex].mintNodeCount;
    if(xintNodeCount) {
        if(mNodeBuff.mintNumNodes < xintNodeCount) {
            if(mNodeBuff.mpaNodeIDs)
                delete [] mNodeBuff.mpaNodeIDs;
            if(mNodeBuff.mpadblX)
                delete [] mNodeBuff.mpadblX;
            if(mNodeBuff.mpadblY)
                delete [] mNodeBuff.mpadblY;
            if(mNodeBuff.mpadblZ)
                delete [] mNodeBuff.mpadblZ;
            mNodeBuff.mintNumNodes = xintNodeCount;
            mNodeBuff.mpaNodeIDs = new UnsignedInt32[xintNodeCount];
            mNodeBuff.mpadblX = new double[xintNodeCount];
            mNodeBuff.mpadblY = new double[xintNodeCount];
            mNodeBuff.mpadblZ = new double[xintNodeCount];
        }
        xpaintNodeIDs = mNodeBuff.mpaNodeIDs;
        xpadblX = mNodeBuff.mpadblX;
        xpadblY = mNodeBuff.mpadblY;
        xpadblZ = mNodeBuff.mpadblZ;

        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mFEModel.mintFEModelEndian);
    
        lpIO->BeginReadBlock(mintFEModelOffset,
            mpaGeoms[xintIndex].mintNodeOffset);
        lpIO->Read(mNodeBuff.mpaNodeIDs, xintNodeCount);
        lpIO->Read(mNodeBuff.mpadblX, xintNodeCount);
        lpIO->Read(mNodeBuff.mpadblY, xintNodeCount);
        lpIO->Read(mNodeBuff.mpadblZ, xintNodeCount);
        lpIO->EndReadBlock();

        delete lpIO;
    }
    else {
        xpaintNodeIDs = NULL;
        xpadblX = xpadblY = xpadblZ = NULL;
    }
}

void CFEModel::ReadElems(UnsignedInt32 xintIndex, UnsignedInt32& xintGeomID,
                         UnsignedInt32& xintNumTypes, SElemData*& xpaElemData)
{
    if(!mpReadFile)
        throw CCubitFile::eFileReadError;
    if(xintIndex >= mFEModel.mintGeometryCount)
        throw CCubitFile::eNotFound;
    if(!mpaGeoms)
        throw CCubitFile::eOrderError;

	xintGeomID = mpaGeoms[xintIndex].mintGeomID;
    xintNumTypes = mpaGeoms[xintIndex].mintElemTypeCount;
    if(xintNumTypes) {
        // Resize the element return buffer if necessary and then set the return
        // pointers to the buffer.
        UnsignedInt32 lintConnGuess = 
            mpaGeoms[xintIndex].mintElemLength / sizeof(UnsignedInt32) -
            mpaGeoms[xintIndex].mintElemCount;
            // - 3 * mpaGeoms[xintIndex].mintElemTypeCount; (exact)
        xpaElemData = AdjustBuffer(xintNumTypes,
            mElemBuff.mintNumTypes, mElemBuff.mpaElemData);
        UnsignedInt32* lpIDs =
            AdjustBuffer(mpaGeoms[xintIndex].mintElemCount,
            mElemBuff.mintNumElems, mElemBuff.mpaElemIDs);
        UnsignedInt32* lpConnect =
            AdjustBuffer(lintConnGuess,
            mElemBuff.mintNumConnect, mElemBuff.mpaElemConnect);

        // Read the element from the file.
        UnsignedInt32 lintNumElems, lintNumConnect;
        UnsignedInt32 lintTotalElems = 0, lintTotalConnect = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mFEModel.mintFEModelEndian);

        lpIO->BeginReadBlock(mintFEModelOffset,
            mpaGeoms[xintIndex].mintElemOffset);
        for(UnsignedInt32 lintType = 0; lintType < xintNumTypes; lintType++) {
            lpIO->Read(&xpaElemData[lintType].mintElemType, 1);
            lpIO->Read(&xpaElemData[lintType].mintElemOrder, 1);
            lpIO->Read(&xpaElemData[lintType].mintElemCount, 1);

            xpaElemData[lintType].mpaElemIDs = lpIDs;
            xpaElemData[lintType].mpaElemConnect = lpConnect;
            lintNumElems = xpaElemData[lintType].mintElemCount;
            lintNumConnect = lintNumElems * xpaElemData[lintType].mintElemOrder;
            // Make sure the total number of elements or connection entries
            // does not exceed what was specified in the geometry table entry.
            lintTotalElems += lintNumElems;
            lintTotalConnect += lintNumConnect;
            if((lintTotalElems > mpaGeoms[xintIndex].mintElemCount) ||
                (lintTotalConnect > lintConnGuess))
                throw CCubitFile::eFileReadError;

            lpIO->Read(xpaElemData[lintType].mpaElemIDs, lintNumElems);
            if(xpaElemData[lintType].mintElemOrder)
                lpIO->Read(xpaElemData[lintType].mpaElemConnect, lintNumConnect);

            lpIDs = &lpIDs[lintNumElems];
            lpConnect = &lpConnect[lintNumConnect];
        }
        lpIO->EndReadBlock();
        delete lpIO;
    }
    else
        xpaElemData = NULL;
}

void CFEModel::ReadGroupIdentity(UnsignedInt32 xintIndex,
                                 UnsignedInt32& xintGroupID,
                                 UnsignedInt32& xintGroupType,
                                 const char*& xpachrGroupName)
{
    if(!mpReadFile)
        throw CCubitFile::eFileReadError;
    if(xintIndex >= mFEModel.mintGroupCount)
        throw CCubitFile::eNotFound;
    if(!mpaGroups)
        throw CCubitFile::eOrderError;

    xintGroupID = mpaGroups[xintIndex].mintGroupID;
    xintGroupType = mpaGroups[xintIndex].mintGroupType;
    if(mGroupMetaData.GetValue(xintGroupID, "NAME", xpachrGroupName) !=
        CCubitFile::eSuccess)
        xpachrGroupName = NULL;
}

void CFEModel::ReadGroupMembers(UnsignedInt32 xintIndex,
                                UnsignedInt32& xintNumTypes,
                                SGroupData*& xpaGroupData)
{
    if(!mpReadFile)
        throw CCubitFile::eFileReadError;
    if(xintIndex >= mFEModel.mintGroupCount)
        throw CCubitFile::eNotFound;
    if(!mpaGroups)
        throw CCubitFile::eOrderError;

    xintNumTypes = mpaGroups[xintIndex].mintMemberTypeCount;
    if(xintNumTypes) {
        // Resize the group return buffer if necessary and then set the return
        // pointers to the buffer.
        xpaGroupData = AdjustBuffer(xintNumTypes,
            mGroupBuff.mintNumTypes, mGroupBuff.mpaGroupData);
        UnsignedInt32* lpIDs =
            AdjustBuffer(mpaGroups[xintIndex].mintMemberCount,
            mGroupBuff.mintNumMembers, mGroupBuff.mpaintMemberIDs);

        // Read the group member data from the file.
        UnsignedInt32 lintNumMembers, lintTotalMembers = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mFEModel.mintFEModelEndian);
        lpIO->BeginReadBlock(mintFEModelOffset,
            mpaGroups[xintIndex].mintMemberOffset);
        for(UnsignedInt32 lintType = 0; lintType < xintNumTypes; lintType++) {
            lpIO->Read(&xpaGroupData[lintType].mintMemberType, 1);
            lpIO->Read(&xpaGroupData[lintType].mintMemberCount, 1);

            xpaGroupData[lintType].mpaMemberIDs = lpIDs;
            lintNumMembers = xpaGroupData[lintType].mintMemberCount;
            // Make sure the total number of group members does not exceed what
            // was specified in the group table entry.
            lintTotalMembers += lintNumMembers;
            if(lintTotalMembers > mpaGroups[xintIndex].mintMemberCount)
                throw CCubitFile::eFileReadError;

            lpIO->Read(xpaGroupData[lintType].mpaMemberIDs, lintNumMembers);

            lpIDs = &lpIDs[lintNumMembers];
        }
        lpIO->EndReadBlock();
        delete lpIO;
    }
    else
        xpaGroupData = NULL;
}

void CFEModel::ReadBlock(UnsignedInt32 xintIndex,
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
{
    if(!mpReadFile)
        throw CCubitFile::eFileReadError;
    if(xintIndex >= mFEModel.mintBlockCount)
        throw CCubitFile::eNotFound;
    if(!mpaBlocks)
        throw CCubitFile::eOrderError;

    xintBlockID = mpaBlocks[xintIndex].mintBlockID;
    xintBlockType = mpaBlocks[xintIndex].mintBlockElementType;
    xintNumTypes = mpaBlocks[xintIndex].mintMemberTypeCount;
    xintBlockColor = mpaBlocks[xintIndex].mintBlockColor;
    xintMixedElemType = mpaBlocks[xintIndex].mintBlockMixedElemType;
    xintDefPyramidType = mpaBlocks[xintIndex].mintBlockDefPyramidType;
    xintMaterialID = mpaBlocks[xintIndex].mintBlockMaterial;
    xintBlockDimension = mpaBlocks[xintIndex].mintBlockDimension;
    xintAttributeOrder = mpaBlocks[xintIndex].mintAttributeOrder;

    if(xintNumTypes) {
        // Resize the block return buffer if necessary and then set the return
        // pointers to the buffer.
        xpaBlockData = AdjustBuffer(xintNumTypes,
            mBlockBuff.mintNumTypes, mBlockBuff.mpaBlockData);
        UnsignedInt32* lpIDs =
            AdjustBuffer(mpaBlocks[xintIndex].mintMemberCount,
            mBlockBuff.mintNumMembers, mBlockBuff.mpaintMemberIDs);
        xpadblAttributes = AdjustBuffer(xintAttributeOrder,
            mBlockBuff.mintAttributeOrder, mBlockBuff.mpadblAttributes);

        // Read the block from the file.
        UnsignedInt32 lintNumMembers, lintTotalMembers = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mFEModel.mintFEModelEndian);
        lpIO->BeginReadBlock(mintFEModelOffset,
            mpaBlocks[xintIndex].mintMemberOffset);
        for(UnsignedInt32 lintType = 0; lintType < xintNumTypes; lintType++) {
            lpIO->Read(&xpaBlockData[lintType].mintMemberType, 1);
            lpIO->Read(&xpaBlockData[lintType].mintMemberCount, 1);

            xpaBlockData[lintType].mpaMemberIDs = lpIDs;
            lintNumMembers = xpaBlockData[lintType].mintMemberCount;
            // Make sure the total number of block members does not exceed what
            // was specified in the block table entry.
            lintTotalMembers += lintNumMembers;
            if(lintTotalMembers > mpaBlocks[xintIndex].mintMemberCount)
                throw CCubitFile::eFileReadError;

            lpIO->Read(xpaBlockData[lintType].mpaMemberIDs, lintNumMembers);

            lpIDs = &lpIDs[lintNumMembers];
        }
        if(xintAttributeOrder)
            lpIO->Read(xpadblAttributes, xintAttributeOrder);
        lpIO->EndReadBlock();
        delete lpIO;
    }
    else {
        xpaBlockData = NULL;
        xpadblAttributes = NULL;
    }
}

void CFEModel::ReadNodeSet(UnsignedInt32 xintIndex,
                           UnsignedInt32& xintNodeSetID,
                           UnsignedInt32& xintColor,
                           UnsignedInt32& xintPointSymbol,
                           UnsignedInt32& xintNumTypes,
                           SNodeSetData*& xpaNodeSetData)
{
    if(!mpReadFile)
        throw CCubitFile::eFileReadError;
    if(xintIndex >= mFEModel.mintNodeSetCount)
        throw CCubitFile::eNotFound;
    if(!mpaNodeSets)
        throw CCubitFile::eOrderError;

    xintNodeSetID = mpaNodeSets[xintIndex].mintNodeSetID;
    xintNumTypes = mpaNodeSets[xintIndex].mintMemberTypeCount;
    xintPointSymbol = mpaNodeSets[xintIndex].mintNodeSetPointSym;
    xintColor = mpaNodeSets[xintIndex].mintNodeSetColor;

    if(xintNumTypes) {
        // Resize the node set return buffer if necessary and then set the return
        // pointers to the buffer.
        xpaNodeSetData = AdjustBuffer(xintNumTypes,
            mNodeSetBuff.mintNumTypes, mNodeSetBuff.mpaNodeSetData);
        UnsignedInt32* lpIDs =
            AdjustBuffer(mpaNodeSets[xintIndex].mintMemberCount,
            mNodeSetBuff.mintNumMembers, mNodeSetBuff.mpaintMemberIDs);

        // Read the node set from the file.
        UnsignedInt32 lintNumMembers, lintTotalMembers = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mFEModel.mintFEModelEndian);
        lpIO->BeginReadBlock(mintFEModelOffset,
            mpaNodeSets[xintIndex].mintMemberOffset);
        for(UnsignedInt32 lintType = 0; lintType < xintNumTypes; lintType++) {
            lpIO->Read(&xpaNodeSetData[lintType].mintMemberType, 1);
            lpIO->Read(&xpaNodeSetData[lintType].mintMemberCount, 1);

            xpaNodeSetData[lintType].mpaMemberIDs = lpIDs;
            lintNumMembers = xpaNodeSetData[lintType].mintMemberCount;
            // Make sure the total number of node set members does not exceed
            // what was specified in the node set table entry.
            lintTotalMembers += lintNumMembers;
            if(lintTotalMembers > mpaNodeSets[xintIndex].mintMemberCount)
                throw CCubitFile::eFileReadError;

            lpIO->Read(xpaNodeSetData[lintType].mpaMemberIDs, lintNumMembers);

            lpIDs = &lpIDs[lintNumMembers];
        }
        lpIO->EndReadBlock();
        delete lpIO;
    }
    else
        xpaNodeSetData = NULL;
}

void CFEModel::ReadSideSet_10(UnsignedInt32 xintIndex,
                           UnsignedInt32& xintSideSetID,
                           UnsignedInt32& xintColor,
                           UnsignedInt32& xintUseShells,
                           UnsignedInt32& xintNumTypes,
                           SSideSetData_10*& xpaSideSetData,
                           UnsignedInt32& xintNumDistFact,
                           double*& xpadblDistribution)
{
    if(!mpReadFile)
        throw CCubitFile::eFileReadError;
    if(xintIndex >= mFEModel.mintSideSetCount)
        throw CCubitFile::eNotFound;
    if(!mpaSideSets)
        throw CCubitFile::eOrderError;

    xintSideSetID = mpaSideSets[xintIndex].mintSideSetID;
    xintNumTypes = mpaSideSets[xintIndex].mintMemberTypeCount;
    xintColor = mpaSideSets[xintIndex].mintSideSetColor;
    xintUseShells = mpaSideSets[xintIndex].mintUseShells;
    xintNumDistFact = mpaSideSets[xintIndex].mintNumDistFact;

    if(xintNumTypes) {
        // Resize the side set return buffer if necessary and then set the return
        // pointers to the buffer.
        xpaSideSetData = AdjustBuffer(xintNumTypes,
            mSideSetBuff_10.mintNumTypes, mSideSetBuff_10.mpaSideSetData);
        UnsignedInt32* lpaintIDs =
            AdjustBuffer(mpaSideSets[xintIndex].mintMemberCount,
            mSideSetBuff_10.mintNumMembersIDs, mSideSetBuff_10.mpaintMemberIDs);
        char* lpachrSense =
            AdjustBuffer(mpaSideSets[xintIndex].mintMemberCount,
            mSideSetBuff_10.mintNumMembersSense8, mSideSetBuff_10.mpachrMemberSense);
        UnsignedInt32* lpaintSense =
            AdjustBuffer(mpaSideSets[xintIndex].mintMemberCount,
            mSideSetBuff_10.mintNumMembersSense32, mSideSetBuff_10.mpaintMemberSense);
        char* lpachrSideNum =
            AdjustBuffer(mpaSideSets[xintIndex].mintMemberCount,
            mSideSetBuff_10.mintNumMembersSideNum, mSideSetBuff_10.mpachrMemberSideNum);
        xpadblDistribution = AdjustBuffer(xintNumDistFact,
            mSideSetBuff_10.mintNumDistFact, mSideSetBuff_10.mpadblDistribution);

        // Read the block from the file.
        UnsignedInt32 lintNumMembers, lintTotalMembers = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mFEModel.mintFEModelEndian);
        lpIO->BeginReadBlock(mintFEModelOffset,
            mpaSideSets[xintIndex].mintMemberOffset);
        for(UnsignedInt32 lintType = 0; lintType < xintNumTypes; lintType++) {
            lpIO->Read(&xpaSideSetData[lintType].mintMemberType, 1);
            lpIO->Read(&xpaSideSetData[lintType].mintMemberCount, 1);
            lpIO->Read(&xpaSideSetData[lintType].mintMemberSenseSize, 1);

            xpaSideSetData[lintType].mpaintMemberIDs = lpaintIDs;
            xpaSideSetData[lintType].mpachrMemberSideNumber = lpachrSideNum;
            lintNumMembers = xpaSideSetData[lintType].mintMemberCount;
            // Make sure the total number of side set members does not exceed
            // what was specified in the side set table entry.
            lintTotalMembers += lintNumMembers;
            if(lintTotalMembers > mpaSideSets[xintIndex].mintMemberCount)
                throw CCubitFile::eFileReadError;

            lpIO->Read(xpaSideSetData[lintType].mpaintMemberIDs, lintNumMembers);
            switch(xpaSideSetData[lintType].mintMemberSenseSize) {
            case CCubitFile::eSideSetSenseNone:
            default:
                xpaSideSetData[lintType].mpachrMemberSense = NULL;
                xpaSideSetData[lintType].mpaintMemberSense = NULL;
                break;
            case CCubitFile::eSideSetSenseByte:
                xpaSideSetData[lintType].mpachrMemberSense = lpachrSense;
                xpaSideSetData[lintType].mpaintMemberSense = NULL;
                lpIO->Read(xpaSideSetData[lintType].mpachrMemberSense,
                    lintNumMembers, 1);
                lpaintSense = &lpaintSense[lintNumMembers];
                break;
            case CCubitFile::eSideSetSenseInt32:
                xpaSideSetData[lintType].mpachrMemberSense = NULL;
                xpaSideSetData[lintType].mpaintMemberSense = lpaintSense;
                lpIO->Read(xpaSideSetData[lintType].mpaintMemberSense,
                    lintNumMembers);
                lpaintSense = &lpaintSense[lintNumMembers];
                break;
            }
            lpIO->Read(xpaSideSetData[lintType].mpachrMemberSideNumber,
                lintNumMembers, 1);

            lpaintIDs = &lpaintIDs[lintNumMembers];
            lpachrSideNum = &lpachrSideNum[lintNumMembers];
        }
        if(xintNumDistFact)
            lpIO->Read(xpadblDistribution, xintNumDistFact);
        lpIO->EndReadBlock();
        delete lpIO;
    }
    else {
        xpaSideSetData = NULL;
        xpadblDistribution = NULL;
    }
}

void CFEModel::ReadSideSet_11(UnsignedInt32 xintIndex,
                           UnsignedInt32& xintSideSetID,
                           UnsignedInt32& xintColor,
                           UnsignedInt32& xintUseShells,
                           UnsignedInt32& xintNumTypes,
                           SSideSetData_11*& xpaSideSetData,
                           UnsignedInt32& xintNumDistFact,
                           double*& xpadblDistribution)
{
    if(!mpReadFile)
        throw CCubitFile::eFileReadError;
    if(xintIndex >= mFEModel.mintSideSetCount)
        throw CCubitFile::eNotFound;
    if(!mpaSideSets)
        throw CCubitFile::eOrderError;

    xintSideSetID = mpaSideSets[xintIndex].mintSideSetID;
    xintNumTypes = mpaSideSets[xintIndex].mintMemberTypeCount;
    xintColor = mpaSideSets[xintIndex].mintSideSetColor;
    xintUseShells = mpaSideSets[xintIndex].mintUseShells;
    xintNumDistFact = mpaSideSets[xintIndex].mintNumDistFact;

    if(xintNumTypes) {
        // Resize the side set return buffer if necessary and then set the return
        // pointers to the buffer.
        xpaSideSetData = AdjustBuffer(xintNumTypes,
            mSideSetBuff_11.mintNumTypes, mSideSetBuff_11.mpaSideSetData);
        UnsignedInt32 count = mSideSetBuff_11.mintNumMembersIDs;
        UnsignedInt32* lpaintTypes =
            AdjustBuffer(mpaSideSets[xintIndex].mintMemberCount,
            mSideSetBuff_11.mintNumMembersIDs, mSideSetBuff_11.mpaintMemberTypes);
        UnsignedInt32* lpaintIDs =
            AdjustBuffer(mpaSideSets[xintIndex].mintMemberCount,
            count, mSideSetBuff_11.mpaintMemberIDs);
        char* lpachrSense =
            AdjustBuffer(mpaSideSets[xintIndex].mintMemberCount,
            mSideSetBuff_11.mintNumMembersSense, mSideSetBuff_11.mpachrMemberSense);
        xpadblDistribution = AdjustBuffer(xintNumDistFact,
            mSideSetBuff_11.mintNumDistFact, mSideSetBuff_11.mpadblDistribution);

        // Read the block from the file.
        UnsignedInt32 lintNumMembers, lintTotalMembers = 0;
        CIOWrapper* lpIO = new CIOWrapper(mpReadFile, mFEModel.mintFEModelEndian);
        lpIO->BeginReadBlock(mintFEModelOffset,
            mpaSideSets[xintIndex].mintMemberOffset);
        for(UnsignedInt32 lintType = 0; lintType < xintNumTypes; lintType++) {
            lpIO->Read(&xpaSideSetData[lintType].mintMemberCount, 1);
            
            lintNumMembers = xpaSideSetData[lintType].mintMemberCount;
            
            xpaSideSetData[lintType].mpaintMemberTypes = lpaintTypes;
            xpaSideSetData[lintType].mpaintMemberIDs = lpaintIDs;
            
            // Make sure the total number of side set members does not exceed
            // what was specified in the side set table entry.
            lintTotalMembers += lintNumMembers;
            if(lintTotalMembers > mpaSideSets[xintIndex].mintMemberCount)
                throw CCubitFile::eFileReadError;

            lpIO->Read(xpaSideSetData[lintType].mpaintMemberTypes, lintNumMembers);

            lpIO->Read(xpaSideSetData[lintType].mpaintMemberIDs, lintNumMembers);

            xpaSideSetData[lintType].mpachrMemberSenses = lpachrSense;
            lpIO->Read(xpaSideSetData[lintType].mpachrMemberSenses,
                lintNumMembers, 1);

            UnsignedInt32 wrt_size=0;
            lpIO->Read(&wrt_size, 1);
        
            UnsignedInt32* lpaintWRT =
                AdjustBuffer(wrt_size, mSideSetBuff_11.mintNumWRTEntities, mSideSetBuff_11.mpaintMemberWRTEntities);

            xpaSideSetData[lintType].mpaintMemberWRTEntities = lpaintWRT;
        
            lpIO->Read(xpaSideSetData[lintType].mpaintMemberWRTEntities, wrt_size);

            lpaintIDs = &lpaintIDs[lintNumMembers];
            lpaintTypes = &lpaintTypes[lintNumMembers];
        }
        if(xintNumDistFact)
            lpIO->Read(xpadblDistribution, xintNumDistFact);
        lpIO->EndReadBlock();
        delete lpIO;
    }
    else {
        xpaSideSetData = NULL;
        xpadblDistribution = NULL;
    }
}

void CFEModel::EndRead()
{
    mpReadFile = NULL;
}

///////////////////////////////////////////////////////////////////////////////
// Meta-Data Functions
///////////////////////////////////////////////////////////////////////////////

CMetaData& CFEModel::GetGeomMetaData()
{
    return mGeomMetaData;
}

CMetaData& CFEModel::GetNodeMetaData()
{
    return mNodeMetaData;
}

CMetaData& CFEModel::GetElemMetaData()
{
    return mElemMetaData;
}

CMetaData& CFEModel::GetGroupMetaData()
{
    return mGroupMetaData;
}

CMetaData& CFEModel::GetBlockMetaData()
{
    return mBlockMetaData;
}

CMetaData& CFEModel::GetNodeSetMetaData()
{
    return mNodeSetMetaData;
}

CMetaData& CFEModel::GetSideSetMetaData()
{
    return mSideSetMetaData;
}
