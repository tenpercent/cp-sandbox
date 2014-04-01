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

 
 Filename      : CubitFileIOWrapper.hpp

 Purpose       : Encapsulates file I/O operations for the Cubit file.
           
 Special Notes :

 Creator       : Will A. Helden

 Creation Date : 02/15/02

 Owner         : Will A. Helden

*******************************************************************************/

#include "CubitFile.hpp"
#include "CubitUtilConfigure.h"

#ifndef CubitFileIOWrapper_HPP
#define CubitFileIOWrapper_HPP

namespace NCubitFile {

class CUBIT_UTIL_EXPORT CIOWrapper {
public:
    CIOWrapper(FILE* xpFile,
        UnsignedInt32 xeSourceEndian = CCubitFile::mintNativeEndian);
    CIOWrapper(UnsignedInt32 swap_endian, 
	FILE* xpFile );
    CIOWrapper(FILE* xpFile,
        UnsignedInt32 xintAbsoluteOffset, UnsignedInt32 xintRelativeOffset);
    virtual ~CIOWrapper();
    
    virtual UnsignedInt32 BeginWriteBlock(UnsignedInt32 xintAbsoluteOffset = 0);
    virtual void BeginRewriteBlock(UnsignedInt32 xintAbsoluteOffset,
        UnsignedInt32 xintRelativeOffset);
    virtual void Write(const UnsignedInt32* xpaintData, UnsignedInt32 xintCount);
    virtual void Write(const char* xpachrData, UnsignedInt32 xintCount,
        UnsignedInt32 xint32bitPadded = 0);
    virtual void Write(const double* xpadblData, UnsignedInt32 xintCount);
    virtual void Write(const char* xpachrData);
    virtual UnsignedInt32 EndWriteBlock();
    
    virtual void BeginReadBlock(UnsignedInt32 xintAbsoluteOffset,
        UnsignedInt32 xintRelativeOffset = 0);
    virtual void Read(UnsignedInt32* xpaintData, UnsignedInt32 xintCount);
    virtual void Read(char* xpachrData, UnsignedInt32 xintCount,
        UnsignedInt32 xint32bitPadded = 0);
    virtual void Read(double* xpadblData, UnsignedInt32 xintCount);
    virtual char* Read();
    virtual void EndReadBlock();
    virtual UnsignedInt32 get_endian() { return mintSwapEndian; }  

private:
    FILE* mpFile;
    UnsignedInt32 mintSwapEndian;
    UnsignedInt32 mintBlockStart;
    UnsignedInt32 mintBlockEnd;
};

template <class T> void SwapEndian(unsigned int xintCount, T* xpT)
{
    int lintToByte;
    unsigned char* lpCurFromByte;
    unsigned char lachrBuffer[sizeof(T)];
    
    unsigned char* lpCurAtom = reinterpret_cast<unsigned char*>(xpT);
    int lintAtom = xintCount;
    while(lintAtom) {
        lintAtom--;
        lintToByte = sizeof(T);
        lpCurFromByte = lpCurAtom;
        while(lintToByte) {
            lintToByte--;
            lachrBuffer[lintToByte] = *lpCurFromByte;
            lpCurFromByte++;
        }
        memcpy(lpCurAtom, lachrBuffer, sizeof(T));
        lpCurAtom += sizeof(T);
    }
}

/*#define SWAPENDIAN(XINTATOMSIZE, XINTCOUNT, XPTARGET) \
    { \
        int lintToByte; \
        unsigned char* lpCurFromByte; \
        unsigned char lachrBuffer[(XINTATOMSIZE)]; \
        \
        unsigned char* lpCurAtom = (unsigned char*)(XPTARGET); \
        int lintAtom = (XINTCOUNT); \
        while(lintAtom) { \
            lintAtom--; \
            lintToByte = (XINTATOMSIZE); \
            lpCurFromByte = lpCurAtom; \
            while(lintToByte) { \
                lintToByte--; \
                lachrBuffer[lintToByte] = *lpCurFromByte; \
                lpCurFromByte++; \
            } \
            memcpy(lpCurAtom, lachrBuffer, (XINTATOMSIZE)); \
            lpCurAtom += (XINTATOMSIZE); \
        } \
    }*/

} // namespace NCubitFile

#endif

