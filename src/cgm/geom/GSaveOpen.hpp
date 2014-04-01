//- Class:          GSaveOpen
//- Description:    Save/Open tool for the geometry engine; For now,
//-                 is just used to get the proper ids.  If importing
//-                 a Cubit file into an existing model, an increment
//-                 is added to the CA_ENTITY_ID id so that everything
//-                 else can be hooked up properly.  The class CubitSaveOpen
//-                 is derived from this class.
//- Owner:          Steve Storm

#ifndef GSAVEOPEN_HPP
#define GSAVEOPEN_HPP

#include "DLIList.hpp"
#include "CubitGeomConfigure.h"

class RefEntity;

class CUBIT_GEOM_EXPORT GSaveOpen 
{
public:
    GSaveOpen();
    ~GSaveOpen();

    static int gso_sets_ids();
    static int get_id_inc( RefEntity *entity );
    static void set_error();
    static void add_error_id( int id );

protected:
   static int gsoSetsIds;

   static int gsoIncBodyId;
   static int gsoIncRefVolumeId;
   static int gsoIncRefFaceId;
   static int gsoIncRefEdgeId;
   static int gsoIncRefVertexId;

   static int gsoErrorCount;
   static DLIList<int> gsoErrorIdList;
};

inline int
GSaveOpen::gso_sets_ids()
{return gsoSetsIds;}

#endif

