#ifndef _GEOMETRYFEATURETOOL_HPP
#define _GEOMETRYFEATURETOOL_HPP

// *** BEGIN INCLUDES *** //
#include "CubitDefines.h"
#include "DLIList.hpp"
#include "CubitGeomConfigure.h"
#include "GeometryFeatureEngine.hpp"
// *** END INCLUDES *** //

// *** BEGIN FORWARD DECLARATIONS *** //
class TopologyBridge;
class TopologyEntity;
class RefEntity;
class RefFace;
class RefEdge;
class Body;
class BasicTopologyEntity;
// *** END FORWARD DECLARATIONS *** //

//! This is a singleton pattern class for the feature functions
class CUBIT_GEOM_EXPORT GeometryFeatureTool
{
public:
    //! singleton pattern class instance interface
    //! \return
    //! the one and only instace
    static GeometryFeatureTool* instance( GeometryFeatureEngine* GFEPtr = NULL );

    //! destructor
    ~GeometryFeatureTool();

    //! add a feature engine to the list
    void add_gfe( GeometryFeatureEngine *gfe_ptr );

    //! Returns the feature type that created the entity
    GeometryFeatureEngine::FeatureType
        feature_type(BasicTopologyEntity* bte);

    //! Get all of the entities related to 'entity_in' by the same source feature
    CubitStatus get_related_by_feature(
        BasicTopologyEntity* bte,
        DLIList<GeometryEntity*>& return_list);

     //! Get the feature types available from this FeatureEngine
    DLIList<GeometryFeatureEngine::FeatureType>
        available_feature_types(BasicTopologyEntity* bte);

private:
    //! returns the feature engine of an entity
    GeometryFeatureEngine* get_engine( TopologyBridge *tb_ptr ) const;

    //! returns the feature engine of an entity
    GeometryFeatureEngine* get_engine( TopologyEntity *te_ptr ) const;

    //! determines if entities belong to the same engine
    CubitBoolean same_feature_engine( DLIList<RefEntity*> &ref_entity_list,
        CubitBoolean check_children ) const;

    //! determines if entities belong to the same engine
    CubitBoolean same_feature_engine( DLIList<TopologyEntity*> &topo_list ) const;

public:


protected:

    //! constructor
    GeometryFeatureTool( GeometryFeatureEngine *GFEPtr );

private:

    //! singleton pattern class instance interface
    static GeometryFeatureTool* instance_;

    //! list of feature engines
    DLIList<GeometryFeatureEngine*> gfeList;

};
#endif

