#ifndef _GEOMETRYFEATUREENGINE_HPP
#define _GEOMETRYFEATUREENGINE_HPP

// *** BEGIN INCLUDES *** //
#include "CubitDefines.h"
#include "CubitGeomConfigure.h"
#include "DLIList.hpp"
// *** END INCLUDES *** //

// *** BEGIN FORWARD DECLARATIONS *** //
class TopologyBridge;
class TopologyEntity;
class GeometryEntity;
// *** END FORWARD DECLARATIONS *** //

//! Defines an interface for GeometryFeatureEngine
class CUBIT_GEOM_EXPORT GeometryFeatureEngine
{
public:

    enum FeatureType
    {
        FEATURE_UNDEFINED,
        FEATURE_HOLE,       
        FEATURE_ROUND,      
        FEATURE_CHAMFER,    
        FEATURE_SLOT ,      
        FEATURE_CUT,
        FEATURE_IMPRINT
    };

    //! virtual destructor
    virtual ~GeometryFeatureEngine() {}

    //! returns CUBIT_TRUE if the entity belongs to this feature engine
    virtual CubitBoolean is_feature_engine( const TopologyBridge* ) = 0;

    //! returns CUBIT_TRUE if the entity belongs to this feature engine
    virtual CubitBoolean is_feature_engine( const TopologyEntity* ) = 0;

    //! Gets the entities related to entity_in by feature 
    /*! For example, passing in a surface created by a round would get all of the surfaces
    related to that round*/
    virtual CubitStatus get_related_by_feature(
        GeometryEntity *entity_in,
        DLIList<GeometryEntity*> &list_out) = 0;

    //! Gets the features available on the entity
    /*! In other words, the feature types available from the geometry kernel*/
#if _MSC_VER <= 1200
    virtual DLIList<FeatureType>
#else
    virtual DLIList<GeometryFeatureEngine::FeatureType> 
#endif
        available_feature_types(GeometryEntity* entity) = 0;

    virtual GeometryFeatureEngine::FeatureType
        feature_type(GeometryEntity* entity) = 0;

protected:

private:

};

#endif

