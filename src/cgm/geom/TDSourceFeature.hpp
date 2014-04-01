//-Class: TDSourceFeature.hpp
#ifndef TD_SOURCE_FEATURE_HPP
#define TD_SOURCE_FEATURE_HPP

#include "ToolData.hpp"
#include "CastTo.hpp"
#include "CubitGeomConfigure.h"
#include "GeometryFeatureEngine.hpp"

//! Source feature tool data
class CUBIT_GEOM_EXPORT TDSourceFeature : public ToolData
{
public:

    TDSourceFeature(GeometryFeatureEngine::FeatureType type_in);

    //!Destructor
    ~TDSourceFeature();

    //! Returns the source feature on this RefEntity
    GeometryFeatureEngine::FeatureType source_feature()
    {return sourceFeature;}

    //! Sets the source feature on this RefEntity
    void source_feature(GeometryFeatureEngine::FeatureType type_in)
    {sourceFeature = type_in;}

    static int is_source_feature(const ToolData* td)
    {return (CAST_TO(td, const TDSourceFeature) != NULL);}
private:

    GeometryFeatureEngine::FeatureType sourceFeature;
};

#endif // TD_SOURCE_FEATURE_HPP
