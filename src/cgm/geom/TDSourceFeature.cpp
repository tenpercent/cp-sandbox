#include "TDSourceFeature.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CoEdge.hpp"
#include "Loop.hpp"
#include "GMem.hpp"

TDSourceFeature::TDSourceFeature(GeometryFeatureEngine::FeatureType type_in)
{
    sourceFeature = type_in;
}


TDSourceFeature::~TDSourceFeature()
{
}
