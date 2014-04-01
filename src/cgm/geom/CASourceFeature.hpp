#ifndef CA_SOURCE_FEATURE_HPP
#define CA_SOURCE_FEATURE_HPP

#include "CubitAttrib.hpp"
#include "DLIList.hpp"
#include "CubitDefines.h"
#include "CADefines.hpp"
#include "GeometryFeatureEngine.hpp"

class RefEntity;

//! Attribute for managing RefEntity source feature types
class CUBIT_GEOM_EXPORT CASourceFeature: public CubitAttrib
{
private:
    //! Temp source feature enum
    GeometryFeatureEngine::FeatureType sourceFeature;

public:
    //! Constructor
    CASourceFeature(RefEntity*);

    //! Constructor
    CASourceFeature(RefEntity*, CubitSimpleAttrib *);

    //! destructor
    virtual ~CASourceFeature();

    //! Returns typeid(CASourceFeature)
    virtual const type_info& entity_type_info() const
    { return typeid(CASourceFeature);}

    //! Actuate
    CubitStatus actuate();

    //! Update
    CubitStatus update();

    //! Reset function, cleans out name lists
    CubitStatus reset();

    //! Returns the simple cubit attribute for this attribute
    CubitSimpleAttrib* cubit_simple_attrib();

    //! Given a cubit attribute string return the FeatureType
    GeometryFeatureEngine::FeatureType 
        string_to_feature_type(CubitString value_in);

    //! Returning a cubit simple attribute string for the input feature type
    CubitString
        feature_type_to_string(GeometryFeatureEngine::FeatureType type_in);

    //! Return the #define attribute type
    int int_attrib_type() {return CA_SOURCE_FEATURE;}

    //! Print this attribute to the output
    void print();

};

//! global CASourceFeature creation function
CubitAttrib* CASourceFeature_creator(RefEntity* entity, CubitSimpleAttrib *p_csa);

#endif

