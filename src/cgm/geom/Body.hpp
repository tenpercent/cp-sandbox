/*-------------------------------------------------------------------------
 * Filename      : Body.hpp
 *
 * Purpose       : This class represents the topological entity, Body,
 *                 which is the highest level abstraction for a complete
 *                 geometric model.
 *
 * Special Notes : Body is a GroupingEntity in the TopologyEntity hierarchy.
 *                 It is also a RefEntity and get a lot of its functionality
 *                 and data from the RefEntity base class for meshing
 *                 purposes.
 *
 *                 The user is provided a handle to Bodies.
 *
 * Creator       : 
 *
 * Creation Date : 
 *
 * Owner         : 
 *-------------------------------------------------------------------------
 */

#ifndef BODY_HPP
#define BODY_HPP

#include "GroupingEntity.hpp"
#include "RefEntity.hpp"
#include "CoVolume.hpp"
#include "RefVolume.hpp"

class CubitTransformMatrix;

class CUBIT_GEOM_EXPORT Body : public GroupingEntity,
             public RefEntity
{
public :

  /*- the factory is allowed to call the (private) constructors */
  friend class RefEntityFactory;

    /* constructors/destructors */

  virtual ~Body() ;
  /*- The destructor. */

    /* topology */
  
  DagType dag_type() const { return DagType::body_type(); }
  const type_info& entity_type_info() const { return typeid(Body); }

  static const char* get_class_name()
     {
       return "Body";
     }
  virtual const char* class_name() const
     {
       return get_class_name();
     }
  
  BodySM* get_body_sm_ptr() const;

  virtual CubitBox bounding_box();
  /*- Returns the bounding box of this body
   */

  virtual CubitVector center_point();
  /*- Return a CubitVector set to the centroid of this body
   */

  CubitBoolean get_mass_props(CubitVector& cofg);
    /* Get certain mass props, cofg is center of gravity */
       
  CubitPointContainment point_containment( CubitVector &point );
    //- Determines whether a point is inside, outside, or on boundary
    //- of a volume.
  
  virtual double measure();
  CubitBoolean is_sheet_body();

#ifdef BOYD14
  CubitStatus get_transforms(CubitTransformMatrix &tfm);
  /*- returns the body transformation matrix
   */
#endif
  
#ifdef BOYD14
  CubitStatus transform_position(CubitVector const& position_vector,
                                 CubitVector & transformed_vector);
  /*R CubitStatus
   *R- CUBIT_SUCCESS/FAILURE
   *I position_vector
   *I- The input vector that needs to be transformed.
   *O transformed_vector
   *O- The transformed CubitVector -- sent in by reference.
   *- This function will transform the input position_vector using
   *- the transformation matrix associated with the underlying geometry
   *- of this Body.  If there is none, then CUBIT_FAILURE is returned
   *- and transformed_vector is left untouched.
   *- Note: The input vector *must* have 3 components.
   */
#endif

  virtual int validate();
  /*- do a measure and api entity check.
   */
  
  virtual CubitString measure_label();
  /*- Functions related to "measuring" the Body
   */

  void copiedFromId( int Id );
  int copiedFromId ();
  /*- Sets/Gets the id of the body that it was copied from */


  virtual void color(int value);
  virtual int color() const;

protected: 

  Body() ;
  /*- The default constructor.
   */

  Body(BodySM* OSMEPtr) ;
  /*- The constructor with a pointer to an other solid model entity.
   */

private:
   
#ifdef BOYD14
  void notify_model(CubitEventType eventType);
  /*R void
   *I eventType
   *I- The type of event that occurred 
   *- This function notifies the owning Model of the event that
   *- occured to this Body (see the EventType enum for the list of
   *- recognized events).
   */
#endif

  int copied_from_body_id;
  
  Body( const Body& );
  void operator=( const Body& );

} ;

#endif

