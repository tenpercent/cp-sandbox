#ifndef INTERMEDIATE_GEOM_ENGINE_HPP
#define INTERMEDIATE_GEOM_ENGINE_HPP

template <class X> class DLIList;
class TopologyBridge;
class CubitTransformMatrix;
class Body;
class BodySM;
class TBOwner;

class IntermediateGeomEngine
{

public:
  
  virtual ~IntermediateGeomEngine() {}

  virtual bool is_composite(TBOwner *bridge_owner) = 0;
  virtual bool is_partition(TBOwner *bridge_owner) = 0;

  virtual int level() const = 0;
  
  virtual void remove_imprint_attributes_after_modify
                                ( DLIList<BodySM*> &old_sms,
                                  DLIList<BodySM*> &new_sms )=0;
  virtual void push_imprint_attributes_before_modify
                                ( DLIList<BodySM*> &body_sms ) = 0;
  virtual CubitStatus export_geometry( DLIList<TopologyBridge*>& geometry_list ) = 0;
  
  virtual CubitStatus import_geometry( DLIList<TopologyBridge*>& geometry_list ) = 0;
  
  virtual void clean_out_deactivated_geometry() = 0;

  virtual void remove_attributes( DLIList<TopologyBridge*> &bridge_list ) = 0;
  virtual void attribute_after_imprinting( DLIList<TopologyBridge*> &new_tbs,
                                                    DLIList<TopologyBridge*> &att_tbs,
                                                    DLIList<BodySM*> &new_sms,
                                                        DLIList<Body*> &old_bodies)=0;
  
  virtual void remove_attributes_from_unmodifed_virtual(DLIList<TopologyBridge*> &bridges) = 0;
  virtual void remove_modified(DLIList<TopologyBridge*>& geometry_list) = 0;
  virtual CubitStatus notify_transform( TopologyBridge* entity,
                                        const CubitTransformMatrix& xform ) = 0;
};

#endif

