//-------------------------------------------------------------------------
// Filename      : RefEntityFactory.hpp
//
// Purpose       : This file produces RefEntity and derived classes; it is
//                 the interface through which RE's should be created;
//                 applications of the geometry module can provide their
//                 own factory, otherwise this factory is used
//
// Creator       : Tim Tautges (in it's new state)
//
// Creation Date : 9/99
//-------------------------------------------------------------------------

#ifndef REFENTITYFACTORY_HPP
#define REFENTITYFACTORY_HPP

#include "CubitObserver.hpp"
#include "CubitDefines.h"
#include "CubitGeomConfigure.h"

class RefVertex;
class RefEdge;
class RefFace;
class RefVolume;
class RefGroup;
class Body;
class RefEntity;

class Point;
class Curve;
class Surface;
class Lump;
class BodySM;

class RefVertex;
class RefEdge;
class RefFace;
class RefGroup;
class RefVolume;
template <class X> class DLIList;
class CubitEntity;

class CUBIT_GEOM_EXPORT RefEntityFactory : public CubitObserver
{
public:


  static RefEntityFactory *instance();
    //- the function used to access the singleton instance

  virtual RefVertex *construct_RefVertex(Point *point = NULL);

  virtual RefEdge *construct_RefEdge(Curve *curve = NULL);

  virtual RefFace *construct_RefFace(Surface *surface = NULL);

  virtual RefVolume *construct_RefVolume(Lump *lump = NULL);

  virtual Body *construct_Body(BodySM *body_sm = NULL);

  virtual RefGroup *construct_RefGroup(const char* name = NULL);
  virtual RefGroup *construct_RefGroup (DLIList<RefEntity*>& entity_list);

  virtual CubitStatus ref_entity_list( char const* keyword,
                                       DLIList<RefEntity*> &entity_list,
                                       const CubitBoolean print_errors = CUBIT_TRUE);
    //- return entities of type keyword from the global list in entity_list;
    //- virtual to allow derived factories to work from their native lists

  virtual void bodies               (DLIList<Body*> &bodies);
  virtual void ref_volumes          (DLIList<RefVolume*> &ref_volumes);
  virtual void ref_groups           (DLIList<RefGroup*> &ref_groups);
  virtual void ref_faces            (DLIList<RefFace*> &ref_faces);
  virtual void ref_edges            (DLIList<RefEdge*> &ref_edges);
  virtual void ref_vertices         (DLIList<RefVertex*> &ref_vertices);

  virtual int num_bodies() const;
  virtual int num_ref_volumes() const;
  virtual int num_ref_groups() const;
  virtual int num_ref_faces() const;
  virtual int num_ref_edges() const;
  virtual int num_ref_vertices() const;

  virtual RefEntity          *get_ref_entity          (const char *type, int id);
  virtual RefEntity          *get_ref_entity          (const type_info& type, int id);
  virtual Body               *get_body                ( int id );
  virtual RefVolume          *get_ref_volume          ( int id );
  virtual RefGroup           *get_ref_group           ( int id );
  virtual RefFace            *get_ref_face            ( int id );
  virtual RefEdge            *get_ref_edge            ( int id );
  virtual RefVertex          *get_ref_vertex          ( int id );

  virtual Body *get_first_body                ();
  virtual RefVolume *get_first_ref_volume          ();
  virtual RefGroup *get_first_ref_group           ();
  virtual RefFace *get_first_ref_face            ();
  virtual RefEdge *get_first_ref_edge            ();
  virtual RefVertex *get_first_ref_vertex          ();

  virtual Body *get_next_body                ();
  virtual RefVolume *get_next_ref_volume          ();
  virtual RefGroup *get_next_ref_group           ();
  virtual RefFace *get_next_ref_face            ();
  virtual RefEdge *get_next_ref_edge            ();
  virtual RefVertex *get_next_ref_vertex          ();

  virtual Body *get_last_body                ();
  virtual RefVolume *get_last_ref_volume          ();
  virtual RefGroup *get_last_ref_group           ();
  virtual RefFace *get_last_ref_face            ();
  virtual RefEdge *get_last_ref_edge            ();
  virtual RefVertex *get_last_ref_vertex          ();

  virtual void add ( RefEntity    *ref_entity );
  virtual void add ( Body         *bodyPtr  );
  virtual void add ( RefVolume    *refVolumePtr  );
  virtual void add ( RefGroup     *refGroupPtr  );
  virtual void add ( RefFace      *refFacePtr  );
  virtual void add ( RefEdge      *refEdgePtr  );
  virtual void add ( RefVertex    *refVertexPtr  );

  virtual void remove (RefEntity    *ref_entity );
  virtual void remove (RefVertex* refVertexPtr);
  virtual void remove (RefEdge*   refEdgePtr);
  virtual void remove (RefFace*   refFacePtr);
  virtual void remove (RefGroup*  refGroupPtr);
  virtual void remove (RefVolume* refVolumePtr);
  virtual void remove (Body*      bodyPtr);

  int next_body_id   ();
  int next_ref_volume_id ();
  int next_ref_group_id  ();
  int next_ref_face_id   ();
  int next_ref_edge_id   ();
  int next_ref_vertex_id ();
  int next_surf_sub_domain_id ();
  int next_curve_sub_domain_id();

  int current_body_id() {return maxBodyId;};
  int current_volume_id() {return maxRefVolumeId;};
  int current_group_id() {return maxRefGroupId;};
  int current_face_id() {return maxRefFaceId;};
  int current_edge_id() {return maxRefEdgeId;};
  int current_vertex_id() {return maxRefVertexId;};

  void incorporate_id (RefEntity *ref_ent);
  
  int maximum_id(const char *entity_type);
  int maximum_id (RefEntity *ref_ent);
    //- get max id of specified type
  
  void maximum_id (const type_info& type, int max_id);
    //- set max id of specified type
  
  void compress_ref_ids(const char *entity_type, int retain_max_id);
    //- take out holes in the id space; make virtual so that factory can work
    //- off of native lists

  static void compress_ids(DLIList<CubitEntity*> &list);
    //- compress ids of the entities in this list; this function should really be
    //- in CubitEntity, but it references label drawing functions, which aren't
    //- accessable from CubitEntity

  void reset_ids();
    //- resets all the id counters to 0
  
  virtual void renumber_geometry_by_properties(CubitBoolean retain_max_ids);
    //- Renumbers all geometry based on geometric properties such as volume,
    //- bounding box, etc.  Used to get consistent numbering, no matter
    //- what order the entities appear in the input file.
  
  CubitStatus notify_observer(CubitObservable *observable,
                              const CubitEvent &observer_event,
                              CubitBoolean from_observable = CUBIT_FALSE);
    //- handle MODEL_ENTITY_DESTRUCTED/MODEL_ENTITY_CONSTRUCTED events

protected: 

  RefEntityFactory(RefEntityFactory *factory = NULL) ;
    //- This is a singleton class, but derived classes need to be able
    //- to instantiate it
  virtual ~RefEntityFactory ();

  static RefEntityFactory *instance_;
    //- the singleton instance
  
  int maxBodyId;
  int maxRefVolumeId;
  int maxRefGroupId;
  int maxRefFaceId;
  int maxRefEdgeId;
  int maxRefVertexId;
  int maxRefCoordSysId;
  int maxSurfSubDomainId;
  int maxCurveSubDomainId;

#ifdef BOYD17 
  DLIList<RefEntity*> refEntityList;
    //- temporary list used return generic versions of global lists
#endif

private:

  DLIList<RefVertex*> *refVertexList;
  DLIList<RefEdge*> *refEdgeList;
  DLIList<RefFace*> *refFaceList;
  DLIList<RefGroup*> *refGroupList;
  DLIList<RefVolume*> *refVolumeList;
  DLIList<Body*> *bodyList;

};

inline RefEntityFactory *RefEntityFactory::instance()
{
  if (instance_ == NULL) {
    instance_ = new RefEntityFactory();
    instance_->register_static_observer(instance_); 
  }
  
  return instance_;
}

#endif

