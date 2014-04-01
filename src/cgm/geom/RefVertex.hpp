#ifndef REF_VERTEX_HPP
#define REF_VERTEX_HPP

#include "BasicTopologyEntity.hpp"

class RefVolume;
class RefFace;
class Point;

class CUBIT_GEOM_EXPORT RefVertex : public BasicTopologyEntity
{
public :
  
  friend class RefEntityFactory;
    //- the factory is allowed to call the (private) constructors

  virtual ~RefVertex() ;
    //- The destructor

  static const char* get_class_name()
     {
       return "Vertex";
     }

  virtual const char* class_name() const
     {
       return get_class_name();
     }
  
  DagType dag_type() const { return DagType::ref_vertex_type(); }
  const type_info& entity_type_info() const { return typeid(RefVertex); }
  
  Point* get_point_ptr() ;
  Point const* get_point_ptr() const ;
    //R Point*
    //R- A pointer to the Point to which the current 
    //R- vertex points. 
    //- This function returns a pointer to the Point
    //- which the current vertex points to.

  CubitVector coordinates () const;
    //R CubitVector
    //R- Contains the coordinate values {x y z} of this RefVertex
    //- Returns the spatial coordinates of this RefVertex.

  virtual CubitVector center_point();
    //- Return the "center" of this RefVertex (its coordinates)

  virtual int dimension() const
      { return 0; }
    //- returns the geometric dimension of RefVertex'es. 

   //int number_mesh_elements(int dimension = -1,
   //                        ElementType type = INVALID_ELEMENT_TYPE);
    //- Returns number of mesh entities owned by this geometry entity
    //- of the specified dimension. If dimension == -1, the highest 
    //- dimension elements are returned.

  void common_ref_edges( RefVertex *other_vertex,
                         DLIList <RefEdge*> &common_ref_edges,
                         RefFace *owning_face = NULL );
    //- returns all the common edges between this and the other vertex.
    //- It is possible for two vertices to be shared by more than
    //- one refedge on a surface.  Use this function rather than
    //- the next for robustness.
    //- The common edges must have this vertex and other_vertex as its
    //- start and end vertices.

  RefEdge *common_ref_edge(RefVertex *other_vertex, 
                           RefFace *owning_face = NULL);
    //- returns an edge sharing the other vertex and owned by owning_face 
    //- (if non-NULL)
    //- The edge must have this vertex and other_vertex as its
    //- start and end vertices.

  CubitBoolean about_spatially_equal( RefVertex* ref_vertex_ptr_2,
                                      double tolerance_factor = 1.0,
                                      CubitBoolean notify_ref_entity =
                                      CUBIT_FALSE );
    //R-CUBIT_TRUE/CUBIT_FALSE
    //I RefVertex*, double, CubitBoolean
    //I- Second RefVertex to compare, Tolerance factor to for GEOMETRY_RESABS,
    //I- and flag for notifying compared RefEntities.
    //O CubitBoolean
    //O- if the two vertices are spatially equal within the GEOMETRY_RESABS*
    //-the tolerance_factor, then CUBIT_TRUE will be returned.  Otherwise
    //-CUBIT_FALSE is returned.

protected :

  RefVertex(Point* pointPtr) ;
    //- Contructor: Sets up a reference vertex from a Point

private :

  void initialize();
    //- Initializes all member data

#ifdef BOYD17 
  int refVertexNbFlag;
    //- various flags and identifiers
#endif

  RefVertex( const RefVertex& );
  void operator=( const RefVertex& );

};

// ********** BEGIN HELPER CLASSES         **********
// ********** END   HELPER CLASSES         **********

// ********** BEGIN INLINE FUNCTIONS       **********
// ********** END INLINE FUNCTIONS         **********
 
// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********
 
// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********
 
#endif

