//-------------------------------------------------------------------------
// Filename      : RefVolume.cpp
//
// Purpose       : This file contains the implementation of the class 
//                 RefVolume. 
//
// Special Notes :
//
// Creator       : 
//
// Creation Date : 
//
// Owner         : 
//-------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "CubitDefines.h"
#include "CubitMessage.hpp"
#include "CubitBox.hpp"

#include "Body.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CoVolume.hpp"
#include "Shell.hpp"

// tools
#include "RefEntityFactory.hpp"
#include "GeometryQueryTool.hpp"

// lists
#include "DLIList.hpp"

#include "Lump.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor with a pointer to a Lump 
//
// Special Notes :
//
// Creator       : Xuchen Liu
//
// Creation Date : 07/11/96
//-------------------------------------------------------------------------
RefVolume::RefVolume(Lump* lumpPtr)
{
   // Set the GeometryEntity pointer
   if (lumpPtr != NULL)
   {
      set_geometry_entity_ptr(lumpPtr) ;

      // Commented-out code from Phil Tuchinsky for FUTURE meshing
      //   from a surface mesh with no solid model.
      //
//    hasGeometry = CUBIT_TRUE;
//
      // initialize the surface mesh entity lists exactly once
//    nodes_boundary  ( surfaceNodes );
//    edges_inclusive ( surfaceEdges );
//    faces_inclusive ( surfaceFaces );
   }
   else
   {
      PRINT_ERROR("In the RefVolume(Lump*) constructor\n");
      PRINT_ERROR("       Input Lump pointer is NULL\n");
      assert(CUBIT_FALSE);
   }
   
   // Initialize the member data
   initialize();
}

//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------
RefVolume::~RefVolume()
{
}

//-------------------------------------------------------------------------
// Purpose       : Return a pointer to the Lump associated with a volume.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------

Lump* RefVolume::get_lump_ptr() 
{
  return CAST_TO(get_geometry_entity_ptr(), Lump);
}

Lump const* RefVolume::get_lump_ptr() const 
{
  return CAST_TO(get_geometry_entity_ptr(), Lump);
}

//-------------------------------------------------------------------------
// Purpose       : Get parent Body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/11/03
//-------------------------------------------------------------------------
Body* RefVolume::get_body_ptr()
{
  SenseEntity* co_volume_ptr = get_first_sense_entity_ptr();
  
  if (!co_volume_ptr                  // no covolume 
    || co_volume_ptr->next_on_bte() ) // multiple covolumes )
  {
    PRINT_ERROR("%s:%d RefVolume %d has other than one parent CoVolume.\n"
                "This is a BUG -- please report it.\n",
                __FILE__, __LINE__, id() );
    return 0;
  }
  
  return dynamic_cast<Body*>(co_volume_ptr->get_grouping_entity_ptr());
}

int RefVolume::genus() 
{
  int i;
  DLIList<RefFace*> faces;
  ref_faces(faces);
  int gs = 0;
  for (i = faces.size(); i > 0; i--)
    gs += faces.get_and_step()->genus();
  
  return 1 - (num_ref_vertices() - num_ref_edges() + faces.size() - gs)/2;
}

CubitVector RefVolume::center_point()
{
  return bounding_box().center();
}

int RefVolume::dimension() const
{
  return 3;
}

CubitStatus RefVolume::mass_properties( CubitVector &centroid, double &volume )
{
  DLIList<Body*> bodies;
  this->bodies( bodies );
  if( bodies.get()->is_sheet_body() )
  {
    centroid.set(0,0,0);
    volume = 0;
    return CUBIT_SUCCESS;
  }
  else
  {
    Lump *lump = get_lump_ptr();
    return lump->mass_properties( centroid, volume );
  }
}

CubitString RefVolume::measure_label()
{
  return "volume";
}

// Draws each surface in the volume
/*
void RefVolume::draw_my_faces ( int color )
{
   DLIList<RefFace*> ref_face_list;
   int i;

   if (color == CUBIT_DEFAULT_COLOR)
      color = this->color();
   ref_faces ( ref_face_list );
   ref_face_list.reset();
   for (i = 0; i < ref_face_list.size(); i++ )
      ref_face_list.get_and_step()->draw(color);
}
*/

//-------------------------------------------------------------------------
// Purpose       : Initializes member data
//
// Special Notes :
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/07/96
//-------------------------------------------------------------------------
void RefVolume::initialize()
{
   // Set the Entity and Block IDs for this new RefVolume
   GeometryEntity* geom_ptr = get_geometry_entity_ptr();
   int saved_id = geom_ptr->get_saved_id();
   if ( !saved_id || RefEntityFactory::instance()->get_ref_volume(saved_id) )
   {
     saved_id =  RefEntityFactory::instance()->next_ref_volume_id();
     geom_ptr->set_saved_id(saved_id);
   }
   entityId = saved_id;

     // read and initialize attributes
   auto_read_cubit_attrib();
   auto_actuate_cubit_attrib();

     // Assign a default entity name
   assign_default_name();   
}

// unmark this face, and
// recursively unmark all other marked faces sharing an edge with this face.
void unmark_component( RefFace *face )
{
  assert( face->marked() );
  face->marked( CUBIT_FALSE );
  DLIList<RefEdge*> edges;
  face->ref_edges( edges );
  DLIList<RefFace*> edge_faces;
  RefEdge *edge;
  for ( int e = edges.size(); e--; ) 
  {
    edge = edges.get_and_step();
    edge_faces.clean_out();
    edge->ref_faces( edge_faces );
    for ( int f = edge_faces.size(); f--; ) 
    {
      face = edge_faces.get_and_step();
      if ( face->marked() )
         unmark_component( face );
    }
  }
}

int RefVolume::num_boundary_components() 
{
  // get all the surfaces
  // march from one to the next, using shared-ref_edge connectivity
  
  int num_components = 0;
  DLIList<RefFace*> face_list;
  ref_faces( face_list );

  // mark them as belonging to this volume and not seen yet.
  int i;
  RefFace *face;
  for( i = face_list.size(); i--; )
  {
    face = face_list.get_and_step();
    face->marked( CUBIT_FALSE );
  }
  
  for ( i = face_list.size(); i--; ) 
  {
    face = face_list.get_and_step();
    assert( !face->marked() );
    face->marked( CUBIT_TRUE );
  }

  for ( i = face_list.size(); i--; ) 
  {
    face = face_list.get_and_step();
    if ( face->marked() ) 
    {
      num_components++;
      unmark_component( face );      
    }
  }
  return num_components;
}

CubitBoolean RefVolume::about_spatially_equal( 
    RefVolume* ref_vol_ptr_2,double tolerance_factor)
{
     // Get rid of the trivial case...
   if( this == ref_vol_ptr_2)
   {
      return CUBIT_TRUE;
   }

   DLIList<RefEdge*> ref_edge_list_1, ref_edge_list_2;
   this->ref_edges( ref_edge_list_1 );
   ref_vol_ptr_2->ref_edges( ref_edge_list_2 );
      
        //compare the size of the two lists.
   if ( ref_edge_list_1.size() != ref_edge_list_2.size() )
     return CUBIT_FALSE;
      
   DLIList<RefFace*> ref_face_list_1, ref_face_list_2;
   this->ref_faces( ref_face_list_1 );
   ref_vol_ptr_2->ref_faces( ref_face_list_2 );
      
        //compare the size of the two lists.
   if ( ref_face_list_1.size() != ref_face_list_2.size() )
     return CUBIT_FALSE;
      
   
     //This compare precedure does the following :
     //   1. Test the bounding boxes of the 2 volumes for equality;
        
   CubitBox box_1 = this->bounding_box();
   CubitBox box_2 = ref_vol_ptr_2->bounding_box();

   CubitBoolean passes_b_box = CUBIT_FALSE;

       // This test checks to see that the min and max vectors of the
       // bounding boxes are within 10% of the length of the bbox diagonal.
       // Note that this assumes the default values of resabs=1e-6 and
       // tolerance_factor=500

   CubitVector tol_vect(
      CUBIT_MIN(box_1.x_range(), box_2.x_range()),
      CUBIT_MIN(box_1.y_range(), box_2.y_range()),
      CUBIT_MIN(box_1.z_range(), box_2.z_range()) );
   tol_vect *= 200.0 * tolerance_factor * GEOMETRY_RESABS;
   if( tol_vect.x() < GEOMETRY_RESABS )
      tol_vect.x(GEOMETRY_RESABS);
   if( tol_vect.y() < GEOMETRY_RESABS )
      tol_vect.y(GEOMETRY_RESABS);
   if( tol_vect.z() < GEOMETRY_RESABS )
      tol_vect.z(GEOMETRY_RESABS);
   
   if((fabs(box_1.minimum().x() - box_2.minimum().x()) < tol_vect.x()) &&
      (fabs(box_1.maximum().x() - box_2.maximum().x()) < tol_vect.x()) &&
      (fabs(box_1.minimum().y() - box_2.minimum().y()) < tol_vect.y()) &&
      (fabs(box_1.maximum().y() - box_2.maximum().y()) < tol_vect.y()) &&
      (fabs(box_1.minimum().z() - box_2.minimum().z()) < tol_vect.z()) &&
      (fabs(box_1.maximum().z() - box_2.maximum().z()) < tol_vect.z())) 
   {
     passes_b_box = CUBIT_TRUE;
   }
   
   return passes_b_box;
}    

//-------------------------------------------------------------------------
// Purpose       : Check if all shells are sheets
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/15/04
//-------------------------------------------------------------------------
CubitBoolean RefVolume::is_sheet()
{
  DLIList<Shell*> shells;
  this->shells(shells);
  while (shells.size())
    if (!shells.pop()->is_sheet())
      return CUBIT_FALSE;
  return CUBIT_TRUE;
}
  
int RefVolume::validate()
{
     //- This function determines whether the entity is valid.
     //- Several types of checks can be done, 
   int error = 0;

     // Perform general RefEntity checks (measure > 0)
   error += RefEntity::validate();

     // Pass through to surface and add in its validation
   Lump *lump = this->get_lump_ptr();

     // check surface ptr
   if (lump != NULL) {
        // Check underlying surface
	   DLIList <TopologyEntity*> bad_entities;
      error += lump->validate(entity_name(), bad_entities);
   } else {
      PRINT_WARNING("\tWARNING: Null underlying volume for %s, (%s %d)\n",
                    entity_name().c_str(), class_name(), id());
      error++;
   }
   return error;
}

