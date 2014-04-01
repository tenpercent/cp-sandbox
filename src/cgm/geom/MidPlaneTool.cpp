//-------------------------------------------------------------------------
// Filename      : MidPlaneTool.hpp
//
// Purpose       : Create a mide surface between the two planar
//                 surfaces passed in by the tool.
//
// Special Notes : 
//
// Creator       : David White
//
// Creation Date : 9/30/00
//-------------------------------------------------------------------------

#include "MidPlaneTool.hpp"
#include "CubitVector.hpp"
#include "CubitPlane.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "Loop.hpp"
#include "Body.hpp"
#include "CubitMessage.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "DLIList.hpp"

CubitStatus MidPlaneTool::sort_surfaces(RefFace* base_face, 
					RefFace* predecessor,
                                        DLIList <RefFace*> &mid_surfaces)
{
  DLIList <Loop*> loop_list;
  DLIList <RefEdge*> loop_edges;
  DLIList <RefFace*> edge_faces;
  RefFace *next_face;
  RefEdge *temp_edge;
  Loop* temp_loop;

  base_face->ordered_loops( loop_list );

  if(loop_list.size() == 1)
  {
    if(surfIndex%2==1)
    {
      mid_surfaces.append(base_face);
    }  
    if ( predecessor == NULL )
    {
      PRINT_ERROR("Bad logic in mid-plane extraction of sorting surfaces.\n");
      return CUBIT_SUCCESS;
    }
    //Now test to see if there are other surfaces attached to this
    //single loop (other than predecessor)
    DLIList <RefFace*> other_surfs;
    CubitStatus tmp_stat = get_other_surfs(base_face, predecessor, other_surfs);
    if ( tmp_stat != CUBIT_SUCCESS )
      return CUBIT_FAILURE;
    int jj;
    for ( jj = other_surfs.size(); jj > 0; jj-- )
    {
      next_face = other_surfs.get_and_step();
      surfIndex++;
      sort_surfaces(next_face,base_face,mid_surfaces);
      surfIndex--;
    }
      
    return CUBIT_SUCCESS;
  }

  if(surfIndex%2==1)
  {
    mid_surfaces.append(base_face);
  }
  loop_list.reset();
  loop_list.remove();
  for ( int ii = loop_list.size(); ii > 0; ii-- )
  {
    temp_loop = loop_list.get_and_step();
    loop_edges.clean_out();
    edge_faces.clean_out();
    temp_loop->ref_edges( loop_edges );
    temp_edge = loop_edges.get();
    temp_edge->ref_faces( edge_faces );
    if( edge_faces.size() > 2 )
    {
      PRINT_ERROR("Bad Logic.\n");
      return CUBIT_FAILURE;
    }
  
    if(edge_faces.size() == 2)
    {
      if(edge_faces.get() == base_face)
      {
        next_face = edge_faces.next();
      } else
      {
        next_face = edge_faces.get();
      }

      surfIndex++;

      sort_surfaces(next_face,base_face,mid_surfaces);

      surfIndex--;

    }
  }  
  return CUBIT_SUCCESS;
}

CubitStatus MidPlaneTool::get_other_surfs(RefFace* base_face, 
					  RefFace* predecessor, 
					  DLIList <RefFace*> &other_surfs)
{
  DLIList <Loop*> loop_list;
  DLIList <RefEdge*> loop_edges;
  DLIList <RefFace*> edge_faces;
  RefFace *temp_face;
  RefEdge *temp_edge;
  Loop* temp_loop;

  base_face->loops( loop_list );

  for ( int ii = loop_list.size(); ii > 0; ii-- )
  {
    temp_loop = loop_list.get_and_step();
    loop_edges.clean_out();
    temp_loop->ref_edges( loop_edges );
    for ( int jj = loop_edges.size(); jj > 0; jj-- )
    {
      temp_edge = loop_edges.get_and_step();
      edge_faces.clean_out();
      temp_edge->ref_faces( edge_faces );
      for ( int kk = edge_faces.size(); kk > 0; kk-- )
      {
	temp_face = edge_faces.get_and_step();
	if ( ( temp_face != base_face ) && ( temp_face != predecessor ) )
	{
	  other_surfs.append_unique( temp_face );
	}
      }
    }
  }
  return CUBIT_SUCCESS;
}      

