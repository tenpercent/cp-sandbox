//- Class: SurfParamTool
//-------------------------------------------------------------------------
// Filename      : SurfParamTool.cpp
//
// Purpose       : This is the generic version of ParamTool when the
//                 geometry engine has a sufficient parameterization.
//                 It uses the tool's existing functions to get transform
//                 between uv and xyz spaces.
//
// Creator       : Christopher Hynes
//
// Creation Date : 7/10/2002
//
// Owner         : Christopher Hynes
//-------------------------------------------------------------------------

#include "SurfParamTool.hpp"
#include "CastTo.hpp"
#include "Surface.hpp"
#include "DLIList.hpp"

//-------------------------------------------------------------------------
// Function:    SurfParamTool
// Description: constructor
// Author:      chynes
// Date:        7/10/2002
//-------------------------------------------------------------------------
SurfParamTool::SurfParamTool(Surface *surf) 
{
//	FacetEvalTool **surf_eval_tool = NULL;
//	DLIList<CubitFacet *> facet_list;

	//- update private variables
	refSurf = surf;
//	refSurf->setup_use_facets(facet_list,0,surf_eval_tool);
//	feTool = *surf_eval_tool;

}

//-------------------------------------------------------------------------
// Function:    SurfParamTool
// Description: deconstructor
// Author:      chynes
// Date:        7/10/2002
//-------------------------------------------------------------------------
SurfParamTool::~SurfParamTool() {}

//===================================================================================
// Function: set_up_space (Public)
// Description: sets up space of flattening
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus SurfParamTool::set_up_space(void) {
	
	//double u, v, s;
	CubitStatus rv = CUBIT_SUCCESS;

	//- calculate uv position for every point
//	TDFacetBoundaryPoint *td_fbp;
//	CubitBoolean on_internal_boundary;
//	DLIList<CubitFacet *> facet_list;
//	DLIList<CubitPoint *> point_list;
//	CubitPoint *pt_ptr;
//	int ii;
//	feTool->get_points(point_list);
//	for(ii = 0; ii < point_list.size() && rv == CUBIT_SUCCESS; ii++) {
//		pt_ptr = point_list.get_and_step();
//		rv = refSurf->u_v_from_position(pt_ptr->coordinates(), u, v);
//		facet_list.clean_out();
//		pt_ptr->facets_on_surf( feTool->tool_id(), facet_list, on_internal_boundary );
//		if (on_internal_boundary)
//		{
//			td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( pt_ptr );
//			if (!td_fbp)
//			{
//				TDFacetBoundaryPoint::add_facet_boundary_point( pt_ptr );
//				td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( pt_ptr );
//				td_fbp->add_surf_facets( facet_list );
//				td_fbp->set_uv( facet_list.get(), u, v);
//			}
//			else
//			{
//				if (td_fbp->set_uv( facet_list.get(),
//					u, v ) != CUBIT_SUCCESS)
//				{
//					td_fbp->add_surf_facets( facet_list );
//					td_fbp->set_uv( facet_list.get(), u, v);
//				}
//			}
//		}
//		else
//		{
//			pt_ptr->set_uv( u, v);
//		}
//	}
//
//	//- calculate distortion [sum(original facet area/uv area)/total facets]
//	CubitFacet *facet_ptr;
//	CubitPoint *p0, *p1, *p2;
//	int jj;
//	double orig_area, uv_area, sum_area, u0, u1, u2, v0, v1, v2;
//	for(ii = 0; ii < point_list.size() && rv == CUBIT_SUCCESS; ii++) {
//		pt_ptr = point_list.get_and_step();
//		facet_list.clean_out();
//		pt_ptr->facets_on_surf(feTool->tool_id(), facet_list, on_internal_boundary);
//		
//		sum_area = 0;
//		for(jj = 0; jj < facet_list.size(); jj++){
//			facet_ptr = facet_list.get_and_step();
//
//			//- calculate original area
//			orig_area = facet_ptr->area();
//
//			//- extract uv data
//			facet_ptr->points(p0, p1, p2);
//			p0->get_uv(facet_ptr, u0, v0);
//			p1->get_uv(facet_ptr, u1, v1);
//			p2->get_uv(facet_ptr, u2, v2);
//
//			//- calculate uv area
//			CubitVector uv0(u0,v0,0.0);
//			CubitVector uv1(u1,v1,0.0);
//			CubitVector uv2(u2,v2,0.0);
//			CubitVector e0(uv0, uv1);
//			CubitVector e1(uv0, uv2);
//			uv_area = (e0*e1).length()*0.5;
//
//			sum_area += sqrt(uv_area/orig_area);
//		}
//	
//		s = sum_area/facet_list.size();
//		pt_ptr->get_uv(facet_list.get(), u, v);
//
//		if (on_internal_boundary)
//		{
//			td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( pt_ptr );
//			if (!td_fbp)
//			{
//				TDFacetBoundaryPoint::add_facet_boundary_point( pt_ptr );
//				td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( pt_ptr );
//				td_fbp->add_surf_facets( facet_list );
//				td_fbp->set_uvs( facet_list.get(), u, v, s);
//			}
//			else
//			{
//				if (td_fbp->set_uvs( facet_list.get(),
//					u, v, s ) != CUBIT_SUCCESS)
//				{
//						td_fbp->add_surf_facets( facet_list );
//						td_fbp->set_uvs( facet_list.get(), u, v, s);
//				}
//			}
//		}
//		else
//		{
//			pt_ptr->set_uvs( u, v, s);
//		}
//				
//	}
	
	return rv; 
}

//===================================================================================
// Function: transform_to_uv (Public)
// Description: same as title, the local sizing will be returned in the z coord 
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus SurfParamTool::transform_to_uv(CubitVector &xyz_location, CubitVector &uv_location) 
{
	double u,v;
	CubitStatus rv = CUBIT_SUCCESS;
//	DLIList<CubitFacet *> facet_list; 
//	CubitFacet *tri_ptr;
//	CubitBoolean *outside = CUBIT_FALSE;
//	CubitVector closest_point;
//	CubitVector area_coord;
//	double compare_tol;
//
//	// find best compare_tol
//	compare_tol = 1e-3*(feTool->bounding_box().diagonal().length());
//
//	// locate point
//	feTool->get_facets(facet_list);
//	rv = FacetEvalTool::project_to_facets(facet_list, 
//										  tri_ptr, 
//										  0, 
//										  compare_tol, 
//										  xyz_location, 
//										  CUBIT_FALSE, 
//										  outside, 
//										  &closest_point, 
//										  NULL);
//	// determine barycentric coordinates for in facet
//	if(rv == CUBIT_SUCCESS) 
//	{
//		FacetEvalTool::facet_area_coordinate(tri_ptr, closest_point, area_coord);
//
//		// extract data
//		double u0, u1, u2, v0, v1, v2, s0, s1, s2;
//		CubitPoint *p0, *p1, *p2;
//		tri_ptr->points(p0, p1, p2);
//		p0->get_uvs(tri_ptr, u0, v0, s0);
//		p1->get_uvs(tri_ptr, u1, v1, s1);
//		p2->get_uvs(tri_ptr, u2, v2, s2);
//
//		rv = refSurf->u_v_from_position(xyz_location, u, v);
//		// determine sizing
//		s = (area_coord.x()*s0) + (area_coord.y()*s1) + (area_coord.z()*s2);
//
//		uv_location.set(u,v,s);
//	}

	rv = refSurf->u_v_from_position(xyz_location, u, v);
	uv_location.set(u,v,1.0);

	return rv;
}

//===================================================================================
// Function: transform_to_xyz (Public)
// Description: same as title
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus SurfParamTool::transform_to_xyz(CubitVector &xyz_location, CubitVector &uv_location) 
{
	xyz_location = refSurf->position_from_u_v(uv_location.x(), uv_location.y());

	return CUBIT_SUCCESS;
}



//EOF
