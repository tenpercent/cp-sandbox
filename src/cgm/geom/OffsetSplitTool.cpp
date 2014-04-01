//-------------------------------------------------------------------------
// Filename      : OffsetSplitTool.cpp
//
// Purpose       :
//
//   Split Surface <id_list> Offset Curve <id_list>
//       Distance <val> [Segment <val>] [Partition] [Blunt] [Preview [Create]]
//
// Special Notes :
//
// Creator       : Sam Showman
//
// Creation Date : 05/10/2005
//-------------------------------------------------------------------------
#include "OffsetSplitTool.hpp"
#include "AnalyticGeometryTool.hpp"
#include "RefFace.hpp"
#include "Point.hpp"
#include "Curve.hpp"
#include "CubitMessage.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryModifyEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "CubitUtil.hpp"
#include "DLIList.hpp"
#include "GfxPreview.hpp"
#include "GMem.hpp"
#include "Curve.hpp"
#include "BodySM.hpp"
#include "TopologyBridge.hpp"
#include "Surface.hpp"
#include "CubitBox.hpp"

// #define OFFSET_SPLIT_DEBUG

double OffsetSplitTool::tolerance = .001;

OffsetSplitTool::OffsetSplitTool()
{
}

CubitStatus OffsetSplitTool::split_surfaces_offset(DLIList<RefFace*> &ref_face_list,
                                                   DLIList<RefEdge*> &edge_list,
                                                   int num_segs,
                                                   double distance_in,
                                                   CubitBoolean divide_flg,
                                                   CubitBoolean blunt_flg,
                                                   CubitBoolean preview_flg,
                                                   CubitBoolean create_ref_edges_flg)
{
    ref_face_list.uniquify_unordered();
    edge_list.uniquify_unordered();

    // Check for valid number of segments
    if( num_segs < 1 )
    {
        PRINT_ERROR( "Number of specified segments must be >= 1\n" );
        return CUBIT_FAILURE;
    }

    if(ref_face_list.size() < 1)
    {
        PRINT_ERROR( "No surfaces specified for splitting\n" );
        return CUBIT_FAILURE;
    }

    if(edge_list.size() < 1)
    {
        PRINT_ERROR( "No curves specified for splitting\n" );
        return CUBIT_FAILURE;
    }

    DLIList<Curve*> offset_curves;
    DLIList<DLIList<Curve*>*> imprint_list;
    for(int seg = 1;seg<num_segs+1;seg++)
    {
        // set the current distance
        double distance = (distance_in/(double)num_segs)*(seg);

        // generate a swept body at each curve
        DLIList<BodySM*> body_list;
        DLIList<Point*> points;
        GeometryModifyEngine* gme = 0;
        GeometryQueryEngine* gqe = 0;
        edge_list.reset();
        int i = 0;
        for(i =edge_list.size();i--;)
        {
            RefEdge* ref_edge = edge_list.get_and_step();
            Curve* cur_curve = ref_edge->get_curve_ptr();

            // add the curve to the class list of source curves
            // will use this latter for dividing surfaces
            sourceCurves.append(cur_curve);

            gqe = cur_curve->get_geometry_query_engine();
            gme = GeometryModifyTool::instance()->get_engine(cur_curve);
            tolerance = gqe->get_sme_resabs_tolerance();

            if(!gme)
            {
                PRINT_ERROR( "No ModifyEngine found!\n" );
                return CUBIT_FAILURE;
            }

            cur_curve->points(points);
            Curve* swept_curve = create_sweep_curve(cur_curve,distance,tolerance,CUBIT_TRUE);
            if(!swept_curve)
            {
                PRINT_WARNING("Failed to create offset geometry for curve %d.\n",ref_edge->id());
                continue;
            }

            CubitVector swpt_loc, source_loc;
            cur_curve->position_from_fraction(0.0,source_loc);
            swept_curve->position_from_fraction(0.0,swpt_loc);

            CubitVector up_vector = swpt_loc-source_loc;

            Surface* section = create_sweep_section(cur_curve,distance,up_vector);
            if(!section)
            {
                PRINT_WARNING("Failed to create offset geometry for curve %d.\n",ref_edge->id());
                // delete the Curve path
                gqe->delete_solid_model_entities(swept_curve);
                continue;
            }

            BodySM* swept_body = create_sweep_body(section,swept_curve);

            // delete the Curve path and Surface section
            gqe->delete_solid_model_entities(swept_curve);
            if(!swept_body)
            {
                PRINT_WARNING("Failed to create offset geometry for curve %d.\n",ref_edge->id());
                continue;
            }

            body_list.append(swept_body);
        }

        if(body_list.size() == 0)
        {
            PRINT_ERROR("Failed to offset any curve.\n");
            return CUBIT_FAILURE;
        }

        DLIList<Point*> blunt_points;
        if(blunt_flg == CUBIT_TRUE)
        {
            for(i = 0; i < points.size(); i++)
            {
                bool remove_point = true;
                CubitVector pnt_loc_0 = points[i]->coordinates();
                int j = 0;
                for(j = 0;j < points.size(); j++)
                {
                    if(i==j)
                        continue;

                    if(points[j] == points[i])
                    {
                        remove_point = false;
                        break;
                    }

                    CubitVector pnt_loc_1 = points[j]->coordinates();

                    if(pnt_loc_1.within_tolerance(pnt_loc_0,tolerance))
                    {
                        remove_point = false;
                        break;
                    }
                }
                if(remove_point)
                    blunt_points.append(points[i]);
            }
        }
        points.uniquify_unordered();
        points -= blunt_points;

        for(i = points.size(); i--;)
        {
            BodySM* shpere_body = gme->sphere(distance);
            CubitVector trans_dist = points.get_and_step()->coordinates();
            gqe->translate(shpere_body,trans_dist);
            body_list.append(shpere_body);
        }

        // create the offset geometry for debuging
#ifdef OFFSET_SPLIT_DEBUG
        for(i=body_list.size();i--;)
        {
            GeometryQueryTool::instance()->make_Body(
                gme->copy_body(body_list.get_and_step()));
        }
#endif

        DLIList<BodySM*> united_body_list;
        if(body_list.size() == 1)
            united_body_list = body_list;
        else
        {
#ifdef BOYD17 
            BodySM* return_body = 0;
#endif
            if(!gme->unite(body_list,united_body_list))//united_body_list))
            {
                if(body_list.size() == 1)
                    PRINT_ERROR( "Offset split failed at the unite step.\n" );
                else
                    PRINT_ERROR( "Offset split failed at the unite step. Try reducing the curve count\n" );

                // delete the solid bodies
                for(i = body_list.size(); i--;)
                    gqe->delete_solid_model_entities(body_list.get_and_step());

                return CUBIT_FAILURE;
            }
        }

        // get the surfaces contained in the body
        DLIList<Surface*> surface_united_list;
        for(i = united_body_list.size(); i--;)
        {
            BodySM* cur_body = united_body_list.get_and_step();
            cur_body->surfaces(surface_united_list);
        }

        ref_face_list.reset();
        imprint_list.reset();
        for(i = ref_face_list.size(); i--;)
        {
            Surface* surf_1 = ref_face_list.get_and_step()->get_surface_ptr();
            DLIList<Curve*> cur_curves;
            for(int j = surface_united_list.size(); j--;)
            {
                DLIList<Curve*> intersection_curves;
                Surface* surf_2 = surface_united_list.get_and_step();
                gme->surface_intersection(surf_1,surf_2,intersection_curves,0.0);

                // remove any very small curves from the intersection graph
                // for some reason ACIS seems to return zero length curves
                for(int k=0;k<intersection_curves.size();k++)
                {
                    Curve* cur_curve = intersection_curves.get_and_step();
                    double len = cur_curve->
                        length_from_u(cur_curve->start_param(),cur_curve->end_param());

                    if(fabs(len) < tolerance)
                    {
                        intersection_curves.remove(cur_curve);
                        gqe->delete_solid_model_entities(cur_curve);
                    }
                }

                // do the division if the user asks for it
                if(divide_flg &&
                    num_segs == seg &&
                    // skip the curves created by sphere surfaces
                    surf_2->geometry_type() != SPHERE_SURFACE_TYPE &&
                    surf_2->geometry_type() != PLANE_SURFACE_TYPE)
                {
                    intersection_curves +=
                        this->create_divide_curves(intersection_curves,surf_1,distance);
                }

                cur_curves += intersection_curves;
            }
            DLIList<Curve*>* surf_inter_list = 0;

            if(ref_face_list.size() > imprint_list.size())
            {
                surf_inter_list = new DLIList<Curve*>;
                imprint_list.append(surf_inter_list);
            }
            else
                surf_inter_list = imprint_list.get_and_step();

            (*surf_inter_list) += cur_curves;
            offset_curves += cur_curves;

            // draw the preview as the offset curves are being calculated
            if(preview_flg && !create_ref_edges_flg)
                draw_preview(cur_curves);
        }

        // delete the solid bodies
        for(i = united_body_list.size(); i--;)
            gqe->delete_solid_model_entities(united_body_list.get_and_step());
    }

    // ok, now do the imprint or preview
    if(preview_flg)
    {
        if(!create_ref_edges_flg)
        {
            //Delete the Curve entities
            for(int i=0;i<offset_curves.size();i++)
            {
                Curve* cur_curve = offset_curves.get_and_step();
                cur_curve->get_geometry_query_engine()->
                    delete_solid_model_entities(cur_curve);
            }
        }
        else
        {
            create_ref_edges(offset_curves);
        }

        // delete the list that contained the curves that would have been imprinted
        for(int i= imprint_list.size();i--;)
            delete imprint_list.get_and_step();
    }
    // Imprint
    else
    {
        int i = 0;
        ref_face_list.reset();
        DLIList<Surface*> surface_list;
        for( i = ref_face_list.size(); i--;)
            surface_list.append(ref_face_list.get_and_step()->get_surface_ptr());

        Body *new_body_ptr;
        if( GeometryModifyTool::instance()->imprint( surface_list,
            imprint_list, new_body_ptr ) == CUBIT_FAILURE )
        {
            //Delete the Curve entities
            for(i= offset_curves.size();i--;)
            {
                Curve* cur_curve = offset_curves.get_and_step();
                cur_curve->get_geometry_query_engine()->
                    delete_solid_model_entities(cur_curve);
            }

            // delete the curves that would have been imprinted
            for( i= imprint_list.size();i--;)
                delete imprint_list.get_and_step();

            return CUBIT_FAILURE;
        }

        //Delete the Curve entities
        for( i= offset_curves.size();i--;)
        {
            Curve* cur_curve = offset_curves.get_and_step();
            cur_curve->get_geometry_query_engine()->
                delete_solid_model_entities(cur_curve);
        }

        // delete the curves that were used for imprinted
        for( i= imprint_list.size();i--;)
            delete imprint_list.get_and_step();
    }

    return CUBIT_SUCCESS;
}

//- Begin Private Functions
CubitStatus OffsetSplitTool::draw_preview(
    DLIList<Curve*> &curve_list,
    int color )
{
    // clear any previous previews
    GfxPreview::clear();

    int i;
    Curve *curve_ptr;
    curve_list.reset();
    for( i=curve_list.size(); i--; )
    {
        curve_ptr = curve_list.get_and_step();
        draw_preview( curve_ptr, CUBIT_FALSE, color );
    }

    GfxPreview::flush();

    return CUBIT_SUCCESS;
}

CubitStatus OffsetSplitTool::draw_preview(
    Curve *curve_ptr,
    CubitBoolean flush,
    int color )
{
    int num_points;
    CubitStatus result;
    GMem g_mem;

    // clear any previous previews
    GfxPreview::clear();

    // get the graphics
    result = curve_ptr->get_geometry_query_engine()->
        get_graphics( curve_ptr, num_points, &g_mem );

    if (result==CUBIT_FAILURE || num_points == 0)
    {
        PRINT_WARNING("Unable to preview a curve\n" );
    }

    // Draw the polyline
    GfxPreview::draw_polyline( g_mem.point_list(), g_mem.pointListCount, color );
    if( flush )
        GfxPreview::flush();

    return CUBIT_SUCCESS;
}


Curve *OffsetSplitTool::create_sweep_curve(
    Curve *curve_in,
    double distance,
    double chord_tol,
    CubitBoolean iterate)
{
    GeometryModifyEngine* gme =
        GeometryModifyTool::instance()->get_engine( curve_in );
    GeometryQueryEngine* gqe = curve_in->get_geometry_query_engine();

    CubitVector start_pos;
    CubitVector end_pos;
    CubitVector dummy_pos;
    CubitVector x_axis;
    CubitVector y_axis;
    CubitVector z_axis;
    CubitVector curvature;

    if(!curve_in->position_from_fraction(0.0,start_pos))
        return 0;

    if(!curve_in->position_from_fraction(1.0,end_pos))
        return 0;

    curve_in->closest_point(start_pos,dummy_pos,&z_axis,&curvature);

    // check to see if the curve is planar
    CubitBoolean is_planar = CUBIT_FALSE;
    DLIList<Surface*> attached_surfs;
    curve_in->surfaces(attached_surfs);
    int i = 0;
    for(i = 0;i < attached_surfs.size();i++)
    {
        Surface* cur_surf = attached_surfs.get_and_step();
        if(cur_surf && 
            cur_surf->geometry_type() == PLANE_SURFACE_TYPE)
        {
            CubitVector center_vec;
            cur_surf->get_point_normal(center_vec,x_axis);
            is_planar = CUBIT_TRUE;
            break;
        }
    }

    if(is_planar)
    {
        y_axis = z_axis*x_axis;
    }
    // If curvature is significant, use it to create the sweep coordinate system
    else if(curvature.length() < CUBIT_DBL_MIN)
    {
        z_axis.orthogonal_vectors(x_axis,y_axis);
    }
    else
    {
        y_axis = curvature;
        x_axis = z_axis*y_axis;
    }

    x_axis.normalize();
    y_axis.normalize();
    z_axis.normalize();
    x_axis *=distance;

    GeometryType geom_type = curve_in->geometry_type();
    if(geom_type == STRAIGHT_CURVE_TYPE) // if the curve is a line just offset it
    {
        // ACIS has problems calculating the intersection curves of surface edges and surfaces
        // so rotate the coordinate system to some angle that is less likely to be used
        x_axis.normalize();
        y_axis.normalize();
        CubitVector move_vec = y_axis - x_axis;
        y_axis = y_axis+(move_vec*.744561789854);
        x_axis = z_axis*y_axis;

        x_axis.normalize();
        y_axis.normalize();
        z_axis.normalize();
        x_axis *=distance;

        // just offset the curve
        Curve* offset_curve = gme->make_Curve(curve_in);
        gqe->translate(offset_curve,x_axis);
        return offset_curve;
    }
    else if(geom_type == ARC_CURVE_TYPE ||
        geom_type == ELLIPSE_CURVE_TYPE ||
        is_planar) // if the curve is a arc just offset it
    {
        // just offset the curve
        Curve* offset_curve = gme->make_Curve(curve_in);
        gqe->translate(offset_curve,x_axis);
        return offset_curve;
    }

#ifdef BOYD17 
    int frac_div = 1;
#endif
    DLIList<CubitVector*> pnt_list;
    CubitVector* loc_0 = new CubitVector(start_pos);
    CubitVector* loc_1 = new CubitVector(end_pos);
    pnt_list.append(loc_0);
    pnt_list.append(loc_1);
    DLIList<double> frac_list;
    frac_list.append(0.0);
    frac_list.append(1.0);
    frac_list.reset();

    // create the curve points with a maximum deviation of 'chord_tol'
    while(1)
    {
        double frac_0 = frac_list.get_and_step();
        double frac_1 = frac_list.get();
        pnt_list.step();

#ifdef BOYD17 
        if(frac_0>frac_1)
            double test_0 =  frac_1 -frac_0;
#endif

        if(frac_0>frac_1)
            break;

        CubitVector tmp_pos_0;
        curve_in->position_from_fraction(frac_0,tmp_pos_0);
        CubitVector tmp_pos_1;
        curve_in->position_from_fraction(frac_1,tmp_pos_1);
        CubitVector tmp_pos_2 = (tmp_pos_0+tmp_pos_1)/2.0;
        CubitVector tmp_pos_3;
        double mid_frac = (frac_0+frac_1)/2.0;
        curve_in->position_from_fraction(mid_frac,tmp_pos_3);

        if(tmp_pos_3.distance_between(tmp_pos_2) > chord_tol)
        {
            CubitVector* loc = new CubitVector(tmp_pos_3);
            pnt_list.back();
            frac_list.back();
            pnt_list.insert(loc);
            frac_list.insert(mid_frac);
            pnt_list.back();
            frac_list.back();
        }
    }

    // generate the offset curve points
    CubitVector last_vec = x_axis;
    pnt_list.reset();
    for(i = 0;i < pnt_list.size();i++)
    {
        CubitVector curv,tangent,dummy;
        CubitVector *cur_pnt = pnt_list.get_and_step();
        curve_in->closest_point(*cur_pnt,dummy,&tangent,&curv);

        // if the radius of curvature is smaller than the distance of the
        // sweep offset then return a NULL curve
        if(curvature.length()> CUBIT_DBL_MIN &&
            1.0/(curv.length()) <= distance)
        {
            // clean out the points
            for(int i =pnt_list.size();i--;)
                delete pnt_list.get_and_step();

            PRINT_WARNING("Radius of curvature is smaller than or equal to %f.\n"
                "Cannot offset curve.\n",distance);
            return 0;
        }

        CubitVector up_vec = tangent*last_vec;
        CubitVector offset_vec = up_vec*tangent;
        offset_vec.normalize();
        offset_vec *= distance;
        cur_pnt->set(*cur_pnt + offset_vec);
        last_vec = offset_vec;
    }

    // create a spline sweep path
    pnt_list.reset();
    Point* start_point = gme->make_Point(*(pnt_list.get_and_back()));
    Point* end_point = gme->make_Point(*(pnt_list.get_and_step()));

    // Create the curve

    Curve* curve_out =
        gme->make_Curve(
        SPLINE_CURVE_TYPE,
        start_point,end_point,pnt_list);

    pnt_list.reset();
    for(i =pnt_list.size();i--;)
        delete pnt_list.get_and_step();

    pnt_list.clean_out();
    return curve_out;
}


Surface* OffsetSplitTool::create_sweep_section(
    Curve* path,
    double distance,
    CubitVector up_vector,
    double fraction)
{
    up_vector.normalize();
    GeometryModifyEngine* gme =
        GeometryModifyTool::instance()->get_engine( path );

    CubitVector start_pos;
    CubitVector dummy_pos;
    CubitVector x_axis;
    CubitVector y_axis;
    CubitVector z_axis;
    if(!path->position_from_fraction(fraction,start_pos))
        return 0;

    path->closest_point(start_pos,dummy_pos,&z_axis);
    y_axis = z_axis*up_vector;
    x_axis = up_vector;

    x_axis.normalize();
    y_axis.normalize();
    z_axis.normalize();

    Point* start_point = gme->make_Point(start_pos+x_axis);
    Point* end_point = gme->make_Point(start_pos-x_axis);
    Point* center_point = gme->make_Point(start_pos);

    Curve* arc_curve_0 =
        gme->create_arc_center_edge(center_point,start_point,end_point,z_axis,distance);

    start_point = gme->make_Point(start_pos+x_axis);
    end_point = gme->make_Point(start_pos-x_axis);
    center_point = gme->make_Point(start_pos);

    Curve* arc_curve_1 =
        gme->create_arc_center_edge(center_point,start_point,end_point,-z_axis,distance);

    DLIList<Curve*> curve_list;
    curve_list.append(arc_curve_0);
    curve_list.append(arc_curve_1);

    Surface* surf_out =
        gme->make_Surface(PLANE_SURFACE_TYPE,curve_list,0,false);

    return surf_out;
}

BodySM* OffsetSplitTool::create_sweep_body(
    Surface* section,
    Curve* path)
{
    GeometryModifyEngine* gme =
        GeometryModifyTool::instance()->get_engine( section );

    DLIList<GeometryEntity*> ent_list;
    ent_list.append(section);
    DLIList<BodySM*> body_list;
    DLIList<Curve*> path_list;
    path_list.append(path);

    gme->sweep_along_curve(ent_list,body_list,path_list);

    if(body_list.size())
        return body_list.get();

    return 0;
}

CubitStatus OffsetSplitTool::create_ref_edges(
    DLIList<Curve*> &curve_list )
{
    int i;
    Curve *curve_ptr;
    curve_list.reset();
    for( i=curve_list.size(); i--; )
    {
        curve_ptr = curve_list.get_and_step();
        GeometryQueryTool::instance()->make_free_RefEdge(curve_ptr);
    }
    return CUBIT_SUCCESS;
}

DLIList<Curve*> OffsetSplitTool::create_divide_curves(
    DLIList<Curve*> &curve_list,
    Surface* target_surface,
    double distance)
{
    DLIList<Curve*> return_list;
    curve_list.reset();
    int curve_list_size = curve_list.size();
    for(int i = curve_list_size; i--;)
    {
        Curve* cur_curve = curve_list.get_and_step();
        GeometryModifyEngine* gme =
            GeometryModifyTool::instance()->get_engine( cur_curve );
        GeometryQueryEngine* gqe =
            cur_curve->get_geometry_query_engine();

        // Find the matching source curve.
        // there are better ways to do this such as using the surface
        // arc intersection but that only seems to work for one engine or another
        // not both
        for(int j = sourceCurves.size(); j--;)
        {
            Curve* source_curve = sourceCurves.get_and_step();
            CubitVector coord1;
            CubitVector coord2;
            double distance_out = 0.0;
            if(!gqe->entity_entity_distance( source_curve, cur_curve,
                coord1, coord2, distance_out ))
                continue;

            if(fabs(distance_out-distance) < tolerance)
            {
                // check three points
                CubitVector pnt_0;
                CubitVector pnt_1;
                CubitVector pnt_2;
                CubitVector pnt_3;
                CubitVector pnt_4;
                CubitVector pnt_5;
                cur_curve->position_from_fraction(0.0,pnt_0);
                source_curve->closest_point(pnt_0,pnt_1);
                cur_curve->position_from_fraction(0.5,pnt_2);
                source_curve->closest_point(pnt_2,pnt_3);
                cur_curve->position_from_fraction(1.0,pnt_4);
                source_curve->closest_point(pnt_4,pnt_5);

                if(fabs(pnt_0.distance_between(pnt_1) - distance) < tolerance &&
                    fabs(pnt_2.distance_between(pnt_3) - distance) < tolerance &&
                    fabs(pnt_4.distance_between(pnt_5) - distance) < tolerance)
                {
                    // ok now just create curves between each end of the surface intersection curves
                    // and the projection to the source curves

                    // first curve
                    Point* start_point_0 = gme->make_Point(pnt_0);
                    Point* end_point_0 = gme->make_Point(pnt_1);
                    CubitVector vert3_coord_0 = start_point_0->coordinates();

                    Curve* new_curve_0 =
                        gme->make_Curve(STRAIGHT_CURVE_TYPE,
                        start_point_0,end_point_0,&vert3_coord_0,CUBIT_FORWARD);

                    // second curve
                    Point* start_point_1 = gme->make_Point(pnt_4);
                    Point* end_point_1 = gme->make_Point(pnt_5);
                    CubitVector vert3_coord_1 = start_point_1->coordinates();

                    Curve* new_curve_1 =
                        gme->make_Curve(STRAIGHT_CURVE_TYPE,
                        start_point_1,end_point_1,&vert3_coord_1,CUBIT_FORWARD);

                    // if the target surface is not a plane then do a projection
                    if(target_surface->geometry_type() != PLANE_SURFACE_TYPE)
                    {
                        DLIList<Curve*> curves_in;

                        if(new_curve_0)
                            curves_in.append(new_curve_0);

                        if(new_curve_1)
                            curves_in.append(new_curve_1);

                        DLIList<Curve*> curves_out;
                        DLIList<Surface*> surface_in;
                        surface_in.append(target_surface);
                        if(gme->project_edges(surface_in,curves_in,curves_out))
                            return_list += curves_out; // add the curves to the return list
                    }
                    else // add the curves to the return list
                    {
                        if(new_curve_0)
                            return_list.append(new_curve_0);

                        if(new_curve_1)
                            return_list.append(new_curve_1);
                    }
                }
            }
        }
    }
    return return_list;
}
