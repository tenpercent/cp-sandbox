//-------------------------------------------------------------------------
// Filename      : OffsetSplitTool.hpp
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
#ifndef OFFSETSPLITTOOL_HPP
#define OFFSETSPLITTOOL_HPP

#include "CubitDefines.h"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CoEdge.hpp"

class RefEdge;
class RefFace;
class GeometryModifyEngine;
class Curve;
class Surface;
template <class X> class DLIList;

class OffsetSplitTool
{
public:

    OffsetSplitTool();
    ~OffsetSplitTool(){}

    CubitStatus split_surfaces_offset(DLIList<RefFace*> &ref_face_list,
        DLIList<RefEdge*> &edge_list,
        int num_segs,
        double distance,
        CubitBoolean divide_flg = CUBIT_FALSE,
        CubitBoolean blunt_flg = CUBIT_FALSE,
        CubitBoolean preview_flg = CUBIT_FALSE,
        CubitBoolean create_ref_edges_flg = CUBIT_FALSE);
    //- Splits surface by offsetting curves.
    //- ref_face_list - list of surfaces to split.
    //- num_segs - the number of segments to create (must be >= 1 );
    //- distance - distance of the offset
    //- divide_flg - divide the surface between each offset curve to create a better
    //-              region for map mesh.
    //- blunt_flg - create a blunt type ending instead of the arc type ending
    //- preview_flg - if CUBIT_TRUE, just draw the curves that will be used to split
    //-               instead of actually splitting.
    //- create_ref_edges_flg - valid only if preview_flg=CUBIT_TRUE.  If CUBIT_TRUE,
    //-                    create RefEdges *instead* of splitting.

private:
    CubitStatus draw_preview(
        DLIList<Curve*> &curve_list,
        int color=CUBIT_BLUE);
    //- Draw the curves

    CubitStatus draw_preview(
        Curve *curve_ptr,
        CubitBoolean flush=CUBIT_TRUE,
        int color=CUBIT_BLUE);
    //- Draw the curves

    Curve *create_sweep_curve(
        Curve *curve_in,
        double distance,
        double chord_tol,
        CubitBoolean iterate = CUBIT_FALSE);
    //- Create a sweep curve offset a distance given the input curve.

    Surface *create_sweep_section(
        Curve* path,
        double distance,
        CubitVector up_vector,
        double fraction = 0.0);
    //- Uses two arcs to create a surface section normal to the input curve

    BodySM* create_sweep_body(
        Surface* section,
        Curve* path);
    //- Sweeps the input surface along the curve path

    CubitStatus create_ref_edges(
        DLIList<Curve*> &curve_list );
    //- Creates RefEdges from the list of Curve entities

    DLIList<Curve*> create_divide_curves(
        DLIList<Curve*> &curve_list,
        Surface* surface_in,
        double distance);
    //- Creates perpendicular curves of length 'distance' at the ends of the
    //- input curves. The curves are projected to the input surface if the 
    //- Surface in not of type 'PLANE_SURFACE_TYPE'.

    DLIList<Curve*> sourceCurves;
    //- List of the Curve entities to offset.

    static double tolerance;
    //- Linear tolerace
};

#endif

