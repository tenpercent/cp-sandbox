#include <math.h>
#include "AppUtil.hpp"
#include "CGMApp.hpp"
#include "GeometryQueryTool.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitObserver.hpp"
//#include "CubitPointData.hpp"
//#include "CubitFacetData.hpp"
#include "DLIList.hpp"
#include "Surface.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "Body.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CubitBox.hpp"

extern "C" {
#include "cgmwrap.h"
}

#define CGMA 1

/* Init */
extern "C" int CGM_Init(void) {
	AppUtil::instance()->startup(0, 0);
	CGMApp::instance()->startup(0, 0);
	return 0;
}
extern "C" CGMmodel CGM_LoadOCCModelFromFile(const char *fname, const char *format) {
	OCCQueryEngine *oqe = OCCQueryEngine::instance();
	GeometryQueryTool *gqt = GeometryQueryTool::instance(oqe);
	gqt->delete_geometry();
	if (gqt->import_solid_model(fname, format)!=CUBIT_SUCCESS)   printf("failed\n");
	if (1) {
	    PRINT_INFO("Number of vertices = %d\n", gqt->num_ref_vertices());
	    PRINT_INFO("Number of edges = %d\n", gqt->num_ref_edges());
	    PRINT_INFO("Number of faces = %d\n", gqt->num_ref_faces());
	    PRINT_INFO("Number of volumes = %d\n", gqt->num_ref_volumes());
	    PRINT_INFO("Number of bodies = %d\n", gqt->num_bodies());
	}

	if (0) {
		DLIList<RefVertex*> verts;
		gqt->ref_vertices(verts);
		int i;
		for (i = 0; i < verts.size(); i++) {
			CubitVector coords = verts.get_and_step()->coordinates();
			PRINT_INFO("Vertex %d: %4.2f, %4.2f, %4.2f.\n", i, coords.x(), coords.y(), coords.z());
		}

		DLIList<RefEdge*> edges;
		gqt->ref_edges(edges);
		for (i = 0; i < edges.size(); i++) {
			RefEdge *edge = edges.get_and_step();
			PRINT_INFO("Edge %d (0x%lx): %4.2f, %4.2f.\n", i, (long)edge,  edge->start_param(), edge->end_param());
		}
		DLIList<RefFace*> faces;
		gqt->ref_faces(faces);
		for (i = 0; i < faces.size(); i++) {
			RefFace *face = faces.get_and_step();
			PRINT_INFO("Face %d (0x%lx):\n", i, (long)face);
		}
	}
	return gqt;
}
extern "C" CGMmodel CGM_LoadOCCModelFromTopoDS(const void *shape) {
	OCCQueryEngine *oqe = OCCQueryEngine::instance();
	GeometryQueryTool *gqt = GeometryQueryTool::instance(oqe);
	gqt->delete_geometry();
	DLIList<TopologyBridge*> bridge_list;
	bridge_list = oqe->populate_topology_bridge(*((TopoDS_Shape*) shape));
	gqt->construct_refentities(bridge_list, NULL);
	if (1) {
	    PRINT_INFO("Number of vertices = %d\n", gqt->num_ref_vertices());
	    PRINT_INFO("Number of edges = %d\n", gqt->num_ref_edges());
	    PRINT_INFO("Number of faces = %d\n", gqt->num_ref_faces());
	    PRINT_INFO("Number of volumes = %d\n", gqt->num_ref_volumes());
	    PRINT_INFO("Number of bodies = %d\n", gqt->num_bodies());
	}
	return gqt;
}
extern "C" CGMmodel CGM_LoadModelFromFile(char *fname) {
	if (CGMA)   return CGM_LoadOCCModelFromFile(fname, "OCC");
	if (!CGMA)  return CGM_LoadOCCModelFromFile(fname, "BREP");
	return NULL;
}

/* Number of <...> */
extern "C" int CGM_NumVertices(CGMmodel model) {
	return ((GeometryQueryTool*)model)->num_ref_vertices();
}
extern "C" int CGM_NumEdges(CGMmodel model) {
	return ((GeometryQueryTool*)model)->num_ref_edges();
}
extern "C" int CGM_NumFaces(CGMmodel model) {
	return ((GeometryQueryTool*)model)->num_ref_faces();
}
extern "C" int CGM_NumVolumes(CGMmodel model) {
	return ((GeometryQueryTool*)model)->num_ref_volumes();
}

/* Get ith <...> from list */
extern "C" CGMvertex CGM_ithVertex(CGMmodel model, int i) {
	DLIList<RefVertex*> list;
	((GeometryQueryTool*)model)->ref_vertices(list);
	return list[i];
}
extern "C" CGMedge CGM_ithEdge(CGMmodel model, int i) {
	DLIList<RefEdge*> list;
	((GeometryQueryTool*)model)->ref_edges(list);
	return list[i];
}
extern "C" CGMface CGM_ithFace(CGMmodel model, int i) {
	DLIList<RefFace*> list;
	((GeometryQueryTool*)model)->ref_faces(list);
	return list[i];
}
extern "C" CGMvolume CGM_ithVolume(CGMmodel model, int i) {
	DLIList<RefVolume*> list;
	((GeometryQueryTool*)model)->ref_volumes(list);
	return list[i];
}
/* Get index in list */
extern "C" int CGM_VertexIndex(CGMmodel model, CGMvertex vertex) {
	DLIList<RefVertex*> list;
	((GeometryQueryTool*)model)->ref_vertices(list);
	return list.where_is_item((RefVertex*)vertex);
}
extern "C" int CGM_EdgeIndex(CGMmodel model, CGMedge edge) {
	DLIList<RefEdge*> list;
	((GeometryQueryTool*)model)->ref_edges(list);
	return list.where_is_item((RefEdge*)edge);
}
extern "C" int CGM_FaceIndex(CGMmodel model, CGMface face) {
	DLIList<RefFace*> list;
	((GeometryQueryTool*)model)->ref_faces(list);
	return list.where_is_item((RefFace*)face);
}
extern "C" int CGM_VolumeIndex(CGMmodel model, CGMvolume volume) {
	DLIList<RefVolume*> list;
	((GeometryQueryTool*)model)->ref_volumes(list);
	return list.where_is_item((RefVolume*)volume);
}
/* Get first <...> */
extern "C" CGMvertex CGM_FirstVertex(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_first_ref_vertex();
}
extern "C" CGMedge CGM_FirstEdge(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_first_ref_edge();
}
extern "C" CGMface CGM_FirstFace(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_first_ref_face();
}
extern "C" CGMvolume CGM_FirstVolume(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_first_ref_volume();
}
/* Get next <...> */
/*///NOT safe
extern "C" CGMvertex CGM_NextVertex(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_next_ref_vertex();
}
extern "C" CGMedge CGM_NextEdge(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_next_ref_edge();
}
extern "C" CGMface CGM_NextFace(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_next_ref_face();
}
extern "C" CGMvolume CGM_NextVolume(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_next_ref_volume();
}
// */
extern "C" CGMvertex CGM_NextVertex(CGMmodel model, int *idx) {
	return CGM_ithVertex(model, (*idx)++);
}
extern "C" CGMedge CGM_NextEdge(CGMmodel model, int *idx) {
	return CGM_ithEdge(model, (*idx)++);
}
extern "C" CGMface CGM_NextFace(CGMmodel model, int *idx) {
	return CGM_ithFace(model, (*idx)++);
}
extern "C" CGMvolume CGM_NextVolume(CGMmodel model, int *idx) {
	return CGM_ithVolume(model, (*idx)++);
}
/* Get last <...> */
extern "C" CGMvertex CGM_LastVertex(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_last_ref_vertex();
}
extern "C" CGMedge CGM_LastEdge(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_last_ref_edge();
}
extern "C" CGMface CGM_LastFace(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_last_ref_face();
}
extern "C" CGMvolume CGM_LastVolume(CGMmodel model) {
	return ((GeometryQueryTool*)model)->get_last_ref_volume();
}
/* XYZ & parametrization */
extern "C" void CGM_GetVertexCoords(CGMvertex vertex, double outxyz[3]) {
	CubitVector p = ((RefVertex*)vertex)->coordinates();
	outxyz[0] = p.x();  outxyz[1] = p.y();  outxyz[2] = p.z();
}
extern "C" void CGM_GetEdgeCoordsFromU(CGMedge edge, double u, double outxyz[3]) {
	CubitVector p;
	((RefEdge*)edge)->position_from_u(u, p);
	outxyz[0] = p.x();  outxyz[1] = p.y();  outxyz[2] = p.z();
}
extern "C" double CGM_GetEdgeUFromCoords(CGMedge edge, double xyz[3]) {
	CubitVector p = xyz;
	return ((RefEdge*)edge)->u_from_position(p);
}
extern "C" void CGM_GetFaceCoordsFromUV(CGMface face, double uv[2], double outxyz[3]) {
	CubitVector p = ((RefFace*)face)->position_from_u_v(uv[0], uv[1]);
	outxyz[0] = p.x();  outxyz[1] = p.y();  outxyz[2] = p.z();
}
extern "C" void CGM_GetFaceUVFromCoords(CGMface face, double xyz[3], double outuv[2]) {
	int r;
	CubitVector p = xyz;
	r = ((RefFace*)face)->u_v_from_position(p, outuv[0], outuv[1]);
	if (r != CUBIT_SUCCESS)  printf("CGM_GetFaceUVFromCoords() Failed! [code = %d]\n", r);
}
/* Edges */
extern "C" int       CGM_NumFacesWithEdge(CGMedge edge) {
	return ((RefEdge*)edge)->num_ref_faces();
}
extern "C" CGMface   CGM_ithFaceWithEdge(CGMedge edge, int i) {
	DLIList<RefFace*> list;
	((RefEdge*)edge)->ref_faces(list);
	return list[i];
}
extern "C" CGMface   CGM_NextFaceWithEdge(CGMedge edge, int *idx) {
	return CGM_ithFaceWithEdge(edge, (*idx)++);
}
extern "C" CGMvertex CGM_EdgeStartVertex(CGMedge edge) {
	return ((RefEdge*)edge)->start_vertex();
}
extern "C" CGMvertex CGM_EdgeEndVertex(CGMedge edge) {
	return ((RefEdge*)edge)->end_vertex();
}
extern "C" void CGM_GetEdgeTangent(CGMedge edge, double location[3], double tangent[3]) {
	CubitVector p = location;
	CubitVector tan;
	((RefEdge*)edge)->tangent(p, tan);
	tangent[0] = tan.x();  tangent[1] = tan.y();  tangent[2] = tan.z();
}
extern "C" void CGM_ClosestPointOnEdge(CGMedge edge, double inxyz[3], double outxyz[3]) {
	CubitVector p = inxyz;
	((RefEdge*)edge)->move_to_curve(p);
	outxyz[0] = p.x();  outxyz[1] = p.y(); outxyz[2] = p.z();
}
extern "C" double CGM_EdgeStartParam(CGMedge edge) {
	return ((RefEdge*)edge)->start_param();
}
extern "C" double CGM_EdgeEndParam(CGMedge edge) {
	return ((RefEdge*)edge)->end_param();
}
extern "C" double CGM_EdgeLength(CGMedge edge) {
	return ((RefEdge*)edge)->measure();
}
extern "C" double CGM_EdgeArcLength(CGMedge edge, double p1[3], double p2[3]) {
	CubitVector v1=p1;
	CubitVector v2=p2;
	return ((RefEdge*)edge)->get_arc_length(v1, v2); 
}
extern "C" double    CGM_EdgeArcLengthParam(CGMedge edge, double u1, double u2) {
	return ((RefEdge*)edge)->length_from_u(u1, u2); 
}
extern "C" int CGM_FaceEdgeOrientation(CGMface face, CGMedge edge) {
	CubitSense orient = ((RefEdge*)edge)->sense(((RefFace*)face));
	if (orient == CUBIT_FORWARD) return 0;
	else if (orient == CUBIT_REVERSED) return 1;
	else return -1;
}
extern "C" int CGM_EdgeIsLinear(CGMedge edge) {
	CubitVector d1, d2;
	if (((RefEdge*)edge)->get_point_direction(d1, d2) == CUBIT_FAILURE) return 0;
	else return 1;
}
/* Faces */
extern "C" int       CGM_NumVolumesWithFace(CGMface face) {
	return ((RefFace*)face)->num_ref_volumes();
}
extern "C" CGMvolume CGM_ithVolumeWithFace(CGMface face, int i) {
	DLIList<RefVolume*> list;
	((RefFace*)face)->ref_volumes(list);
	return list[i];
}
extern "C" CGMvolume CGM_NextVolumeWithFace(CGMface face, int *idx) {
	return CGM_ithVolumeWithFace(face, (*idx)++);
}
extern "C" int CGM_NumLoops(CGMface face) {
	return ((RefFace*)face)->number_of_Loops();
}
extern "C" CGMloop CGM_ithLoop(CGMface face, int i) {
	DLIList<DLIList<RefEdge*>*> list;
	((RefFace*)face)->ref_edge_loops(list);
	return list[i];
}
extern "C" CGMloop CGM_NextLoop(CGMface face, int *idx) {
	return CGM_ithLoop(face, (*idx)++);
}
extern "C" int CGM_NumEdgesInLoop(CGMloop loop) {
	return ((DLIList<RefEdge*>*)loop)->size();
}
extern "C" CGMedge CGM_ithEdgeInLoop(CGMloop loop, int i) {
	return (*((DLIList<RefEdge*>*)loop))[i];
}
extern "C" CGMedge CGM_NextEdgeInLoop(CGMloop loop, int *idx) {
	return CGM_ithEdgeInLoop(loop, (*idx)++);
}
extern "C" int CGM_NumEdgesInFace(CGMface face) {
	return ((RefFace*)face)->num_ref_edges();
}
extern "C" CGMedge CGM_ithEdgeInFace(CGMface face, int i) {
	DLIList<RefEdge*> list;
	((RefFace*)face)->ref_edges(list);
	return list[i];
}
extern "C" CGMedge CGM_NextEdgeInFace(CGMface face, int *idx) {
	return CGM_ithEdgeInFace(face, (*idx)++);
}
extern "C" void CGM_GetFaceNormal(CGMface face, double inoutposition[3], double normal[3]) {
	CubitVector p = inoutposition;
	CubitVector n = ((RefFace*)face)->normal_at(p);
	normal[0] = n.x();  normal[1] = n.y(); normal[2] = n.z();
}
extern "C" void CGM_ClosestPointOnFace(CGMface face, double inxyz[3], double outxyz[3]) {
	CubitVector p = inxyz;
	((RefFace*)face)->move_to_surface(p);
	outxyz[0] = p.x();  outxyz[1] = p.y(); outxyz[2] = p.z();
}
extern "C" void CGM_GetFaceParamRange(CGMface face, double *umin, double *umax, double *vmin, double *vmax) {
	((RefFace*)face)->get_param_range_U(*umin, *umax);
	((RefFace*)face)->get_param_range_V(*vmin, *vmax);
}
/*///undefined reference to `RefFace::uv_derivitives(double, double, CubitVector&, CubitVector&)
extern "C" void CGM_FaceDerivatives(CGMface face, double u, double v, double du[3], double dv[3]) {
	CubitVector pdu;
	CubitVector pdv;
	((RefFace*)face)->uv_derivitives(u, v, pdu, pdv);
	du[0] = pdu.x();  du[1] = pdu.y();  du[2] = pdu.z();
	dv[0] = pdv.x();  dv[1] = pdv.y();  dv[2] = pdv.z();
}
// */
extern "C" void CGM_FaceCurvatures(CGMface face, double xyz[3], double *c1, double *c2) {
	CubitVector p = xyz;
	((RefFace*)face)->get_principal_curvatures(p, *c1, *c2);
}
extern "C" double CGM_FaceArea(CGMface face) {
	return ((RefFace*)face)->area();
}
extern "C" int CGM_VolumeFaceOrientation(CGMvolume volume, CGMface face) {
	CubitSense orient = ((RefFace*)face)->sense(((RefVolume*)volume));
	if (orient == CUBIT_FORWARD) return 0;
	else if (orient == CUBIT_REVERSED) return 1;
	else return -1;
}
/*///undefined reference to `RefFace::get_geometry_sense() 
extern "C" int CGM_GeneralFaceOrientation(CGMface face) {
	CubitSense orient = ((RefFace*)face)->get_geometry_sense();
	if (orient == CUBIT_FORWARD) return 0;
	else if (orient == CUBIT_REVERSED) return 1;
	else return -1;
}
// */
extern "C" int CGM_FaceIsPlanar(CGMface face) {
	if (((RefFace*)face)->is_planar()) return 1;
	else return 0;
}
extern "C" double CGM_FaceIsPeriodicU(CGMface face) {
	double period=0.0;
	if (((RefFace*)face)->is_periodic_in_U(period))  return period;
	else return 0.0;
}
extern "C" double CGM_FaceIsPeriodicV(CGMface face) {
	double period=0.0;
	if (((RefFace*)face)->is_periodic_in_V(period))  return period;
	else return 0.0;
}

/* Volumes */
/*///no shells yet
extern "C" int       CGM_NumShells(CGMvolume volume);
extern "C" CGMshell  CGM_ithShellInVolume(CGMvolume volume, int i);
extern "C" CGMshell  CGM_NextShellInVolume(CGMvolume volume, int *idx);
extern "C" int       CGM_NumFacesInShell(CGMshell shell);
extern "C" CGMface   CGM_ithFaceInShell(CGMshell shell, int i);
extern "C" CGMface   CGM_NextFaceInShell(CGMshell shell, int *idx);
// */
extern "C" int       CGM_NumFacesInVolume(CGMvolume volume) {
	return ((RefVolume*)volume)->num_ref_faces();
}
extern "C" CGMface   CGM_ithFaceInVolume(CGMvolume volume, int i) {
	DLIList<RefFace*> list;
	((RefVolume*)volume)->ref_faces(list);
	return list[i];
}
extern "C" CGMface   CGM_NextFaceInVolume(CGMvolume volume, int *idx) {
	return CGM_ithFaceInVolume(volume, (*idx)++);
}
extern "C" double    CGM_Volume(CGMvolume volume) {
	return ((RefVolume*)volume)->measure();
}
extern "C" void      CGM_VolumeBoundingBox(CGMvolume volume, double min[3], double max[3]) {
	CubitBox box=((RefVolume*)volume)->bounding_box();
	CubitVector vmin = box.minimum();
	CubitVector vmax = box.maximum();
	min[0]=vmin.x(); min[1]=vmin.y(); min[2]=vmin.z();
	max[0]=vmax.x(); max[1]=vmax.y(); max[2]=vmax.z();
}
extern "C" void      CGM_ModelBoundingBox(CGMmodel model, double min[3], double max[3]) {
	CubitBox box=((GeometryQueryTool*)model)->model_bounding_box();
	CubitVector vmin = box.minimum();
	CubitVector vmax = box.maximum();
	min[0]=vmin.x(); min[1]=vmin.y(); min[2]=vmin.z();
	max[0]=vmax.x(); max[1]=vmax.y(); max[2]=vmax.z();
}

/* Angles */ 
extern "C" double CGM_AngleOnFace(CGMmodel model, CGMface face, CGMedge edge1, CGMedge edge2) {
	return ((GeometryQueryTool*)model)->geometric_angle(((RefEdge*)edge1), ((RefEdge*)edge2), ((RefFace*)face));
}
extern "C" double CGM_AngleBetweenFaces(CGMmodel model, CGMface face1, CGMface face2, CGMedge edge) {
	return ((GeometryQueryTool*)model)->surface_angle(((RefFace*)face1), ((RefFace*)face2), ((RefEdge*)edge));
}
