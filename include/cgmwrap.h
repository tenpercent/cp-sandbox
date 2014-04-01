#ifndef CGMWRAP_H
#define CGMWRAP_H
typedef void* CGMmodel;
typedef void* CGMvertex;
typedef void* CGMedge;
typedef void* CGMcoedge;
typedef void* CGMloop;
typedef void* CGMface;
typedef void* CGMvolume;

/* Init */
int       CGM_Init(void);
CGMmodel  CGM_LoadOCCModelFromFile(const char *fname, const char *format);
CGMmodel  CGM_LoadOCCModelFromTopoDS(const void *shape);
CGMmodel  CGM_LoadModelFromFile(char *fname);

/* Number of <...> */
int       CGM_NumVertices(CGMmodel model);
int       CGM_NumEdges(CGMmodel model);
int       CGM_NumFaces(CGMmodel model);
int       CGM_NumVolumes(CGMmodel model);

/* Get ith <...> from global list */
CGMvertex CGM_ithVertex(CGMmodel model, int i);
CGMedge   CGM_ithEdge(CGMmodel model, int i);
CGMface   CGM_ithFace(CGMmodel model, int i);
CGMvolume CGM_ithVolume(CGMmodel model, int i);

/* Get index in global list */
int       CGM_VertexIndex(CGMmodel model, CGMvertex vertex);
int       CGM_EdgeIndex(CGMmodel model, CGMedge edge);
int       CGM_FaceIndex(CGMmodel model, CGMface face);
int       CGM_VolumeIndex(CGMmodel model, CGMvolume volume);

/* Get first <...> */
CGMvertex CGM_FirstVertex(CGMmodel model);
CGMedge   CGM_FirstEdge(CGMmodel model);
CGMface   CGM_FirstFace(CGMmodel model);
CGMvolume CGM_FirstVolume(CGMmodel model);

/* Get next <...> */
CGMvertex CGM_NextVertex(CGMmodel model, int *idx);
CGMedge   CGM_NextEdge(CGMmodel model, int *idx);
CGMface   CGM_NextFace(CGMmodel model, int *idx);
CGMvolume CGM_NextVolume(CGMmodel model, int *idx);

/* Get last <...> */
CGMvertex CGM_LastVertex(CGMmodel model);
CGMedge   CGM_LastEdge(CGMmodel model);
CGMface   CGM_LastFace(CGMmodel model);
CGMvolume CGM_LastVolume(CGMmodel model);

/* XYZ & parametrization */
void      CGM_GetVertexCoords(CGMvertex vertex, double outxyz[3]);
void      CGM_GetEdgeCoordsFromU(CGMedge edge, double u, double outxyz[3]);
double    CGM_GetEdgeUFromCoords(CGMedge edge, double xyz[3]);
void      CGM_GetFaceCoordsFromUV(CGMface face, double uv[2], double outxyz[3]);
void      CGM_GetFaceUVFromCoords(CGMface face, double xyz[3], double outuv[2]);

/* Edges */
int       CGM_NumFacesWithEdge(CGMedge edge);
CGMface   CGM_ithFaceWithEdge(CGMedge edge, int i);
CGMface   CGM_NextFaceWithEdge(CGMedge edge, int *idx);
CGMvertex CGM_EdgeStartVertex(CGMedge edge);
CGMvertex CGM_EdgeEndVertex(CGMedge edge);
void      CGM_GetEdgeTangent(CGMedge edge, double location[3], double tangent[3]);
void      CGM_ClosestPointOnEdge(CGMedge edge, double inxyz[3], double outxyz[3]);
double    CGM_EdgeStartParam(CGMedge edge);
double    CGM_EdgeEndParam(CGMedge edge);
double    CGM_EdgeLength(CGMedge edge);
double    CGM_EdgeArcLength(CGMedge edge, double p1[3], double p2[3]);
double    CGM_EdgeArcLengthParam(CGMedge edge, double u1, double u2);
int       CGM_FaceEdgeOrientation(CGMface face, CGMedge edge);
int       CGM_EdgeIsLinear(CGMedge edge);

/* Faces */
int       CGM_NumVolumesWithFace(CGMface face);
CGMvolume CGM_ithVolumeWithFace(CGMface face, int i);
CGMvolume CGM_NextVolumeWithFace(CGMface face, int *idx);
int       CGM_NumLoopsInFace(CGMface face);
CGMloop   CGM_ithLoopInFace(CGMface face, int i);
CGMloop   CGM_NextLoopInFace(CGMface face, int *idx);
int       CGM_NumEdgesInLoop(CGMloop loop);
CGMedge   CGM_ithEdgeInLoop(CGMloop loop, int i);
CGMedge   CGM_NextEdgeInLoop(CGMloop loop, int *idx);
int       CGM_NumEdgesInFace(CGMface face);
CGMedge   CGM_ithEdgeInFace(CGMface face, int i);
CGMedge   CGM_NextEdgeInFace(CGMface face, int *idx);
void      CGM_GetFaceNormal(CGMface face, double inoutposition[3], double normal[3]);
void      CGM_ClosestPointOnFace(CGMface face, double inxyz[3], double outxyz[3]);
void      CGM_GetFaceParamRange(CGMface face, double *umin, double *umax, double *vmin, double *vmax);
void      CGM_FaceDerivatives(CGMface face, double u, double v, double *du, double *dv);
void      CGM_FaceCurvatures(CGMface face, double xyz[3], double *c1, double *c2);
double    CGM_FaceArea(CGMface face);
int       CGM_VolumeFaceOrientation(CGMvolume volume, CGMface face);
int       CGM_GeometryFaceOrientation(CGMface face);
int       CGM_FaceIsPlanar(CGMface face);
double    CGM_FaceIsPeriodicU(CGMface face);
double    CGM_FaceIsPeriodicV(CGMface face);

/* Volumes */
int       CGM_NumFacesInVolume(CGMvolume volume);
CGMface   CGM_ithFaceInVolume(CGMvolume volume, int i);
CGMface   CGM_NextFaceInVolume(CGMvolume volume, int *idx);
double    CGM_Volume(CGMvolume volume);
void      CGM_VolumeBoundingBox(CGMvolume volume, double min[3], double max[3]);
void      CGM_ModelBoundingBox(CGMmodel model, double min[3], double max[3]);

/* Angles */
double    CGM_AngleOnFace(CGMmodel model, CGMface face, CGMedge edge1, CGMedge edge2);
double    CGM_AngleBetweenFaces(CGMmodel model, CGMface face1, CGMface face2, CGMedge edge);
#endif


