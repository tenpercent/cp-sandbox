//- Class:       AnalyticGeometryTool
//- Description: This class performs calculations on analytic geometry
//               (points, lines, arcs, planes, polygons).  Capabilities 
//               include vector and point math, matrix operations, 
//               measurements, intersections and comparison/containment checks.
// 
//- Owner: Steve Storm, Caterpillar Inc.
//- Checked by:


#include <float.h>
#include "AnalyticGeometryTool.hpp"
#include "CubitMessage.hpp"
#include "CubitBox.hpp"
#include "CubitPlane.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include <math.h>
#include "CastTo.hpp"

#include <fstream>
using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
using std::ostream;

static double AGT_IDENTITY_MTX_4X4[4][4] = { {1.0, 0.0, 0.0, 0.0},
                                             {0.0, 1.0, 0.0, 0.0},
                                             {0.0, 0.0, 1.0, 0.0},
                                             {0.0, 0.0, 0.0, 1.0} };

static double AGT_IDENTITY_MTX_3X3[3][3] = { {1.0, 0.0, 0.0},
                                             {0.0, 1.0, 0.0},
                                             {0.0, 0.0, 1.0} };

AnalyticGeometryTool* AnalyticGeometryTool::instance_ = 0;

Point2 operator+ (const Point2& p, const Point2& q)
{
    Point2 add;
    add.x = p.x + q.x;
    add.y = p.y + q.y;
    return add;
}

Point2 operator- (const Point2& p, const Point2& q)
{
    Point2 sub;
    sub.x = p.x - q.x;
    sub.y = p.y - q.y;
    return sub;
}

Point2 operator* (double t, const Point2& p)
{
    Point2 prod;
    prod.x = t*p.x;
    prod.y = t*p.y;
    return prod;
}

Point2 operator* (const Point2& p, double t)
{
    Point2 prod;
    prod.x = t*p.x;
    prod.y = t*p.y;
    return prod;
}

Point2 operator- (const Point2& p)
{
    Point2 neg;
    neg.x = -p.x;
    neg.y = -p.y;
    return neg;
}

inline double AnalyticGeometryTool::Dot (const Point2& p, const Point2& q)
{
    return double(p.x*q.x + p.y*q.y);
}

Point3 operator+ (const Point3& p, const Point3& q)
{
    Point3 add;
    add.x = p.x + q.x;
    add.y = p.y + q.y;
    add.z = p.z + q.z;
    return add;
}

Point3 operator- (const Point3& p, const Point3& q)
{
    Point3 sub;
    sub.x = p.x - q.x;
    sub.y = p.y - q.y;
    sub.z = p.z - q.z;
    return sub;
}

Point3 operator* (double t, const Point3& p)
{
    Point3 prod;
    prod.x = t*p.x;
    prod.y = t*p.y;
    prod.z = t*p.z;
    return prod;
}

Point3 operator* (const Point3& p, double t)
{
    Point3 prod;
    prod.x = t*p.x;
    prod.y = t*p.y;
    prod.z = t*p.z;
    return prod;
}

Point3 operator- (const Point3& p)
{
    Point3 neg;
    neg.x = -p.x;
    neg.y = -p.y;
    neg.z = -p.z;
    return neg;
}

// Method: instance
// provides access to the unique model for this execution.
// sets up this instance on first access
AnalyticGeometryTool* AnalyticGeometryTool::instance()
{
  if (instance_ == 0) {
    instance_ = new AnalyticGeometryTool;
  }
  return instance_;
}

AnalyticGeometryTool::AnalyticGeometryTool() 
{
   agtEpsilon = 1e-8;
}

AnalyticGeometryTool::~AnalyticGeometryTool()
{
}

//***************************************************************************
// Double numbers
//***************************************************************************


void AnalyticGeometryTool::round_near_val( double &val )
{
   if( dbl_eq( val, 0.0 ) )
      val = 0.0;
   else if( dbl_eq( val, 1.0 ) )
      val  = 1.0;
   else if( dbl_eq( val, -1.0 ) )
      val = -1.0;
}

//***************************************************************************
// Matrices & Transforms
//***************************************************************************
void AnalyticGeometryTool::transform_pnt( double m[4][4], 
                                         double pin[3],
                                         double pout[3] )
{
    double p[3];  // working buffer
    
    // Check if transformation can occur 
    if (!m) {
       if (pin && pout)
	  copy_pnt(pin, pout);
       return;
    }
    
    // Perform transformation    
    p[0] = m[0][0] * pin[0] + m[1][0] * pin[1] + m[2][0] * pin[2] + m[3][0];
    p[1] = m[0][1] * pin[0] + m[1][1] * pin[1] + m[2][1] * pin[2] + m[3][1];
    p[2] = m[0][2] * pin[0] + m[1][2] * pin[1] + m[2][2] * pin[2] + m[3][2];
    
    // Copy work buffer to out point    
    copy_pnt(p,pout);
} 

void AnalyticGeometryTool::transform_vec( double m3[3][3],
                                         double vin[3],
                                         double vout[3] )
{
   double v[3];    // working buffer
   
   // Determine if transformation can occur 
   if (!m3) {      
      if (vin && vout)
         copy_vec(vin, vout);
      return;
   }
   
   // Perform transformation 
   v[0] = m3[0][0] * vin[0] + m3[1][0] * vin[1] + m3[2][0] * vin[2];
   v[1] = m3[0][1] * vin[0] + m3[1][1] * vin[1] + m3[2][1] * vin[2];
   v[2] = m3[0][2] * vin[0] + m3[1][2] * vin[1] + m3[2][2] * vin[2];
   
   // Copy work buffer to vector out    
   copy_pnt(v,vout);
}

void AnalyticGeometryTool::transform_vec( double m4[4][4], 
                                         double vin[3],
                                         double vout[3] )
{
   double v[3];    // working buffer
   
   // Determine if transformation can occur 
   if (!m4) {      
      if (vin && vout)
         copy_vec(vin, vout);
      return;
   }
   
   // Perform transformation
   
   v[0] = m4[0][0] * vin[0] + m4[1][0] * vin[1] + m4[2][0] * vin[2];
   v[1] = m4[0][1] * vin[0] + m4[1][1] * vin[1] + m4[2][1] * vin[2];
   v[2] = m4[0][2] * vin[0] + m4[1][2] * vin[1] + m4[2][2] * vin[2];
   
   // Copy work buffer to vector out    
   copy_pnt(v,vout);
}

void AnalyticGeometryTool::transform_line( double rot_mtx[4][4],
                                           double origin[3], double axis[3] )
{
   double end_pnt[3]; // Find arbitrary end point on line
   next_pnt( origin, axis, 10.0, end_pnt );
   
   transform_pnt( rot_mtx, origin, origin );
   transform_pnt( rot_mtx, end_pnt, end_pnt );
   
   axis[0] = end_pnt[0] - origin[0];
   axis[1] = end_pnt[1] - origin[1];
   axis[2] = end_pnt[2] - origin[2];
}

void AnalyticGeometryTool::transform_line( double rot_mtx[4][4],
                                          CubitVector &origin, CubitVector &axis )
{
   CubitVector end_point; // Find arbitrary end point on line
   origin.next_point( axis, 10.0, end_point );
   
   double start_pnt[3], end_pnt[3];
   copy_pnt( origin, start_pnt );
   copy_pnt( end_point, end_pnt );
   
   transform_pnt( rot_mtx, start_pnt, start_pnt );
   transform_pnt( rot_mtx, end_pnt, end_pnt );
   
   axis.x( end_pnt[0] - start_pnt[0] );
   axis.y( end_pnt[1] - start_pnt[1] );
   axis.z( end_pnt[2] - start_pnt[2] );
   
   origin.set( start_pnt );
}

void AnalyticGeometryTool::copy_mtx( double from[3][3],double to[3][3] )
{   
   // Determine if identity matrix needed 
   if (!from) // copy in the identity matrix 
      memcpy(to, AGT_IDENTITY_MTX_3X3, sizeof(double)*9); 
   else // Copy from to to 
      memcpy(to, from, sizeof(double)*9);   
}

void AnalyticGeometryTool::copy_mtx( double from[4][4], double to[4][4] )
{   
   // determine if identity matrix needed 
   if (!from) // copy in the identity matrix 
      memcpy(to, AGT_IDENTITY_MTX_4X4, sizeof(double)*16);  
   else // copy from to to 
      memcpy(to, from, sizeof(double)*16);   
}

void AnalyticGeometryTool::copy_mtx( double from[4][4], double to[3][3] )
{   
   size_t dbl3;
   
   dbl3 = sizeof(double) * 3;
   
   // Determine if identity matrix needed 
   if (!from) // Copy in the identity matrix 
      memcpy(to, AGT_IDENTITY_MTX_3X3, sizeof(double)*9); 
   else { // Copy each upper left element of from to to 
      memcpy(to[0], from[0], dbl3);
      memcpy(to[1], from[1], dbl3);
      memcpy(to[2], from[2], dbl3);
   }   
}

void AnalyticGeometryTool::copy_mtx(double from[3][3], double to[4][4],
                                    double* origin )
{   
   size_t dbl3;
   
   dbl3 = sizeof(double) * 3;
   
   // Determine if identity matrix needed 
   if (!from) { // Copy in the identity matrix 
   
      memcpy(to, AGT_IDENTITY_MTX_4X4, sizeof(double)*16);
      
      if (origin)
         memcpy(to[3], origin, dbl3);
   }   
	     
   else { // Copy each upper element of from to to 
      
      memcpy(to[0], from[0], dbl3);
      memcpy(to[1], from[1], dbl3);
      memcpy(to[2], from[2], dbl3);
      
      to[0][3] = 0.0;
      to[1][3] = 0.0;
      to[2][3] = 0.0;
      to[3][3] = 1.0;
      
      if (origin)
      {
         memcpy(to[3], origin, dbl3);
//         to[3][0] = origin[0];
//         to[3][1] = origin[1];
//         to[3][2] = origin[2];
      }
      else
         memcpy(to[3], AGT_IDENTITY_MTX_4X4[3], sizeof(double)*3);
   } 
}

void AnalyticGeometryTool::create_rotation_mtx( double theta, double v[3],
                                               double mtx3x3[3][3] )
{        
   double coeff1;
   double coeff2;
   double coeff3;
   double v_unit[3];
   
   if (!mtx3x3)
      return;
   
   coeff1 = cos(theta);
   coeff2 = (1.0l - coeff1);
   coeff3 = sin(theta);
   
   unit_vec(v, v_unit);
   
   mtx3x3[0][0] = coeff1 + coeff2 * (v_unit[0] * v_unit[0]);
   mtx3x3[0][1] = coeff2 * v_unit[1] * v_unit[0] + coeff3 * v_unit[2];
   mtx3x3[0][2] = coeff2 * v_unit[2] * v_unit[0] - coeff3 * v_unit[1];
   
   mtx3x3[1][0] = coeff2 * v_unit[1] * v_unit[0] - coeff3 * v_unit[2];
   mtx3x3[1][1] = coeff1 + coeff2 * (v_unit[1] * v_unit[1]);
   mtx3x3[1][2] = coeff2 * v_unit[1] * v_unit[2] + coeff3 * v_unit[0];
   
   mtx3x3[2][0] = coeff2 * v_unit[2] * v_unit[0] + coeff3 * v_unit[1];
   mtx3x3[2][1] = coeff2 * v_unit[2] * v_unit[1] - coeff3 * v_unit[0];
   mtx3x3[2][2] = coeff1 + coeff2 * (v_unit[2] * v_unit[2]);
}

void AnalyticGeometryTool::create_rotation_mtx( double theta, double v[3],
                                               double mtx4x4[4][4] )
{        
   double coeff1;
   double coeff2;
   double coeff3;
   double v_unit[3];
   
   if (!mtx4x4)
      return;
   
   coeff1 = cos(theta);
   coeff2 = (1.0l - coeff1);
   coeff3 = sin(theta);

   unit_vec(v, v_unit);
   
   mtx4x4[0][0] = coeff1 + coeff2 * (v_unit[0] * v_unit[0]);
   mtx4x4[0][1] = coeff2 * v_unit[1] * v_unit[0] + coeff3 * v_unit[2];
   mtx4x4[0][2] = coeff2 * v_unit[2] * v_unit[0] - coeff3 * v_unit[1];
   
   mtx4x4[1][0] = coeff2 * v_unit[1] * v_unit[0] - coeff3 * v_unit[2];
   mtx4x4[1][1] = coeff1 + coeff2 * (v_unit[1] * v_unit[1]);
   mtx4x4[1][2] = coeff2 * v_unit[1] * v_unit[2] + coeff3 * v_unit[0];
   
   mtx4x4[2][0] = coeff2 * v_unit[2] * v_unit[0] + coeff3 * v_unit[1];
   mtx4x4[2][1] = coeff2 * v_unit[2] * v_unit[1] - coeff3 * v_unit[0];
   mtx4x4[2][2] = coeff1 + coeff2 * (v_unit[2] * v_unit[2]);
   
   mtx4x4[0][3] = 0.0;
   mtx4x4[1][3] = 0.0;
   mtx4x4[2][3] = 0.0;
   mtx4x4[3][3] = 1.0;
   mtx4x4[3][0] = 0.0;
   mtx4x4[3][1] = 0.0;
   mtx4x4[3][2] = 0.0;  
}

void AnalyticGeometryTool::add_origin_to_rotation_mtx( double rot_mtx[4][4], 
                                                      double origin[3] )
{  
   double tmp_mtx[4][4];

   // Translate to origin
   double t[4][4]; 
   memcpy(t, AGT_IDENTITY_MTX_4X4, sizeof(double)*16);
   //PRINT_INFO( "Rotation matrix, before origin: \n" );
   //print_mtx( rot_mtx );
   t[3][0]=-origin[0]; t[3][1]=-origin[1]; t[3][2]=-origin[2];
   mult_mtx( t, rot_mtx, tmp_mtx );
   
   //PRINT_INFO( "Origin times rotation: \n" );
   //print_mtx( tmp_mtx );
   
   // Translate back
   t[3][0]=origin[0]; t[3][1]=origin[1]; t[3][2]=origin[2];
   mult_mtx( tmp_mtx, t, rot_mtx );
   
   //PRINT_INFO( "Rotation x -origin: \n" );
   //print_mtx( rot_mtx );
}

void AnalyticGeometryTool::identity_mtx( double mtx3x3[3][3] )
{ 
   memcpy(mtx3x3, AGT_IDENTITY_MTX_3X3, sizeof(double)*9);
}

void AnalyticGeometryTool::identity_mtx( double mtx4x4[4][4] )
{ 
   memcpy(mtx4x4, AGT_IDENTITY_MTX_4X4, sizeof(double)*16);
}

void AnalyticGeometryTool::mtx_to_angs( double mtx3x3[3][3],
                                       double &ax, double &ay,
                                       double &az )
{   
//   METHOD:
//   o Rotate x-vector onto xz plane
//   o  Check xp dotted into y  
//   o   If xp dot y is zero  ==> az = 0 (x-vector already in xz plane)
//   o   Otherwise, compute rotation of vector into xz plane to acquire *az     
//   o    Use atan2 (on x-vector) to get *az
//   o    Rotate the system about z 
//   o Use atan2 function (on x-vector in xz plane) to determine *ay 
//   o Rotate the system about y      
//   o Compute ax using y-vector
//   o Resultant angles are negated (to reverse above procedure)
 
   double x[3];                    // x-axis vector 
   double y[3];                    // y-axis vector 
   double z[3];                    // z-axis vector 
   double ar[3][3];                // Rotation matrix 
   double sinr,cosr;               // Used for atan2 function 
   double work_sys[3][3];          // Temporary holder for system 
   double *xp = work_sys[0];       // x-axis vector of system: x-primed 
   double *yp = work_sys[1];       // y-axis vector of system: y-primed 
   
   x[0] = 1.0; x[1] = 0.0; x[2] = 0.0;
   y[0] = 0.0; y[1] = 1.0; y[2] = 0.0;
   z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;
   
   if (!mtx3x3)
      return;
   
   // Copy matrix over to work csys 
   copy_mtx(mtx3x3,work_sys);
      
   // Check xp dotted into y 
   //   If xp dot y is zero  ==> az = 0 
   
   // Otherwise, compute rotation of vector into xz plane to acquire *az 
   
   if (dbl_eq(dot_vec(xp,y), 0.0))
      
      az = 0.0;
	 
   else {
      
      /* 
      Compute *az - rotate xp to xz-plane about z-axis          
               y    xp                                           
               |  /                                              
               | /  negative angle about z (negative of atan2)   
                -----x              (use RH rule about z)        
                 \                                               
                  \ positive angle about z (negative of atan2)         
                   xp                                            
      */
	 
      sinr = dot_vec(xp,y);
      cosr = dot_vec(xp,x);
      az = - atan2(sinr, cosr);
      
      // Rotate the system about z 
      create_rotation_mtx(az,z,ar);
      rotate_system(ar,work_sys,work_sys);
      
   }
   
      /* 
      Compute *ay - rotate xp to x-axis about y-axis            
              z    xp                                           
              |  /                                              
              | /  positive angle about y (positive of atan2)   
               -----x     (use RH rule about y)                 
                \                                               
                 \ negative angle about y (positive of atan2)   
                  xp                                            
      */
   
   sinr = dot_vec(xp,z);
   cosr = dot_vec(xp,x);
   ay = atan2(sinr, cosr);
   
   // Rotate the system about y       
   create_rotation_mtx(ay,y,ar);      
   rotate_system(ar,work_sys,work_sys);
  
   /*
      Compute *ax - rotate yp to y-axis about x-axis            
              z    yp                                           
              |  /                                              
              | /  negative angle about x (negative of atan2)   
               -----y     (use RH rule about x)                 
                \                                               
                 \ positive angle about x (negative of atan2)   
                  yp                                            
   */
   
   sinr = dot_vec(yp,z);
   cosr = dot_vec(yp,y);
   ax = atan2(sinr,cosr);     // Negative of negative - see below 
   
   // Negate above angles for rotation of the system back to original 
   az = -(az);
   ay = -(ay);
   
   // Make sure near zero angles are actually zero 
   if (dbl_eq(ax, 0.0))
      ax = 0.0;
   
   if (dbl_eq(ay, 0.0))
      ay = 0.0;   
}

void AnalyticGeometryTool::mtx_to_angs( double mtx4x4[4][4],
                                       double &ax, double &ay,
                                       double &az )
{   
   double work_sys[3][3];

   if(!mtx4x4)
      return;

   copy_mtx(mtx4x4,work_sys);
   mtx_to_angs( work_sys, ax, ay, az );
}

void AnalyticGeometryTool::rotate_system( double mtx[3][3], 
                                         double sys_in[3][3],
                                         double sys_out[3][3] )
{
   double sys_tmp[3][3];
   double *p_sys_tmp;
   
   // Check to see if rotating in place 
   if (sys_in == sys_out) {
      copy_mtx( sys_in, sys_tmp );
      p_sys_tmp = (double *)sys_tmp; 
   }
   else 
      // Have sys_tmp point at outgoing memory location 
      p_sys_tmp = (double *)sys_out;
   
   
   // X-vector 
   p_sys_tmp[0] =   mtx[0][0] * sys_in[0][0]
                  + mtx[1][0] * sys_in[0][1]
                  + mtx[2][0] * sys_in[0][2];
   p_sys_tmp[1] =   mtx[0][1] * sys_in[0][0]
                  + mtx[1][1] * sys_in[0][1]
                  + mtx[2][1] * sys_in[0][2];
   p_sys_tmp[2] =   mtx[0][2] * sys_in[0][0]
                  + mtx[1][2] * sys_in[0][1]
                  + mtx[2][2] * sys_in[0][2];   
   
   // Y-vector 
   p_sys_tmp[3] =   mtx[0][0] * sys_in[1][0]
                  + mtx[1][0] * sys_in[1][1]
                  + mtx[2][0] * sys_in[1][2];
   p_sys_tmp[4] =   mtx[0][1] * sys_in[1][0]
                  + mtx[1][1] * sys_in[1][1]
                  + mtx[2][1] * sys_in[1][2];
   p_sys_tmp[5] =   mtx[0][2] * sys_in[1][0]
                  + mtx[1][2] * sys_in[1][1]
                  + mtx[2][2] * sys_in[1][2];
   
   // Z-vector 
   p_sys_tmp[6] =   mtx[0][0] * sys_in[2][0]
                  + mtx[1][0] * sys_in[2][1]
                  + mtx[2][0] * sys_in[2][2];
   p_sys_tmp[7] =   mtx[0][1] * sys_in[2][0]
                  + mtx[1][1] * sys_in[2][1]
                  + mtx[2][1] * sys_in[2][2];
   p_sys_tmp[8] =   mtx[0][2] * sys_in[2][0]
                  + mtx[1][2] * sys_in[2][1]
                  + mtx[2][2] * sys_in[2][2];
   
    // Copy sys_tmp to sys_out if rotating in place 
   if (sys_in == sys_out) {
      memcpy(sys_out, sys_tmp, sizeof(double)*9); 
   }
}

void AnalyticGeometryTool::rotate_system( double mtx[4][4], 
                                         double sys_in[4][4],
                                         double sys_out[4][4] )
{
   double sys_tmp[4][4];
   double* p_sys_tmp;
   
   // Check to see if rotating in place 
   if (sys_in == sys_out) {
      copy_mtx( sys_in, sys_tmp );
      p_sys_tmp = (double *)sys_tmp; 
   }
   else 
      // Have p_sys_tmp point at outgoing memory location 
      p_sys_tmp = (double *)sys_out;
   
   
   // X-vector 
   p_sys_tmp[0] =   mtx[0][0] * sys_in[0][0]
                  + mtx[1][0] * sys_in[0][1]
                  + mtx[2][0] * sys_in[0][2];
   p_sys_tmp[1] =   mtx[0][1] * sys_in[0][0]
                  + mtx[1][1] * sys_in[0][1]
                  + mtx[2][1] * sys_in[0][2];
   p_sys_tmp[2] =   mtx[0][2] * sys_in[0][0]
                  + mtx[1][2] * sys_in[0][1]
                  + mtx[2][2] * sys_in[0][2];   
   
   // Y-vector 
   p_sys_tmp[4] =   mtx[0][0] * sys_in[1][0]
                  + mtx[1][0] * sys_in[1][1]
                  + mtx[2][0] * sys_in[1][2];
   p_sys_tmp[5] =   mtx[0][1] * sys_in[1][0]
                  + mtx[1][1] * sys_in[1][1]
                  + mtx[2][1] * sys_in[1][2];
   p_sys_tmp[6] =   mtx[0][2] * sys_in[1][0]
                  + mtx[1][2] * sys_in[1][1]
                  + mtx[2][2] * sys_in[1][2];
   
   // Z-vector 
   p_sys_tmp[8] =   mtx[0][0] * sys_in[2][0]
                  + mtx[1][0] * sys_in[2][1]
                  + mtx[2][0] * sys_in[2][2];
   p_sys_tmp[9] =   mtx[0][1] * sys_in[2][0]
                  + mtx[1][1] * sys_in[2][1]
                  + mtx[2][1] * sys_in[2][2];
   p_sys_tmp[10] =  mtx[0][2] * sys_in[2][0]
                  + mtx[1][2] * sys_in[2][1]
                  + mtx[2][2] * sys_in[2][2];
   
   // Maintain the origin 
   p_sys_tmp[3] = sys_in[0][3];
   p_sys_tmp[7] = sys_in[1][3];
   p_sys_tmp[11] = sys_in[2][3];
   p_sys_tmp[15] = sys_in[3][3];   
   p_sys_tmp[12] = sys_in[3][0];
   p_sys_tmp[13] = sys_in[3][1];    
   p_sys_tmp[14] = sys_in[3][2];   
   
    // Copy sys_tmp to sys_out if rotating in place 
   if (sys_in == sys_out) {
      memcpy(sys_out, sys_tmp, sizeof(double)*16);
   }
}

double AnalyticGeometryTool::det_mtx( double m[3][3] ) 
{
   double determinant;
   
   if (!m)
      return (0.0);
   
   determinant = m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
               + m[0][1]*(m[1][2]*m[2][0]-m[1][0]*m[2][2])
               + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
   
   return (determinant);
}  

void AnalyticGeometryTool::mult_mtx( double a[3][3],double b[3][3], 
                                    double d[3][3] ) 
{       
   double c[3][3];    // working buffer 
      
   if( a != 0 && b != 0 ) {   // a & b are valid 
            
      c[0][0] = (  a[0][0] * b[0][0] + a[0][1] * b[1][0] 
                 + a[0][2] * b[2][0]);
      c[1][0] = (  a[1][0] * b[0][0] + a[1][1] * b[1][0] 
                 + a[1][2] * b[2][0]);
      c[2][0] = (  a[2][0] * b[0][0] + a[2][1] * b[1][0] 
                 + a[2][2] * b[2][0]);
      
      c[0][1] = (  a[0][0] * b[0][1] + a[0][1] * b[1][1] 
                 + a[0][2] * b[2][1]);
      c[1][1] = (  a[1][0] * b[0][1] + a[1][1] * b[1][1] 
                 + a[1][2] * b[2][1]);
      c[2][1] = (  a[2][0] * b[0][1] + a[2][1] * b[1][1] 
                 + a[2][2] * b[2][1]);
      
      c[0][2] = (  a[0][0] * b[0][2] + a[0][1] * b[1][2] 
                 + a[0][2] * b[2][2]);
      c[1][2] = (  a[1][0] * b[0][2] + a[1][1] * b[1][2]
                 + a[1][2] * b[2][2]);
      c[2][2] = (  a[2][0] * b[0][2] + a[2][1] * b[1][2]
                 + a[2][2] * b[2][2]);
      
      copy_mtx(c, d);
      
   }
   else if (a) {                 // b equals 0     
      copy_mtx(a, d);	 
   }
   else if (b) {                  // a equals 0 
      copy_mtx(b, d);	 
   }
   else {                         // a & b equal 0 
      
      copy_mtx(AGT_IDENTITY_MTX_3X3, d);	 
   }
}

void AnalyticGeometryTool::mult_mtx( double a[4][4], 
                                    double b[4][4],
                                    double d[4][4] ) 
{       
   double c[4][4];    // working buffer 
      
   if( a != 0 && b != 0 ) {   // a & b are valid 
            
      c[0][0] = (  a[0][0] * b[0][0] + a[0][1] * b[1][0] 
                 + a[0][2] * b[2][0] + a[0][3] * b[3][0]);
      c[1][0] = (  a[1][0] * b[0][0] + a[1][1] * b[1][0] 
                 + a[1][2] * b[2][0] + a[1][3] * b[3][0]);
      c[2][0] = (  a[2][0] * b[0][0] + a[2][1] * b[1][0] 
                 + a[2][2] * b[2][0] + a[2][3] * b[3][0]);
      c[3][0] = (  a[3][0] * b[0][0] + a[3][1] * b[1][0] 
                 + a[3][2] * b[2][0] + a[3][3] * b[3][0]);
      
      c[0][1] = (  a[0][0] * b[0][1] + a[0][1] * b[1][1] 
                 + a[0][2] * b[2][1] + a[0][3] * b[3][1]);
      c[1][1] = (  a[1][0] * b[0][1] + a[1][1] * b[1][1] 
                 + a[1][2] * b[2][1] + a[1][3] * b[3][1]);
      c[2][1] = (  a[2][0] * b[0][1] + a[2][1] * b[1][1] 
                 + a[2][2] * b[2][1] + a[2][3] * b[3][1]);
      c[3][1] = (  a[3][0] * b[0][1] + a[3][1] * b[1][1] 
                 + a[3][2] * b[2][1] + a[3][3] * b[3][1]);
      
      c[0][2] = (  a[0][0] * b[0][2] + a[0][1] * b[1][2] 
                 + a[0][2] * b[2][2] + a[0][3] * b[3][2]);
      c[1][2] = (  a[1][0] * b[0][2] + a[1][1] * b[1][2]
                 + a[1][2] * b[2][2] + a[1][3] * b[3][2]);
      c[2][2] = (  a[2][0] * b[0][2] + a[2][1] * b[1][2]
                 + a[2][2] * b[2][2] + a[2][3] * b[3][2]);
      c[3][2] = (  a[3][0] * b[0][2] + a[3][1] * b[1][2]
                 + a[3][2] * b[2][2] + a[3][3] * b[3][2]);
      
      c[0][3] = (  a[0][0] * b[0][3] + a[0][1] * b[1][3]
                 + a[0][2] * b[2][3] + a[0][3] * b[3][3]);
      c[1][3] = (  a[1][0] * b[0][3] + a[1][1] * b[1][3]
                 + a[1][2] * b[2][3] + a[1][3] * b[3][3]);
      c[2][3] = (  a[2][0] * b[0][3] + a[2][1] * b[1][3]
                 + a[2][2] * b[2][3] + a[2][3] * b[3][3]);
      c[3][3] = (  a[3][0] * b[0][3] + a[3][1] * b[1][3]
                 + a[3][2] * b[2][3] + a[3][3] * b[3][3]);
      
      copy_mtx(c, d);	 	 
   }
   else if (a) {                 // b equals 0     
      copy_mtx(a, d);	 
   }
   else if (b) {                  // a equals 0 
      copy_mtx(b, d);	 
   }
   else {                         // a & b equal 0 
      
      copy_mtx(AGT_IDENTITY_MTX_4X4, d);	 
   }
}

CubitStatus AnalyticGeometryTool::inv_mtx_adj( double mtx[3][3],
                                              double inv_mtx[3][3] )
{
   int i,i1,i2,j,j1,j2;
   double work_mtx[3][3];
   double determinant;
   
   // Check for null input 
   if (!mtx) {
      copy_mtx(AGT_IDENTITY_MTX_3X3, inv_mtx);
      return CUBIT_SUCCESS;
   }      
   
   // Calculate determinant 
   determinant = det_mtx(mtx);
   
   // Check for singularity 
   if (dbl_eq(determinant,0.0))
      return CUBIT_FAILURE;
   
   // Get work matrix (allow inverting in place) 
   copy_mtx(mtx, work_mtx);
   
   // Inverse is adjoint matrix divided by determinant 
   for (i=1; i<4; i++) {
      
      i1 = (i % 3) + 1;  i2 = (i1 % 3) + 1;
      
      for (j=1; j<4; j++) {
         
         j1 = (j % 3) + 1;  j2 = (j1 % 3) + 1;
         
         inv_mtx[j-1][i-1] = (work_mtx[i1-1][j1-1]*work_mtx[i2-1][j2-1] - 
            work_mtx[i1-1][j2-1]*work_mtx[i2-1][j1-1]) 
            / determinant;
         
      }
   }
   return CUBIT_SUCCESS;
}

CubitStatus AnalyticGeometryTool::inv_trans_mtx( double transf[4][4], 
                                                double inv_transf[4][4] )
{
   double scale_sq;
   double inv_sq_scale;
   double tmp_transf[4][4]; // For temporary storage of incoming matrix
   double *p_tmp_transf = NULL;

   // If input transform is 0 set output to identity matrix
   if (!transf) {
      copy_mtx( AGT_IDENTITY_MTX_4X4, inv_transf );
      return CUBIT_SUCCESS;
   }
   
   // Obtain the matrix scale
   scale_sq = transf[0][0]*transf[0][0] + transf[0][1]*transf[0][1] +
              transf[0][2]*transf[0][2];
   
   // Check for singular matrix
   if (scale_sq < (.000000001 * .000000001))
      return CUBIT_FAILURE;
   
   // Need the inverse scale squared
   inv_sq_scale = 1.0 / scale_sq;   
   
   // Check to see if inverting in place
   if (transf == inv_transf) {
      copy_mtx( transf, tmp_transf );
      p_tmp_transf = (double *)tmp_transf;
   }
   else
      p_tmp_transf = (double *)transf;
   
   // The X vector
   inv_transf[0][0] = p_tmp_transf[0] * inv_sq_scale;
   inv_transf[1][0] = p_tmp_transf[1] * inv_sq_scale;
   inv_transf[2][0] = p_tmp_transf[2] * inv_sq_scale;
   
   // The Y vector
   inv_transf[0][1] = p_tmp_transf[4] * inv_sq_scale;
   inv_transf[1][1] = p_tmp_transf[5] * inv_sq_scale;
   inv_transf[2][1] = p_tmp_transf[6] * inv_sq_scale;
   
   // The Z vector 
   inv_transf[0][2] = p_tmp_transf[8] * inv_sq_scale;
   inv_transf[1][2] = p_tmp_transf[9] * inv_sq_scale;
   inv_transf[2][2] = p_tmp_transf[10] * inv_sq_scale;
   
   // Column 4 
   inv_transf[0][3] = 0.0;
   inv_transf[1][3] = 0.0;
   inv_transf[2][3] = 0.0;
   
   // The X origin 
   inv_transf[3][0] = -inv_sq_scale * (  p_tmp_transf[0] * p_tmp_transf[12]
				       + p_tmp_transf[1] * p_tmp_transf[13]
				       + p_tmp_transf[2] * p_tmp_transf[14]);
   
   // The Y origin 
   inv_transf[3][1] = -inv_sq_scale * (  p_tmp_transf[4] * p_tmp_transf[12]
				       + p_tmp_transf[5] * p_tmp_transf[13]
				       + p_tmp_transf[6] * p_tmp_transf[14]);
   
   // The Z origin 
   inv_transf[3][2] = -inv_sq_scale * (  p_tmp_transf[8] * p_tmp_transf[12]
				       + p_tmp_transf[9] * p_tmp_transf[13]
				       + p_tmp_transf[10] * p_tmp_transf[14]);
   
   // This is always one 
   inv_transf[3][3] = 1.0;
   
   return CUBIT_SUCCESS;
}

void AnalyticGeometryTool::vecs_to_mtx( double xvec[3], 
                                       double yvec[3],
                                       double zvec[3],
                                       double matrix[3][3] )
{      
   if (xvec) 
      copy_pnt(xvec, matrix[0]);
   else
      copy_pnt(AGT_IDENTITY_MTX_3X3[0], matrix[0]);
   
   if (yvec)
      copy_pnt(yvec, matrix[1]);
   else
      copy_pnt(AGT_IDENTITY_MTX_3X3[1], matrix[1]);
   
   if (zvec)
      copy_pnt(zvec, matrix[2]);
   else
      copy_pnt(AGT_IDENTITY_MTX_3X3[2], matrix[2]);
}   

void AnalyticGeometryTool::vecs_to_mtx( double xvec[3], 
                                       double yvec[3],
                                       double zvec[3], 
                                       double origin[3],
                                       double matrix[4][4] )
{      
   if (xvec) 
      copy_pnt(xvec, matrix[0]);
   else
      copy_pnt(AGT_IDENTITY_MTX_3X3[0], matrix[0]);
   
   if (yvec)
      copy_pnt(yvec, matrix[1]);
   else
      copy_pnt(AGT_IDENTITY_MTX_3X3[1], matrix[1]);
   
   if (zvec)
      copy_pnt(zvec, matrix[2]);
   else
      copy_pnt(AGT_IDENTITY_MTX_3X3[2], matrix[2]);
   
   if( origin )
      copy_pnt(origin, matrix[3]);
   else
   {
      matrix[3][0] = 0.0;
      matrix[3][1] = 0.0;
      matrix[3][2] = 0.0;
   }

   matrix[0][3] = 0.0;
   matrix[1][3] = 0.0;
   matrix[2][3] = 0.0;
   matrix[3][3] = 1.0;    
}

void AnalyticGeometryTool::print_mtx( double mtx[3][3] )
{
   PRINT_INFO( "%f %f %f\n", mtx[0][0], mtx[0][1], mtx[0][2] );
   PRINT_INFO( "%f %f %f\n", mtx[1][0], mtx[1][1], mtx[1][2] );
   PRINT_INFO( "%f %f %f\n", mtx[2][0], mtx[2][1], mtx[2][2] );
}

void AnalyticGeometryTool::print_mtx( double mtx[4][4] )
{
   PRINT_INFO( "%f %f %f %f\n", mtx[0][0], mtx[0][1], mtx[0][2], mtx[0][3] );
   PRINT_INFO( "%f %f %f %f\n", mtx[1][0], mtx[1][1], mtx[1][2], mtx[1][3] );
   PRINT_INFO( "%f %f %f %f\n", mtx[2][0], mtx[2][1], mtx[2][2], mtx[2][3] );
   PRINT_INFO( "%f %f %f %f\n", mtx[3][0], mtx[3][1], mtx[3][2], mtx[3][3] );
}

//***************************************************************************
// 3D Points
//***************************************************************************
void AnalyticGeometryTool::copy_pnt( double pnt_in[3], double pnt_out[3] )
{   
   if (pnt_in == pnt_out)
      return;
   
   if (pnt_out == NULL)
      return;
   
   if (pnt_in == NULL) {
      pnt_out[0] = 0.0;
      pnt_out[1] = 0.0;
      pnt_out[2] = 0.0;
      return;
   }
   
   // Simply copy first point into second point    
   memcpy(pnt_out, pnt_in, sizeof(double)*3);
}

void AnalyticGeometryTool::copy_pnt( double pnt_in[3], CubitVector &cubit_vec )
{
   cubit_vec.set( pnt_in );
}
   
void AnalyticGeometryTool::copy_pnt( CubitVector &cubit_vec, double pnt_out[3] )
{
   pnt_out[0] = cubit_vec.x();
   pnt_out[1] = cubit_vec.y();
   pnt_out[2] = cubit_vec.z();
}


CubitBoolean AnalyticGeometryTool::pnt_eq( double pnt1[3],double pnt2[3] )
{
   double x = pnt2[0] - pnt1[0];  // difference in the x direction 
   double y = pnt2[1] - pnt1[1];  // difference in the y direction 
   double z = pnt2[2] - pnt1[2];  // difference in the z direction 
  
   return (dbl_eq(sqrt(x*x + y*y + z*z), 0.0));
}


CubitStatus AnalyticGeometryTool::mirror_pnt( double pnt[3], 
                                             double pln_orig[3],
                                             double pln_norm[3],
                                             double m_pnt[3])
{
   double int_pnt[3],
      vec[3];
   
   // Find intersection of point and plane 
   if (int_pnt_pln(pnt, pln_orig, pln_norm, int_pnt)) {
      // If intersection is on the plane, return 
      copy_pnt(pnt, m_pnt);
      return CUBIT_FAILURE;
   }
   
   // Find vector from pnt to int_pnt 
   get_vec(pnt, int_pnt, vec);
   
   // Traverse twice the length of vec in vec direction 
   m_pnt[0] = pnt[0] + 2.0 * vec[0];
   m_pnt[1] = pnt[1] + 2.0 * vec[1];
   m_pnt[2] = pnt[2] + 2.0 * vec[2];
   
   return CUBIT_SUCCESS;   
}


CubitStatus AnalyticGeometryTool::next_pnt( double str_pnt[3], 
                                           double vec_dir[3], 
                                           double len,
                                           double new_pnt[3])
{        
   double uv[3];  // unit vector representation of vector direction 
   
   // unitize specified direction 
   if (!unit_vec(vec_dir,uv)) {
      copy_pnt(str_pnt, new_pnt);
      return CUBIT_FAILURE;
   }
   
   // determine next point in space 
   
   new_pnt[0] = str_pnt[0] + (len * uv[0]);     
   new_pnt[1] = str_pnt[1] + (len * uv[1]);     
   new_pnt[2] = str_pnt[2] + (len * uv[2]); 
   
   return CUBIT_SUCCESS;
}


//***************************************************************************
// 3D Vectors
//***************************************************************************
CubitStatus AnalyticGeometryTool::unit_vec( double vin[3], double vout[3] )
{
   double rmag; // holds magnitude of vector 
      
   // Calculate vector magnitude 
   rmag = sqrt(vin[0]*vin[0] + vin[1]*vin[1] + vin[2]*vin[2]);
   
   // check if vector has a magnitude - catch division by zero 
   
   if (dbl_eq(rmag, 0.0)) {
      if (vin != vout)
	 copy_pnt(vin, vout);
      return CUBIT_FAILURE;
   }
   
   // divide each element of the vector by the magnitude 
      
   vout[0] = vin[0] / rmag;
   vout[1] = vin[1] / rmag;
   vout[2] = vin[2] / rmag;
      
   return CUBIT_SUCCESS;
}

double AnalyticGeometryTool::dot_vec( double uval[3], double vval[3] )
{
   double dot_val;
   
   // perform dot calculation = v[0]*u[0] + v[1]*u[1] + v[1]*u[1] 
	    
   dot_val = uval[0]*vval[0] + uval[1]*vval[1] + uval[2]*vval[2];
      
   return(dot_val);   
}

void AnalyticGeometryTool::cross_vec( double uval[3], double vval[3],
                                     double cross[3] )
{      
   // determine cross product of the two vectors 
   
   cross[0] = uval[1] * vval[2] - uval[2] * vval[1];
   cross[1] = uval[2] * vval[0] - uval[0] * vval[2];
   cross[2] = uval[0] * vval[1] - uval[1] * vval[0];   
}

void AnalyticGeometryTool::cross_vec_unit( double uval[3], double vval[3],
                                          double cross[3] )
{      
   // determine cross product of the two vectors 
   cross_vec(uval,vval,cross);
   
   // convert to unit vector 
   unit_vec(cross,cross);
}

void AnalyticGeometryTool::orth_vecs( double uvect[3], double vvect[3], 
                                     double wvect[3] )
{
   double x[3];
   unsigned short i = 0;
   unsigned short imin = 3;
   double rmin = 1.0E20;
   unsigned short iperm1[3];
   unsigned short iperm2[3];
   unsigned short cont_flag = 1;
   double vec[3];
   
   // Initialize perm flags
   iperm1[0] = 1; iperm1[1] = 2; iperm1[2] = 0;
   iperm2[0] = 2; iperm2[1] = 0; iperm2[2] = 1;
   
   // unitize vector 
   
   unit_vec(uvect,vec);
   
   while (i<3 && cont_flag ) {
      if (dbl_eq(vec[i], 0.0)) {
         vvect[i] = 1.0;
         vvect[iperm1[i]] = 0.0;
         vvect[iperm2[i]] = 0.0;
         cont_flag = 0;
      }
      
      if (fabs(vec[i]) < rmin) {
         imin = i;
         rmin = fabs(vec[i]);
      }
      ++i;
   }
   
   if (cont_flag) {
      x[imin] = 1.0;
      x[iperm1[imin]] = 0.0;
      x[iperm2[imin]] = 0.0;
      
      // determine cross product 
      cross_vec_unit(vec,x,vvect);
   }
   
   // cross vector to determine last orthogonal vector 
   cross_vec(vvect,vec,wvect);
}

double AnalyticGeometryTool::mag_vec( double vec[3] )
{
   return (sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]));
}

CubitStatus AnalyticGeometryTool::get_vec ( double str_pnt[3], 
                                           double stp_pnt[3],
                                           double vector_out[3] )
{  
   // Make sure we can create a vector 
   if (pnt_eq(str_pnt, stp_pnt)) {
      vector_out[0] = 0.0;
      vector_out[1] = 0.0;
      vector_out[2] = 0.0;
      return CUBIT_FAILURE;
   }
   
   // determine vector by subtracting starting point from stopping point 
   
  vector_out[0] = stp_pnt[0] - str_pnt[0];
  vector_out[1] = stp_pnt[1] - str_pnt[1];
  vector_out[2] = stp_pnt[2] - str_pnt[2];  
  
  return CUBIT_SUCCESS;
}

CubitStatus AnalyticGeometryTool::get_vec_unit( double str_pnt[3], 
                                               double stp_pnt[3],
                                               double uv_out[3] )
{
  // determine vector between points 
   if (!get_vec(str_pnt,stp_pnt,uv_out))
      return CUBIT_FAILURE;
  
  // unitize vector 
  if (!unit_vec(uv_out,uv_out))
     return CUBIT_FAILURE;
     
  return CUBIT_SUCCESS;  
}

void AnalyticGeometryTool::mult_vecxconst( double constant,
                                          double vec[3],
                                          double vec_out[3] )
{
   // multiply each element of the vector by the constant
   vec_out[0] = constant * vec[0];
   vec_out[1] = constant * vec[1];
   vec_out[2] = constant * vec[2];   
}


void AnalyticGeometryTool::reverse_vec( double vin[3],double vout[3] )
{
   // Multiply the vector components by a -1.0    
   mult_vecxconst(-1.0, vin, vout);      
}

double AnalyticGeometryTool::angle_vec_vec( double v1[3],double v2[3] )
{   
  double denom, dot, cosang, sinang, acrsb, angle;
  double crossed_vec[3];
     
  // For accuracy, use cosine for large angles, sine for small angles 
  denom = sqrt((v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2])
              *(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]));
  
  // Check for a zero length vector 
  if (dbl_eq(denom, 0.0))
     return (0.0);
  
  dot = dot_vec(v1, v2);
  
  cosang = dot/denom;
  
  if (1.0 - fabs(cosang) < 0.01) {
     cross_vec(v1, v2, crossed_vec);
     acrsb = mag_vec(crossed_vec);
     sinang = acrsb/denom;
     if (cosang > 0.0) 
        angle = asin(sinang);
     else 
        angle = AGT_PI - asin(sinang);
  }
  else 
     angle = acos(cosang);
  
  return (angle); 
}

//***************************************************************************
// Distances
//***************************************************************************
double AnalyticGeometryTool::dist_pnt_pnt( double pnt1[3], double pnt2[3] )
{  
  double x = pnt2[0] - pnt1[0];  // difference in the x direction
  double y = pnt2[1] - pnt1[1];  // difference in the y direction
  double z = pnt2[2] - pnt1[2];  // difference in the z direction
  
  // return the distance   
  return(sqrt(x*x + y*y + z*z));
}



double AnalyticGeometryTool::dist_pln_pln( double pln_1_orig[3],
                                          double pln_1_norm[3], 
                                          double pln_2_orig[3],
                                          double pln_2_norm[3],
                                          AgtSide *side,
                                          AgtOrientation *orien,
                                          unsigned short *status )   
{
   double distance;
   double int_pnt[3];
   double vec[3];
   
   // Check to see if planes are parallel
   if (is_vec_par(pln_1_norm, pln_2_norm)) {    
      
      // Set successful status
      if (status)
         *status = 1;
      
      // Calculate perpendicular line plane intersection on plane_2 from 
      // pln_1_origin 
      int_ln_pln(pln_1_orig, pln_1_norm, pln_2_orig, pln_2_norm, 
                 int_pnt);
      
      // Find distance between pln_1_origin and this intersection pnt
      distance = dist_pnt_pnt(pln_1_orig, int_pnt);
      
      // Get side if required
      if (side) {
         
         if (dbl_eq(distance, 0.0))
            
            *side = AGT_ON_PLANE;
         
         else {
            
            // Get vector to intersection point
            get_vec(pln_1_orig, int_pnt, vec);
            
            // Compare angles
            if (dbl_eq(angle_vec_vec(vec, pln_1_norm), 0.0))
               *side = AGT_POS_SIDE;
            else
               *side = AGT_NEG_SIDE;
            
         }
      }
      
      // Get orientation if required
      if (orien) {
         
         // Compare surface normals
         if (dbl_eq(angle_vec_vec(pln_1_norm, pln_2_norm), 0.0))
            *orien = AGT_SAME_DIR;
         else
            *orien = AGT_OPP_DIR;
      }      
      
   }
   
   else {
      
      if (status)
         *status = 0;
      
      if (side)
         *side = AGT_CROSS;
      
      distance = 0.0;
      
   }
   
   return (distance);
}



//***************************************************************************
// Intersections
//***************************************************************************
CubitStatus AnalyticGeometryTool::int_ln_pln( double ln_orig[3],
                                             double ln_vec[3],
                                             double pln_orig[3],
                                             double pln_norm[3],
                                             double int_pnt[3] )
{   
   double denom;
   double t;
   
   // Set parametric eqns of line equal to parametric eqn of plane & solve
   // for t 
   denom = pln_norm[0]*ln_vec[0] + pln_norm[1]*ln_vec[1] + 
           pln_norm[2]*ln_vec[2];
   
   if (dbl_eq(denom, 0.0))
      return CUBIT_FAILURE;
   
   t = (pln_norm[0]*(pln_orig[0]-ln_orig[0]) + 
      pln_norm[1]*(pln_orig[1]-ln_orig[1]) + 
      pln_norm[2]*(pln_orig[2]-ln_orig[2]))/denom;
   
   // Substitute t back into equations of line to get xyz 
   int_pnt[0] = ln_orig[0] + ln_vec[0]*t;
   int_pnt[1] = ln_orig[1] + ln_vec[1]*t;
   int_pnt[2] = ln_orig[2] + ln_vec[2]*t;
      
   return CUBIT_SUCCESS;   
}

int AnalyticGeometryTool::int_ln_ln( double p1[3], double v1[3],
                                    double p2[3], double v2[3], 
                                    double int_pnt1[3], double int_pnt2[3] )
{   
   double norm[3];
   double pln1_norm[3];
   double pln2_norm[3];
   
   // Cross the two vectors to get a normal vector 
   cross_vec(v1, v2, norm);
   
   if (dbl_eq(mag_vec(norm), 0.0))
      return 0;
   
   // Cross v2 & normal to get normal to plane 2 
   cross_vec(v2, norm, pln2_norm);
   
   // Find intersection of pln2_norm and vector 1 - this is int_pnt1 
   int_ln_pln(p1, v1, p2, pln2_norm, int_pnt1);
   
   // Cross v1 & normal to get normal to plane 1 
   cross_vec(v1, norm, pln1_norm);
   
   // Find intersection of pln2_norm and vector 1 - this is int_pnt2 
   int_ln_pln(p2, v2, p1, pln1_norm, int_pnt2);

   // Check to see if the intersection points are the same & return 
   if (pnt_eq(int_pnt1, int_pnt2)) 
      return 1;
   else
      return 2;   
}


int AnalyticGeometryTool::int_pnt_pln( double pnt[3],
                                      double pln_orig[3], 
                                      double pln_norm[3], 
                                      double pln_int[3] )   
{
   // Calculate line plane intersection w/plane normal as line vector 
   int_ln_pln(pnt, pln_norm, pln_orig, pln_norm, pln_int);
   
   // Check to see if point is on the plane 
   if (pnt_eq(pln_int, pnt))
      return 1;
   else
      return 0; 
}



//***************************************************************************
// Comparison/Containment Tests
//***************************************************************************
CubitBoolean AnalyticGeometryTool::is_vec_par( double vec_1[3],
                                              double vec_2[3] )
{   
   double cross[3];
   
   // Get cross product & see if its magnitude is zero
   cross_vec(vec_1, vec_2, cross);
   
   if (dbl_eq(mag_vec(cross), 0.0))
      return CUBIT_TRUE;
   else
      return CUBIT_FALSE;
}

CubitBoolean AnalyticGeometryTool::is_vec_perp( double vec_1[3],double vec_2[3])
{  
   // Check angle between vectors
   if (dbl_eq(angle_vec_vec(vec_1, vec_2), AGT_PI_DIV_2))
      return CUBIT_TRUE;
   else
      return CUBIT_FALSE;   
}

CubitBoolean AnalyticGeometryTool::is_vecs_same_dir( double vec_1[3], 
                                                     double vec_2[3] )
{         
   // Check to see if angle between vectors can be considered zero   
   if (dbl_eq(angle_vec_vec(vec_1, vec_2), 0.0))
      return CUBIT_TRUE;
   else
      return CUBIT_FALSE;   
}


CubitBoolean AnalyticGeometryTool::is_pnt_on_ln( double pnt[3],
                                                double ln_orig[3],
                                                double ln_vec[3] )
{  
   double vec[3];
   
   // Compare pnt and line origin
   if (pnt_eq(pnt, ln_orig))
      return CUBIT_TRUE;   
   
   // Get a vector from line origin to the point      
   get_vec(ln_orig, pnt, vec);
   
   // If this vector is parallel with line vector, point is on the line
   if (is_vec_par(vec, ln_vec))
      return CUBIT_TRUE;
   else
      return CUBIT_FALSE;
}

CubitBoolean AnalyticGeometryTool::is_pnt_on_ln_seg( double pnt[3],
                                                    double end1[3],
                                                    double end2[3] )
{ 
//   METHOD:
//    o Use parametric equations of line
//    
//         x = x1 + t(x2 - x1)
//         y = y1 + t(y2 - y1)
//         z = z1 + t(z2 - z1)
//	 
//    o Note:  two other method's were considered: 
//       1) Comparing sum of distance of point to both end points to the
//          line length.
//       2) Checking to see if area of a triangle with the vertices is zero
//       
//     Using parametric equations is more efficient in many cases.
  double t1 = 0.0,
    t2 = 0.0,
    t3 = 0.0,
    neg_range,
    pos_range;
  
  unsigned short flg1 = 0,
    flg2 = 0,
    flg3 = 0;
  
  neg_range = 0.0 - agtEpsilon;
  pos_range = 1.0 + agtEpsilon;
  
  if (fabs(end2[0] - end1[0]) < agtEpsilon)
  {
    if (fabs(pnt[0] - end1[0]) < agtEpsilon)
      flg1 = 1;
    else
      return CUBIT_FALSE;
  }
  else
  {
    t1 = (pnt[0] - end1[0])/(end2[0] - end1[0]);
    
    if (t1<neg_range || t1>pos_range)
      return CUBIT_FALSE;
  }
  
  if (fabs(end2[1] - end1[1]) < agtEpsilon)
  {
    if (fabs(pnt[1] - end1[1]) < agtEpsilon)
      flg2 = 1;
    else
      return CUBIT_FALSE;
  }
  else
  {
    t2 = (pnt[1] - end1[1])/(end2[1] - end1[1]);
    
    if (t2<neg_range || t2>pos_range)
      return CUBIT_FALSE;
  }
   
  if (fabs(end2[2] - end1[2]) < agtEpsilon)
  {
    if (fabs(pnt[2] - end1[2]) < agtEpsilon)
      flg3 = 1;
    else
      return CUBIT_FALSE;
  }
  else
  {
    t3 = (pnt[2] - end1[2])/(end2[2] - end1[2]);
    
    if (t3<neg_range || t3>pos_range)
      return CUBIT_FALSE;
  }
  
    // If any 2 flags are 1, point is on the line,
    // otherwise, check remaining T's for equality
  
  if (flg1)
  {
      // Here, flg1 = 1
    
    if (flg2)
    {
        // Here, flg1 = 1
        // Here, flg2 = 1
      return CUBIT_TRUE;
    }
    else
    {
        // Here, flg1 = 1
        // Here, flg2 = 0
      
      if (flg3)
          // Here, flg1 = 1
          // Here, flg2 = 0
          // Here, flg3 = 1
        return CUBIT_TRUE;
      else
      {
          // Here, flg1 = 1
          // Here, flg2 = 0
          // Here, flg3 = 0
        if (dbl_eq(t2, t3))
          return CUBIT_TRUE;
        else
          return CUBIT_FALSE;
      }
    }
  }
  else
  {
      // Here, flg1 = 0
    if (flg2)
    {
        // Here, flg1 = 0
        // Here, flg2 = 1
      if (flg3)
        return CUBIT_TRUE;
          // Here, flg1 = 0
          // Here, flg2 = 1
          // Here, flg3 = 1
      else
      {
          // Here, flg1 = 0
          // Here, flg2 = 1
          // Here, flg3 = 0
        if (dbl_eq(t1, t3))
          return CUBIT_TRUE;
        else
          return CUBIT_FALSE;
      }
    }
    else
    {
        // Here, flg1 = 0
        // Here, flg2 = 0
      if (flg3)
      {
          // Here, flg1 = 0
          // Here, flg2 = 0
          // Here, flg3 = 1
        if (dbl_eq(t1, t2))
          return CUBIT_TRUE;
        else
          return CUBIT_FALSE;
      }
      else
      {
          // Here, flg1 = 0
          // Here, flg2 = 0
          // Here, flg3 = 0
        if (dbl_eq(t1, t2) && dbl_eq(t1, t3))
          return CUBIT_TRUE;
        else
          return CUBIT_FALSE;
      }
    }
  }
  
    // This would be a programmer's error if we got to this point
//  return CUBIT_FALSE;
}

CubitBoolean AnalyticGeometryTool::is_pnt_on_pln( double pnt[3],
                                                 double pln_orig[3],
                                                 double pln_norm[3] )
{  
   double result;
   
   // See if point satisfies parametric equation of plane
   
   result = pln_norm[0] * (pnt[0] - pln_orig[0]) +
            pln_norm[1] * (pnt[1] - pln_orig[1]) +
            pln_norm[2] * (pnt[2] - pln_orig[2]);
   
   if (dbl_eq(result, 0.0))
      return CUBIT_TRUE;
   else
      return CUBIT_FALSE;
}

CubitBoolean AnalyticGeometryTool::is_ln_on_pln( double ln_orig[3],
                                                double ln_vec[3],
                                                double pln_orig[3],
                                                double pln_norm[3] )
{     
   
   // Check to see if line origin is on the plane.  If so, check to see if
   // line vector is perpendicular to plane normal.
   
   if (is_pnt_on_pln(ln_orig, pln_orig, pln_norm) &&
       is_vec_perp(ln_vec, pln_norm))      
      return CUBIT_TRUE;      
   else      
      return CUBIT_FALSE;
}



CubitBoolean AnalyticGeometryTool::is_pln_on_pln( double pln_orig1[3],
                                                 double pln_norm1[3],
                                                 double pln_orig2[3],
                                                 double pln_norm2[3] )
{  
   // If 1st plane origin is on the 2nd plane and the normals are
   // parallel they are coincident
   if( is_vec_par( pln_norm1, pln_norm2 ) &&
      is_pnt_on_pln( pln_orig1, pln_orig2, pln_norm2 ) )
      return CUBIT_TRUE;
   else
      return CUBIT_FALSE;
}

//***************************************************************************
// Arcs/Circles
//***************************************************************************
void AnalyticGeometryTool::setup_arc( double vec_1[3], double vec_2[3], 
                                     double origin[3], double start_angle, 
                                     double end_angle, double radius,
                                     AGT_3D_Arc &arc )
{
   copy_pnt( vec_1, arc.Vec1 );
   copy_pnt( vec_2, arc.Vec2 );
   copy_pnt( origin, arc.Origin );
   arc.StartAngle = start_angle;
   arc.EndAngle = end_angle;
   arc.Radius = radius;
}

void AnalyticGeometryTool::setup_arc( CubitVector& vec_1, CubitVector& vec_2, 
                                     CubitVector origin, double start_angle, 
                                     double end_angle, double radius,
                                     AGT_3D_Arc &arc )
{
   vec_1.get_xyz( arc.Vec1 );
   vec_2.get_xyz( arc.Vec2 );
   origin.get_xyz( arc.Origin );
   arc.StartAngle = start_angle;
   arc.EndAngle = end_angle;
   arc.Radius = radius;
}

void AnalyticGeometryTool::get_arc_xyz( AGT_3D_Arc &arc, double param, double pnt[3] )
{
   double Tp;
   
   // Un-normalized parameter
   Tp = arc.StartAngle * ( 1.0 - param ) + arc.EndAngle * param;
   
   // Solve for XYZ
   pnt[0] = arc.Radius * ( cos( Tp ) * arc.Vec1[0] +
                           sin( Tp ) * arc.Vec2[0] ) + 
                           arc.Origin[0];
   
   pnt[1] = arc.Radius * ( cos( Tp ) * arc.Vec1[1] +
                           sin( Tp ) * arc.Vec2[1] ) + 
                           arc.Origin[1];
   
   pnt[2] = arc.Radius * ( cos( Tp ) * arc.Vec1[2] +
                           sin( Tp ) * arc.Vec2[2] ) + 
                           arc.Origin[2];
}

void AnalyticGeometryTool::get_arc_xyz( AGT_3D_Arc &arc, double param, CubitVector& pnt )
{
   double Tp;
   
   // Un-normalized parameter
   Tp = arc.StartAngle * ( 1.0 - param ) + arc.EndAngle * param;
   
   // Solve for XYZ
   pnt.x( arc.Radius * ( cos( Tp ) * arc.Vec1[0] +
                         sin( Tp ) * arc.Vec2[0] ) + 
                         arc.Origin[0] );
   
   pnt.y( arc.Radius * ( cos( Tp ) * arc.Vec1[1] +
                         sin( Tp ) * arc.Vec2[1] ) + 
                         arc.Origin[1] );
   
   pnt.z( arc.Radius * ( cos( Tp ) * arc.Vec1[2] +
                         sin( Tp ) * arc.Vec2[2] ) + 
                         arc.Origin[2] );
}

int 
AnalyticGeometryTool::get_num_circle_tess_pnts( double radius, 
                                                double len_tol )
{
  double cmin, cmax;

  double c = 2*CUBIT_PI*radius; // Circumference

  // Find the number of points required for the given accuracy.  Use
  // a bisection method.
  int nmin = 8, nmax = 100;
  cmin = 2.0*nmin*radius*sin(CUBIT_PI/nmin); // Circumference of circle using segments
  cmax = 2.0*nmax*radius*sin(CUBIT_PI/nmax);

  if( dbl_eq( cmin, c ) )
    return nmin;

  double old_epsilon = set_epsilon( len_tol );

  // Find an n that is more than accurate enough
  while( !dbl_eq( cmax, c ) )
  {
    nmin = nmax;
    nmax = nmin * 10;
    cmin = 2.0*nmin*radius*sin(CUBIT_PI/nmin);
    cmax = 2.0*nmax*radius*sin(CUBIT_PI/nmax);
  }
  
  // Biscect until the minimum number of segments satisfying
  // the tolerance is found.
  int n;
  while( 1 )
  {
    n = (nmin + nmax)/2;
    double cn = 2.0*n*radius*sin(CUBIT_PI/n);
    if( dbl_eq( cn, c ) )
    {
      // Go lower
      nmax = n;
    }
    else
    {
      // Go higher
      nmin = n;
    }
    if( nmax-nmin < 2 )
      break;
  }
  set_epsilon( old_epsilon );

  return nmax;
}

//***************************************************************************
// Miscellaneous
//***************************************************************************
void AnalyticGeometryTool::get_pln_orig_norm( double A, double B, double C,
                                             double D, double pln_orig[3], 
                                             double pln_norm[3] )
{
   double x = 0.0, y = 0.0, z = 0.0;

   // Try to have origin aligned with one of the principal axes
   if( !dbl_eq( C, 0.0 ) )
      z = -D/C;
   else if (!dbl_eq( A, 0.0 ) )
      x = -D/A;
   else if (!dbl_eq( B, 0.0 ) )
      y = -D/B;

   pln_orig[0] = x;
   pln_orig[1] = y;
   pln_orig[2] = z;

   if( pln_norm )
   {
      pln_norm[0] = A;
      pln_norm[1] = B;
      pln_norm[2] = C;
   }
}

void AnalyticGeometryTool::get_box_corners( double box_min[3], 
                                           double box_max[3], 
                                           double c[8][3] )
{
   // Left-Bottom-Front       // Left-Top-Front
   c[0][0] = box_min[0];      c[1][0] = box_min[0];
   c[0][1] = box_min[1];      c[1][1] = box_max[1];
   c[0][2] = box_max[2];      c[1][2] = box_max[2];

   // Right-Top-Front         // Right-Bottom-Front
   c[2][0] = box_max[0];      c[3][0] = box_max[0];
   c[2][1] = box_max[1];      c[3][1] = box_min[1];
   c[2][2] = box_max[2];      c[3][2] = box_max[2];

   // Left-Bottom-Back        // Left-Top-Back
   c[4][0] = box_min[0];      c[5][0] = box_min[0];
   c[4][1] = box_min[1];      c[5][1] = box_max[1];
   c[4][2] = box_min[2];      c[5][2] = box_min[2];

   // Right-Top-Back          // Right-Bottom-Back
   c[6][0] = box_max[0];      c[7][0] = box_max[0];
   c[6][1] = box_max[1];      c[7][1] = box_min[1];
   c[6][2] = box_min[2];      c[7][2] = box_min[2];
   
}

CubitStatus 
AnalyticGeometryTool::min_pln_box_int_corners( const CubitPlane& plane,
                                              const CubitBox& box,
                                              int extension_type,
                                              double extension,
                                              CubitVector& p1, CubitVector& p2,
                                              CubitVector& p3, CubitVector& p4,
                                              CubitBoolean silent )
{
   CubitVector box_min = box.minimum();
   CubitVector box_max = box.maximum();
   
   CubitVector plane_norm = plane.normal();
   
   double box_min_pnt[3], box_max_pnt[3], pln_norm[3];
   box_min.get_xyz( box_min_pnt ); box_max.get_xyz( box_max_pnt );
   plane_norm.get_xyz( pln_norm );
   
   double pnt1[3], pnt2[3], pnt3[3], pnt4[3];
   
   if( min_pln_box_int_corners( pln_norm, plane.coefficient(),
      box_min_pnt, box_max_pnt, 
      extension_type, extension, 
      pnt1, pnt2, pnt3, pnt4, silent ) == CUBIT_FAILURE )
         return CUBIT_FAILURE;
   
   p1.set( pnt1 ); p2.set( pnt2 ); p3.set( pnt3 ); p4.set( pnt4 );
   
   return CUBIT_SUCCESS;
}

CubitStatus 
AnalyticGeometryTool::min_pln_box_int_corners( CubitVector& vec1,
                                              CubitVector& vec2,
                                              CubitVector& vec3,
                                              CubitVector& box_min,
                                              CubitVector& box_max,
                                              int extension_type,
                                              double extension,
                                              CubitVector& p1, CubitVector& p2, 
                                              CubitVector& p3, CubitVector& p4,
                                              CubitBoolean silent )
{
   CubitPlane plane;
   if( plane.mk_plane_with_points( vec1, vec2, vec3 ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

   CubitVector plane_norm = plane.normal();
   double coefficient = plane.coefficient();

   double plane_norm3[3];
   double box_min3[3];
   double box_max3[3];

   box_min.get_xyz( box_min3 );
   box_max.get_xyz( box_max3 );
   plane_norm.get_xyz( plane_norm3 );

   double p1_3[3], p2_3[3], p3_3[3], p4_3[3];
   p1.get_xyz( p1_3 );
   p2.get_xyz( p2_3 );
   p3.get_xyz( p3_3 );
   p4.get_xyz( p4_3 );
   
   if( min_pln_box_int_corners( plane_norm3, coefficient, box_min3, box_max3, 
      extension_type, extension, p1_3, p2_3, p3_3, p4_3, silent ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

   p1.set( p1_3 );
   p2.set( p2_3 );
   p3.set( p3_3 );
   p4.set( p4_3 );

   return CUBIT_SUCCESS;
}

CubitStatus 
AnalyticGeometryTool::min_pln_box_int_corners( double pln_norm[3], 
                                              double pln_coeff,
                                              double box_min[3],
                                              double box_max[3],
                                              int extension_type,
                                              double extension,
                                              double p1[3], double p2[3],
                                              double p3[3], double p4[3],
                                              CubitBoolean silent )
{
   int i;
   double cubit2pln_mtx[4][4],
      pln2cubit_mtx[4][4];
   double pln_orig[3];

   double A = pln_norm[0];
   double B = pln_norm[1];
   double C = pln_norm[2];
   double D = pln_coeff;

   //PRINT_INFO( "A=%0.4lf, B=%0.4lf, C=%0.4lf, D=%0.4lf\n", A, B, C, D );

   get_pln_orig_norm( A, B, C, D, pln_orig );

//   PRINT_INFO( "Plane Orig = %0.4lf, %0.4lf, %0.4lf\n", pln_orig[0],
//      pln_orig[1], pln_orig[2] );

   // Find intersections of edges with plane.  Add to unique
   // array.  At most there are 6 intersections...
   double int_array[6][3];
   int num_int = 0;
   num_int = get_plane_bbox_intersections( box_min, box_max, pln_orig, pln_norm, int_array );

   //attempt to adjust bounding box to x,y,z intercepts of plane
   if( num_int == 0 )
   {
     //Stretch bounding box so that plane will fit, for sure
     //get x,y,z intercepts
     double x_intercept = 0;
     double y_intercept = 0;
     double z_intercept = 0;
     if( !dbl_eq( A, 0.0 ) )
       x_intercept = -D/A;
     if( !dbl_eq( B, 0.0 ) )
       y_intercept = -D/B;
     if( !dbl_eq( C, 0.0 ) )
       z_intercept = -D/C;
    
    //adjust box 
    if( x_intercept < box_min[0] )
      box_min[0] = x_intercept;
    else if( x_intercept > box_max[0] )
      box_max[0] = x_intercept;

    if( y_intercept < box_min[1] )
      box_min[1] = y_intercept;
    else if( y_intercept > box_max[1] )
      box_max[1] = y_intercept;
      
    if( z_intercept < box_min[2] )
      box_min[2] = z_intercept;
    else if( z_intercept > box_max[2] )
      box_max[2] = z_intercept;

     num_int = get_plane_bbox_intersections( box_min, box_max, pln_orig, pln_norm, int_array );
   }

   if( num_int == 0 )
   {
     if( silent == CUBIT_FALSE )
       PRINT_ERROR( "Plane does not intersect the bounding box\n"
       "       Can't find 4 corners of plane\n" );
     return CUBIT_FAILURE;
   }
   if( num_int < 3 )
   {
     if( silent == CUBIT_FALSE )
       PRINT_ERROR( "Plane intersects the bounding box at only %d locations\n"
       "      Can't calculate 4 corners of plane\n", num_int );
      return CUBIT_FAILURE;
   }

   // Transform pnts to plane coordinate system
   double pln_x[3], pln_y[3];
   orth_vecs( pln_norm, pln_x, pln_y );
   vecs_to_mtx( pln_x, pln_y, pln_norm, pln_orig, pln2cubit_mtx );
   inv_trans_mtx( pln2cubit_mtx, cubit2pln_mtx );

   double int_arr_pln[6][3];
   for( i=0; i<num_int; i++ )
      transform_pnt( cubit2pln_mtx, int_array[i], int_arr_pln[i] );

   // Place into format for mimimal box calculation
   Point2 pt[6];
   for ( i=0; i<num_int; i++ )
   {
      pt[i].x = int_arr_pln[i][0];
      pt[i].y = int_arr_pln[i][1];
      if( !dbl_eq( int_arr_pln[i][2], 0.0 ) )
      {
        if( silent == CUBIT_FALSE )
          PRINT_ERROR( "in AnalyticGeometryTool::min_box_pln_int_corners\n"
          "       Transform to plane wrong\n" );
         return CUBIT_FAILURE;
      }
   }
   
   // Find rectangle with minimal area to surround the points
   // (this is definitely overkill esp. for the simple cases.....)
   OBBox2 minimal = MinimalBox2( num_int, pt );
   
   // Strip out results
   double old_epsilon = set_epsilon( 1e-10 );
   double centroid[3];
   centroid[0] = minimal.center.x;
   centroid[1] = minimal.center.y;
   centroid[2] = 0.0;
   round_near_val( centroid[0] ); // Makes near -1, 0, 1 values -1, 0, 1
   round_near_val( centroid[1] );
   transform_pnt( pln2cubit_mtx, centroid, centroid );
   
   double x_axis[3];
   x_axis[0] = minimal.axis[0].x;
   x_axis[1] = minimal.axis[0].y;
   x_axis[2] = 0.0;
   round_near_val( x_axis[0] );
   round_near_val( x_axis[1] );
   transform_vec( pln2cubit_mtx, x_axis, x_axis );
   
   double y_axis[3];
   y_axis[0] = minimal.axis[1].x;
   y_axis[1] = minimal.axis[1].y;
   y_axis[2] = 0.0;
   round_near_val( y_axis[0] );
   round_near_val( y_axis[1] );
   transform_vec( pln2cubit_mtx, y_axis, y_axis );

   set_epsilon( old_epsilon );
   
   double dist_x;
   double dist_y;
   double extension_distance = 0.0;
   if( extension_type == 1 ) // Percentage (of 1/2 diagonal)
   {
      double diag_len = sqrt(  minimal.extent[0]*minimal.extent[0]
                             + minimal.extent[1]*minimal.extent[1] );
      extension_distance = diag_len*extension/100.0;
   }
   else if( extension_type == 2 ) // Absolute distance in x and y
      extension_distance = extension;

   dist_x = minimal.extent[0] + extension_distance;
   dist_y = minimal.extent[1] + extension_distance;
   
   next_pnt( centroid, x_axis, -dist_x, p1 );
   next_pnt( p1, y_axis, -dist_y, p1 );
   
   next_pnt( centroid, x_axis, -dist_x, p2 );
   next_pnt( p2, y_axis, dist_y, p2 );
   
   next_pnt( centroid, x_axis, dist_x, p3 );
   next_pnt( p3, y_axis, dist_y, p3 );
   
   next_pnt( centroid, x_axis, dist_x, p4 );
   next_pnt( p4, y_axis, -dist_y, p4 );
   
   return CUBIT_SUCCESS;
}

int AnalyticGeometryTool::get_plane_bbox_intersections( double box_min[3],
                                                           double box_max[3],
                                                           double pln_orig[3],
                                                           double pln_norm[3],
                                                           double int_array[6][3])
{

   // Fill in an array with all 8 box corners
   double corner[8][3];
   get_box_corners( box_min, box_max, corner );
 
   // Get 12 edges of the box
   double ln_start[12][3], ln_end[12][3];
   copy_pnt( corner[0], ln_start[0] );  copy_pnt( corner[1], ln_end[0] );
   copy_pnt( corner[1], ln_start[1] );  copy_pnt( corner[2], ln_end[1] );
   copy_pnt( corner[2], ln_start[2] );  copy_pnt( corner[3], ln_end[2] );
   copy_pnt( corner[3], ln_start[3] );  copy_pnt( corner[0], ln_end[3] ); 
   copy_pnt( corner[4], ln_start[4] );  copy_pnt( corner[5], ln_end[4] );
   copy_pnt( corner[5], ln_start[5] );  copy_pnt( corner[6], ln_end[5] );
   copy_pnt( corner[6], ln_start[6] );  copy_pnt( corner[7], ln_end[6] );
   copy_pnt( corner[7], ln_start[7] );  copy_pnt( corner[4], ln_end[7] );
   copy_pnt( corner[0], ln_start[8] );  copy_pnt( corner[4], ln_end[8] );
   copy_pnt( corner[1], ln_start[9] );  copy_pnt( corner[5], ln_end[9] );
   copy_pnt( corner[2], ln_start[10] ); copy_pnt( corner[6], ln_end[10] );
   copy_pnt( corner[3], ln_start[11] ); copy_pnt( corner[7], ln_end[11] );
   
   double ln_vec[3];
   double int_pnt[3];
   int num_int = 0;
   int i, j, found;
   for( i=0; i<12; i++ )
   {
      get_vec_unit( ln_start[i], ln_end[i], ln_vec );
      if( int_ln_pln( ln_start[i], ln_vec, pln_orig, pln_norm, int_pnt ) )
      {
         // Only add if on the bounded line segment
         if( is_pnt_on_ln_seg( int_pnt, ln_start[i], ln_end[i] ) )
         {
            // Only add if unique
            found = 0;
            for( j=0; j<num_int; j++ )
            {
               if( pnt_eq( int_pnt, int_array[j] ) )
               {
                  found = 1;
                  break;
               }
            }
            if( !found )
            {
               copy_pnt( int_pnt, int_array[num_int] );
               num_int++;
            }
         }
      }
   }
  return num_int;
}

CubitStatus
AnalyticGeometryTool::get_tight_bounding_box( DLIList<CubitVector*> &point_list,                                              
                                              CubitVector &center,
                                              CubitVector axes[3],
                                              CubitVector &extension )
{
   int num_pnts = point_list.size();
   if( num_pnts == 0 )
      return CUBIT_FAILURE;
   Point3 *pnt_arr = new Point3[num_pnts];

   int i;
   point_list.reset();
   CubitVector *vec;
   for( i=0; i<num_pnts; i++ )
   {
      vec = point_list.get_and_step();

      pnt_arr[i].x = vec->x();
      pnt_arr[i].y = vec->y();
      pnt_arr[i].z = vec->z();
   }

   OBBox3 minimal = MinimalBox3 (point_list.size(), pnt_arr);

   //PRINT_INFO( "MinBox center = %lf, %lf, %lf\n", minimal.center.x, minimal.center.y, minimal.center.z );
   //PRINT_INFO( "MinBox axis 1 = %lf, %lf, %lf\n", minimal.axis[0].x, minimal.axis[0].y, minimal.axis[0].z );
   //PRINT_INFO( "MinBox axis 2 = %lf, %lf, %lf\n", minimal.axis[1].x, minimal.axis[1].y, minimal.axis[1].z );
   //PRINT_INFO( "MinBox axis 3 = %lf, %lf, %lf\n", minimal.axis[2].x, minimal.axis[2].y, minimal.axis[2].z );
   //PRINT_INFO( "MinBox extent = %lf, %lf, %lf\n", minimal.extent[0], minimal.extent[1], minimal.extent[2] );

   center.set(minimal.center.x, minimal.center.y, minimal.center.z);
   axes[0].set(minimal.axis[0].x, minimal.axis[0].y, minimal.axis[0].z);
   axes[1].set(minimal.axis[1].x, minimal.axis[1].y, minimal.axis[1].z);
   axes[2].set(minimal.axis[2].x, minimal.axis[2].y, minimal.axis[2].z);
   extension.set(minimal.extent[0], minimal.extent[1], minimal.extent[2]);

   delete [] pnt_arr;

   return CUBIT_SUCCESS;
}

CubitStatus
AnalyticGeometryTool::get_tight_bounding_box( DLIList<CubitVector*> &point_list,                                              
                                              CubitVector& p1, CubitVector& p2,
                                              CubitVector& p3, CubitVector& p4,
                                              CubitVector& p5, CubitVector& p6,
                                              CubitVector& p7, CubitVector& p8)
{
   int num_pnts = point_list.size();
   if( num_pnts == 0 )
      return CUBIT_FAILURE;
   Point3 *pnt_arr = new Point3[num_pnts];

   int i;
   point_list.reset();
   CubitVector *vec;
   for( i=0; i<num_pnts; i++ )
   {
      vec = point_list.get_and_step();

      pnt_arr[i].x = vec->x();
      pnt_arr[i].y = vec->y();
      pnt_arr[i].z = vec->z();
   }

   OBBox3 minimal = MinimalBox3 (point_list.size(), pnt_arr);

   //PRINT_INFO( "MinBox center = %lf, %lf, %lf\n", minimal.center.x, minimal.center.y, minimal.center.z );
   //PRINT_INFO( "MinBox axis 1 = %lf, %lf, %lf\n", minimal.axis[0].x, minimal.axis[0].y, minimal.axis[0].z );
   //PRINT_INFO( "MinBox axis 2 = %lf, %lf, %lf\n", minimal.axis[1].x, minimal.axis[1].y, minimal.axis[1].z );
   //PRINT_INFO( "MinBox axis 3 = %lf, %lf, %lf\n", minimal.axis[2].x, minimal.axis[2].y, minimal.axis[2].z );
   //PRINT_INFO( "MinBox extent = %lf, %lf, %lf\n", minimal.extent[0], minimal.extent[1], minimal.extent[2] );

   CubitVector center(minimal.center.x, minimal.center.y, minimal.center.z);
   CubitVector x(minimal.axis[0].x, minimal.axis[0].y, minimal.axis[0].z);
   CubitVector y(minimal.axis[1].x, minimal.axis[1].y, minimal.axis[1].z);
   CubitVector z(minimal.axis[2].x, minimal.axis[2].y, minimal.axis[2].z);
   CubitVector extent(minimal.extent[0], minimal.extent[1], minimal.extent[2]);

   center.next_point( -x, extent.x(), p1 ); p1.next_point( -y, extent.y(), p1 );
   p1.next_point( z, extent.z(), p1 );

   center.next_point( -x, extent.x(), p2 ); p2.next_point( y, extent.y(), p2 );
   p2.next_point( z, extent.z(), p2 );

   center.next_point( x, extent.x(), p3 ); p3.next_point( y, extent.y(), p3 );
   p3.next_point( z, extent.z(), p3 );

   center.next_point( x, extent.x(), p4 ); p4.next_point( -y, extent.y(), p4 );
   p4.next_point( z, extent.z(), p4 );

   center.next_point( -x, extent.x(), p5 ); p5.next_point( -y, extent.y(), p5 );
   p5.next_point( -z, extent.z(), p5 );

   center.next_point( -x, extent.x(), p6 ); p6.next_point( y, extent.y(), p6 );
   p6.next_point( -z, extent.z(), p6 );

   center.next_point( x, extent.x(), p7 ); p7.next_point( y, extent.y(), p7 );
   p7.next_point( -z, extent.z(), p7 );

   center.next_point( x, extent.x(), p8 ); p8.next_point( -y, extent.y(), p8 );
   p8.next_point( -z, extent.z(), p8 );

   delete pnt_arr;

   return CUBIT_SUCCESS;
}

CubitStatus 
AnalyticGeometryTool::min_cyl_box_int( double radius,
                                       CubitVector& axis,
                                       CubitVector& center,
                                       CubitBox& box,
                                       int extension_type,
                                       double extension,
                                       CubitVector &start,
                                       CubitVector &end,
                                       int num_tess_pnts )
                                       
{
   CubitVector box_min = box.minimum();
   CubitVector box_max = box.maximum();
   
   double box_min_pnt[3], box_max_pnt[3], axis_vec[3], center_pnt[3];
   box_min.get_xyz( box_min_pnt ); box_max.get_xyz( box_max_pnt );
   axis.get_xyz( axis_vec ); center.get_xyz( center_pnt );
   
   double start_pnt[3], end_pnt[3];
   
   if( min_cyl_box_int( radius, axis_vec, center_pnt, 
                        box_min_pnt, box_max_pnt, 
                        extension_type, extension, 
                        start_pnt, end_pnt, num_tess_pnts )
                        == CUBIT_FAILURE )
     return CUBIT_FAILURE;
   
   start.set( start_pnt ); end.set( end_pnt );
   
   return CUBIT_SUCCESS;
}

CubitStatus 
AnalyticGeometryTool::min_cyl_box_int( double radius, double axis[3], 
                                      double center[3], 
                                      double box_min[3], double box_max[3], 
                                      int extension_type, double extension, 
                                      double start[3], double end[3], 
                                      int num_tess_pnts )
{
  double cyl_z[3];
  unit_vec( axis, cyl_z );

  //PRINT_INFO( "Axis = %f, %f, %f\n", cyl_z[0], cyl_z[1], cyl_z[2] );
  //PRINT_INFO( "Center = %f, %f, %f\n", center[0], center[1], center[2] );

  // Find transformation matrix to take a point into cylinder's
  // coordinate system
  double cubit2cyl_mtx[4][4], cyl2cubit_mtx[4][4];
  double cyl_x[3], cyl_y[3];
  orth_vecs( cyl_z, cyl_x, cyl_y );
  vecs_to_mtx( cyl_x, cyl_y, cyl_z, center, cyl2cubit_mtx );
  inv_trans_mtx( cyl2cubit_mtx, cubit2cyl_mtx );

  // Setup the circle
  double vec_1[3], vec_2[3];
  orth_vecs( cyl_z, vec_1, vec_2 );
  AGT_3D_Arc arc;
  setup_arc( vec_1, vec_2, center, 0.0, 2.0 * CUBIT_PI, radius, arc );

  // Setup the planes of the box
  double pln_norm[6][3], pln_orig[6][3];
  // Front
  pln_orig[0][0] = 0.0; pln_orig[0][1] = 0.0; pln_orig[0][2] = box_max[2];
  pln_norm[0][0] = 0.0; pln_norm[0][1] = 0.0; pln_norm[0][2] = 1.0;
  // Left
  pln_orig[1][0] = box_min[0]; pln_orig[1][1] = 0.0; pln_orig[1][2] = 0.0;
  pln_norm[1][0] = -1.0; pln_norm[1][1] = 0.0; pln_norm[1][2] = 0.0;
  // Top
  pln_orig[2][0] = 0.0; pln_orig[2][1] = box_max[1]; pln_orig[2][2] = 0.0;
  pln_norm[2][0] = 0.0; pln_norm[2][1] = 1.0; pln_norm[2][2] = 0.0;
  // Right
  pln_orig[3][0] = box_max[0]; pln_orig[3][1] = 0.0; pln_orig[3][2] = 0.0;
  pln_norm[3][0] = 1.0; pln_norm[3][1] = 0.0; pln_norm[3][2] = 0.0;
  // Bottom
  pln_orig[4][0] = 0.0; pln_orig[4][1] = box_min[1]; pln_orig[4][2] = 0.0;
  pln_norm[4][0] = 0.0; pln_norm[4][1] = -1.0; pln_norm[4][2] = 0.0;
  // Back
  pln_orig[5][0] = 0.0; pln_orig[5][1] = 0.0; pln_orig[5][2] = box_min[2];
  pln_norm[5][0] = 0.0; pln_norm[5][1] = 0.0; pln_norm[5][2] = -1.0;

  double z; // Intersection along cylinder's axis
  double min_z = CUBIT_DBL_MAX, max_z = -CUBIT_DBL_MAX;
  
  double t = 0.0, dt;
  dt = 1.0/(double)num_tess_pnts;
  double pnt[3];
  double int_pnt[3];
  double box_tol = 1e-14;
  double box_min_0 = box_min[0]-box_tol;
  double box_min_1 = box_min[1]-box_tol;
  double box_min_2 = box_min[2]-box_tol;
  double box_max_0 = box_max[0]+box_tol;
  double box_max_1 = box_max[1]+box_tol;
  double box_max_2 = box_max[2]+box_tol;

  int i,j;
  for( i=0; i<num_tess_pnts; i++ )
  {
    get_arc_xyz( arc, t, pnt );

    for( j=0; j<6; j++ )
    {
      // Evaluate the intersection at this point
      if( int_ln_pln( pnt, cyl_z, pln_orig[j], pln_norm[j], int_pnt )
        == CUBIT_FAILURE )
        continue;

      // Throw-out if intersection not on physical box
      if( int_pnt[0] < box_min_0 || int_pnt[1] < box_min_1 ||
          int_pnt[2] < box_min_2 || int_pnt[0] > box_max_0 ||
          int_pnt[1] > box_max_1 || int_pnt[2] > box_max_2 ) 
        continue;

      // Find min/max cylinder z on box so far
      // z-distance (in cylinder coordinate system)
      z = cubit2cyl_mtx[0][2]*int_pnt[0] + cubit2cyl_mtx[1][2]*int_pnt[1] + 
        cubit2cyl_mtx[2][2]*int_pnt[2] + cubit2cyl_mtx[3][2];

      if( z < min_z ) min_z = z;
      if( z > max_z ) max_z = z;
      
    }

    t += dt;

  }

  // Check the 8 corners of the box - they are likely min/max's.
  double box_corners[8][3];
  get_box_corners( box_min, box_max, box_corners );
  for( i=0; i<8; i++ )
  {
    // Get the corner in the cylinder csys
    transform_pnt( cubit2cyl_mtx, box_corners[i], pnt );
    // If the pnt is within the circle's radius, check it's z-coord
    // (distance from center)
    if(  sqrt( pnt[0]*pnt[0] + pnt[1]*pnt[1] ) <= radius+box_tol )
    {
      if( pnt[2] < min_z ) min_z = pnt[2];
      if( pnt[2] > max_z ) max_z = pnt[2];
    }
  }

  if( min_z == CUBIT_DBL_MAX || max_z == -CUBIT_DBL_MAX )
  {
    PRINT_ERROR( "Unable to find cylinder/box intersection\n" );
    return CUBIT_FAILURE;
  }

  if( min_z == max_z )
  {
    PRINT_ERROR( "Unable to find cylinder/box intersection\n" );
    return CUBIT_FAILURE;
  }

  //PRINT_INFO( "Min dist = %f\n", min_z );
  //PRINT_INFO( "Max dist = %f\n", max_z );

  // Find the start and end of the cylinder
  next_pnt( center, cyl_z, min_z, start );
  next_pnt( center, cyl_z, max_z, end );

  //PRINT_INFO( "Start = %f, %f, %f\n", start[0], start[1], start[2] );
  //PRINT_INFO( "End = %f, %f, %f\n", end[0], end[1], end[2] );

  // Extend start and end, if necessary
  if( extension_type > 0 )
  {
    double ext_distance = 0.0;
    if( extension_type == 1 ) // percentage
    {
      double cyl_length = dist_pnt_pnt( start, end );
      ext_distance = extension/100.0 * cyl_length;
    }
    else
      ext_distance = extension;

    next_pnt( end, cyl_z, ext_distance, end );
    reverse_vec( cyl_z, cyl_z );
    next_pnt( start, cyl_z, ext_distance, start );
  }

  return CUBIT_SUCCESS;
}

// MAGIC SOFTWARE - see .hpp file
// FILE: minbox2.cpp
//---------------------------------------------------------------------------
double AnalyticGeometryTool::Area (int N, Point2* pt, double angle)
{
    double cs = cos(angle), sn = sin(angle);

    double umin = +cs*pt[0].x+sn*pt[0].y, umax = umin;
    double vmin = -sn*pt[0].x+cs*pt[0].y, vmax = vmin;
    for (int i = 1; i < N; i++)
    {
        double u = +cs*pt[i].x+sn*pt[i].y;
        if ( u < umin )
            umin = u;
        else if ( u > umax )
            umax = u;

        double v = -sn*pt[i].x+cs*pt[i].y;
        if ( v < vmin )
            vmin = v;
        else if ( v > vmax )
            vmax = v;
    }

    double area = (umax-umin)*(vmax-vmin);
    return area;
}
//---------------------------------------------------------------------------
void AnalyticGeometryTool::MinimalBoxForAngle (int N, Point2* pt, double angle, 
                                               OBBox2& box)
{
    double cs = cos(angle), sn = sin(angle);


    double umin = +cs*pt[0].x+sn*pt[0].y, umax = umin;
    double vmin = -sn*pt[0].x+cs*pt[0].y, vmax = vmin;
    for (int i = 1; i < N; i++)
    {
        double u = +cs*pt[i].x+sn*pt[i].y;
        if ( u < umin )
            umin = u;
        else if ( u > umax )
            umax = u;

        double v = -sn*pt[i].x+cs*pt[i].y;
        if ( v < vmin )
            vmin = v;
        else if ( v > vmax )
            vmax = v;
    }

    double umid = 0.5*(umax+umin);
    double vmid = 0.5*(vmax+vmin);
    box.center.x = umid*cs-vmid*sn;
    box.center.y = umid*sn+vmid*cs;
    box.axis[0].x = cs;
    box.axis[0].y = sn;
    box.axis[1].x = -sn;
    box.axis[1].y = cs;
    box.extent[0] = 0.5*(umax-umin);
    box.extent[1] = 0.5*(vmax-vmin);
}
//---------------------------------------------------------------------------
OBBox2 AnalyticGeometryTool::MinimalBox2 (int N, Point2* pt)
{
    OBBox2 box;

    // bracket a minimum for angles in [-pi,pi]
    double angle, area;
    int imin = 0;
    double areaMin = Area(N,pt,angleMin);

    int i;
    for (i = 1; i <= maxPartition; i++)
    {
        angle = angleMin+i*angleRange/maxPartition;
        area = Area(N,pt,angle);
        if ( area < areaMin )
        {
            imin = i;
            areaMin = area;
        }
    }

    double angle0 = angleMin+(imin-1)*angleRange/maxPartition;
    double area0 = Area(N,pt,angle0);
    double angle1 = angleMin+(imin+1)*angleRange/maxPartition;
    double area1 = Area(N,pt,angle1);
    
    // use inverse parabolic interpolation to find the minimum
    for (i = 0; i <= invInterp; i++)
    {
        double angleMid, areaMid;

        // test for convergence (do not change these parameters)
        const double epsilon = 1e-08, tol = 1e-04;
//        const double omtol = 1.0-tol;
        if ( fabs(angle1-angle0) <= 2*tol*fabs(angle)+epsilon )
            break;

        // compute vertex of interpolating parabola
        double dangle0 = angle0-angle, dangle1 = angle1-angle;
        double darea0 = area0-area, darea1 = area1-area;
        double temp0 = dangle0*darea1, temp1 = dangle1*darea0;
        double delta = temp1-temp0;
        if ( fabs(delta) < epsilon )
           break;

        angleMid = angle+0.5*(dangle1*temp1-dangle0*temp0)/(temp1-temp0);

        // update bracket
        if ( angleMid < angle )
        {
            areaMid = Area(N,pt,angleMid);
            if ( areaMid <= area )
            {
                angle1 = angle;
                area1 = area;
                angle = angleMid;
                area = areaMid;
            }
            else
            {
                angle0 = angleMid;
                area0 = areaMid;
            }
        }
        else if ( angleMid > angle )
        {
            areaMid = Area(N,pt,angleMid);
            if ( areaMid <= area )
            {
                angle0 = angle;
                area0 = area;
                angle = angleMid;
                area = areaMid;
            }
            else
            {
                angle1 = angleMid;
                area1 = areaMid;
            }
        }
        else
        {
            // bracket middle already vertex of parabola
            break;
        }
    }

    MinimalBoxForAngle(N,pt,angle,box);
    return box;
}
//---------------------------------------------------------------------------

#ifdef MINBOX2_TEST

#include <stdlib.h>

void main ()
{
    const int N = 128;
    Point2 pt[N];

    for (int i = 0; i < N; i++)
    {
        pt[i].x = rand()/double(RAND_MAX);
        pt[i].y = rand()/double(RAND_MAX);
    }

    OBBox2 minimal = MinimalBox2(N,pt);
}

#endif

#if 1
// FILE: minbox3.cpp
//---------------------------------------------------------------------------
void AnalyticGeometryTool::MatrixToAngleAxis (double** R, double& angle, double axis[3])
{
    // Let (x,y,z) be the unit-length axis and let A be an angle of rotation.
    // The rotation matrix is R = I + sin(A)*P + (1-cos(A))*P^2 where
    // I is the identity and
    //
    //       +-        -+
    //   P = |  0 +z -y |
    //       | -z  0 +x |
    //       | +y -x  0 |
    //       +-        -+
    //
    // Some algebra will show that
    //
    //   cos(A) = (trace(R)-1)/2  and  R - R^t = 2*sin(A)*P
    //
    // In the event that A = pi, R-R^t = 0 which prevents us from extracting
    // the axis through P.  Instead note that R = I+2*P^2 when A = pi, so
    // P^2 = (R-I)/2.  The diagonal entries of P^2 are x^2-1, y^2-1, and z^2-1.
    // We can solve these for axis (x,y,z).  Because the angle is pi, it does
    // not matter which sign you choose on the square roots.

    double trace = R[0][0]+R[1][1]+R[2][2];
    double cs = 0.5*(trace-1.0);
    if ( -1 < cs )
    {
        if ( cs < 1 )
            angle = acos(cs);
        else
            angle = 0;
    }
    else
    {
        angle = AGT_PI;
    }

    axis[0] = R[1][2]-R[2][1];
    axis[1] = R[2][0]-R[0][2];
    axis[2] = R[0][1]-R[1][0];
    double length = sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
    const double epsilon = 1e-06;
    if ( length > epsilon )
    {
        axis[0] /= length;
        axis[1] /= length;
        axis[2] /= length;
    }
    else  // angle is 0 or pi
    {
        if ( angle > 1.0 )  // any number strictly between 0 and pi works
        {
            // angle must be pi
            axis[0] = sqrt(0.5*(1.0+R[0][0]));
            axis[1] = sqrt(0.5*(1.0+R[1][1]));
            axis[2] = sqrt(0.5*(1.0+R[2][2]));

            // determine signs of axis components
            double test[3];
            test[0] = R[0][0]*axis[0]+R[0][1]*axis[1]+R[0][2]*axis[2]-axis[0];
            test[1] = R[1][0]*axis[0]+R[1][1]*axis[1]+R[1][2]*axis[2]-axis[1];
            test[2] = R[2][0]*axis[0]+R[2][1]*axis[1]+R[2][2]*axis[2]-axis[2];
            length = test[0]*test[0]+test[1]*test[1]+test[2]*test[2];
            if ( length < epsilon )
                return;

            axis[1] = -axis[1];
            test[0] = R[0][0]*axis[0]+R[0][1]*axis[1]+R[0][2]*axis[2]-axis[0];
            test[1] = R[1][0]*axis[0]+R[1][1]*axis[1]+R[1][2]*axis[2]-axis[1];
            test[2] = R[2][0]*axis[0]+R[2][1]*axis[1]+R[2][2]*axis[2]-axis[2];
            length = test[0]*test[0]+test[1]*test[1]+test[2]*test[2];
            if ( length < epsilon )
                return;

            axis[2] = -axis[2];
            test[0] = R[0][0]*axis[0]+R[0][1]*axis[1]+R[0][2]*axis[2]-axis[0];
            test[1] = R[1][0]*axis[0]+R[1][1]*axis[1]+R[1][2]*axis[2]-axis[1];
            test[2] = R[2][0]*axis[0]+R[2][1]*axis[1]+R[2][2]*axis[2]-axis[2];
            length = test[0]*test[0]+test[1]*test[1]+test[2]*test[2];
            if ( length < epsilon )
                return;

            axis[1] = -axis[1];
            test[0] = R[0][0]*axis[0]+R[0][1]*axis[1]+R[0][2]*axis[2]-axis[0];
            test[1] = R[1][0]*axis[0]+R[1][1]*axis[1]+R[1][2]*axis[2]-axis[1];
            test[2] = R[2][0]*axis[0]+R[2][1]*axis[1]+R[2][2]*axis[2]-axis[2];
            length = test[0]*test[0]+test[1]*test[1]+test[2]*test[2];
            if ( length < epsilon )
                return;
        }
        else
        {
            // Angle is zero, matrix is the identity, no unique axis, so
            // return (1,0,0) for as good a guess as any.
            axis[0] = 1.0;
            axis[1] = 0.0;
            axis[2] = 0.0;
        }
    }
}
//---------------------------------------------------------------------------
void AnalyticGeometryTool::AngleAxisToMatrix (double angle, double axis[3], double R[3][3])
{
    double cs = cos(angle), sn = sin(angle);
    double length = sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
    double x = axis[0]/length;
    double y = axis[1]/length;
    double z = axis[2]/length;
    double omc = 1.0-cs;
    double x2 = x*x, y2 = y*y, z2 = z*z;
    double xy = x*y, xz = x*z, yz = y*z;
    double snx = sn*x, sny = sn*y, snz = sn*z;
    
    R[0][0] = 1.0-omc*(y2+z2);
    R[0][1] = +snz+omc*xy;
    R[0][2] = -sny+omc*xz;
    R[1][0] = -snz+omc*xy;
    R[1][1] = 1.0-omc*(x2+z2);
    R[1][2] = +snx+omc*yz;
    R[2][0] = +sny+omc*xz;
    R[2][1] = -snx+omc*yz;
    R[2][2] = 1.0-omc*(x2+y2);
}
//---------------------------------------------------------------------------
double AnalyticGeometryTool::Volume (int N, Point3* pt, double angle[3])
{
    double cs0 = cos(angle[0]), sn0 = sin(angle[0]);
    double cs1 = cos(angle[1]), sn1 = sin(angle[1]);
    double axis[3] = { cs0*sn1, sn0*sn1, cs1 };
    double rot[3][3];
    AngleAxisToMatrix(angle[2],axis,rot);

    double min[3] =
    {
        rot[0][0]*pt[0].x+rot[1][0]*pt[0].y+rot[2][0]*pt[0].z,
        rot[0][1]*pt[0].x+rot[1][1]*pt[0].y+rot[2][1]*pt[0].z,
        rot[0][2]*pt[0].x+rot[1][2]*pt[0].y+rot[2][2]*pt[0].z
    };

    double max[3] = { min[0], min[1], min[2] };

    for (int i = 1; i < N; i++)
    {
        double test[3] =
        {
            rot[0][0]*pt[i].x+rot[1][0]*pt[i].y+rot[2][0]*pt[i].z,
            rot[0][1]*pt[i].x+rot[1][1]*pt[i].y+rot[2][1]*pt[i].z,
            rot[0][2]*pt[i].x+rot[1][2]*pt[i].y+rot[2][2]*pt[i].z
        };

        if ( test[0] < min[0] )
            min[0] = test[0];
        else if ( test[0] > max[0] )
            max[0] = test[0];

        if ( test[1] < min[1] )
            min[1] = test[1];
        else if ( test[1] > max[1] )
            max[1] = test[1];

        if ( test[2] < min[2] )
            min[2] = test[2];
        else if ( test[2] > max[2] )
            max[2] = test[2];
    }

    double volume = (max[0]-min[0])*(max[1]-min[1])*(max[2]-min[2]);
    return volume;
}
//---------------------------------------------------------------------------
void AnalyticGeometryTool::MinimalBoxForAngles (int N, Point3* pt, double angle[3],
                                 OBBox3& box)
{
    double cs0 = cos(angle[0]), sn0 = sin(angle[0]);
    double cs1 = cos(angle[1]), sn1 = sin(angle[1]);
    double axis[3] = { cs0*sn1, sn0*sn1, cs1 };
    double rot[3][3];
    AngleAxisToMatrix(angle[2],axis,rot);

    double min[3] =
    {
        rot[0][0]*pt[0].x+rot[1][0]*pt[0].y+rot[2][0]*pt[0].z,
        rot[0][1]*pt[0].x+rot[1][1]*pt[0].y+rot[2][1]*pt[0].z,
        rot[0][2]*pt[0].x+rot[1][2]*pt[0].y+rot[2][2]*pt[0].z
    };

    double max[3] = { min[0], min[1], min[2] };

    for (int i = 1; i < N; i++)
    {
        double test[3] =
        {
            rot[0][0]*pt[i].x+rot[1][0]*pt[i].y+rot[2][0]*pt[i].z,
            rot[0][1]*pt[i].x+rot[1][1]*pt[i].y+rot[2][1]*pt[i].z,
            rot[0][2]*pt[i].x+rot[1][2]*pt[i].y+rot[2][2]*pt[i].z
        };

        if ( test[0] < min[0] )
            min[0] = test[0];
        else if ( test[0] > max[0] )
            max[0] = test[0];

        if ( test[1] < min[1] )
            min[1] = test[1];
        else if ( test[1] > max[1] )
            max[1] = test[1];

        if ( test[2] < min[2] )
            min[2] = test[2];
        else if ( test[2] > max[2] )
            max[2] = test[2];
    }

    double mid[3] =
    {
        0.5*(max[0]+min[0]), 0.5*(max[1]+min[1]), 0.5*(max[2]+min[2])
    };

    box.center.x = mid[0]*rot[0][0]+mid[1]*rot[0][1]+mid[2]*rot[0][2];
    box.center.y = mid[0]*rot[1][0]+mid[1]*rot[1][1]+mid[2]*rot[1][2];
    box.center.z = mid[0]*rot[2][0]+mid[1]*rot[2][1]+mid[2]*rot[2][2];
    box.axis[0].x = rot[0][0];
    box.axis[0].y = rot[1][0];
    box.axis[0].z = rot[2][0];
    box.axis[1].x = rot[0][1];
    box.axis[1].y = rot[1][1];
    box.axis[1].z = rot[2][1];
    box.axis[2].x = rot[0][2];
    box.axis[2].y = rot[1][2];
    box.axis[2].z = rot[2][2];
    box.extent[0] = 0.5*(max[0]-min[0]);
    box.extent[1] = 0.5*(max[1]-min[1]);
    box.extent[2] = 0.5*(max[2]-min[2]);
}
//---------------------------------------------------------------------------
void AnalyticGeometryTool::GetInterval (double A[3], double D[3], double& tmin, double& tmax)
{
    //static const double angle_min[3] = { -AGT_PI, 0.0, 0.0 };
    //static const double angle_max[3] = {  AGT_PI,  AGT_PI,  AGT_PI };
    //The pgCC compiler running on solars and cross-compiling for
    //janus has a bug such that it dies if the initialization
    //of angle_min is mixed symbolic and literal constants, so
    //define a symbolid constant for 0.0.
    // -- J.Kraftcheck, 06/05/2001
    static const double ZERO = 0.0;
    static const double angle_min[3] = { -AGT_PI, ZERO, ZERO };
    static const double angle_max[3] = {  AGT_PI, AGT_PI, AGT_PI };

    tmin = -DBL_MAX;
    tmax = +DBL_MAX;

    for (int i = 0; i < 3; i++)
    {
        const double epsilon = 1e-08;
        double b0 = angle_min[i]-A[i];
        double b1 = angle_max[i]-A[i];

        double inv, tmp;
        if ( D[i] > epsilon )
        {
            inv = 1.0/D[i];
            tmp = inv*b0;
            if ( tmp > tmin )
                tmin = tmp;
            tmp = inv*b1;
            if ( tmp < tmax )
                tmax = tmp;
        }
        else if ( D[i] < -epsilon )
        {
            inv = 1.0/D[i];
            tmp = inv*b0;
            if ( tmp < tmax )
                tmax = tmp;
            tmp = inv*b1;
            if ( tmp > tmin )
                tmin = tmp;
        }
    }

    if( tmin == -DBL_MAX || tmax == DBL_MAX ) // Added by SRS 10-6-2000
    { 
       //PRINT_WARNING( "tmin/tmax not set\n" );
       tmin = -1.0;
       tmax = 1.0;
    }
}
//---------------------------------------------------------------------------
void AnalyticGeometryTool::Combine (double result[3], double A[3], double t, double D[3])
{
    for (int i = 0; i < 3; i++)
        result[i] = A[i]+t*D[i];
}
//---------------------------------------------------------------------------
double AnalyticGeometryTool::MinimizeOnInterval (int N, Point3* pt, double A[3], double D[3])
{
    // compute intersection of line A+t*D with domain of function
    double tmin, tmax;
    GetInterval(A,D,tmin,tmax);
    double tran = tmax-tmin;
    double angle[3];

    if( tran == 0.0 ) // Added by SRS 10-6-2000
       tran = 1.0;

    // bracket a minimum for angles in [A+tmin*D,A+tmax*D]
    double t = 0.0;
    double volumeMin = Volume(N,pt,A);
    double volume;

    const int max_partition = 64;
    int i, imin;
    for (i = 0, imin = -1; i <= max_partition; i++)
    {
        t = tmin+i*tran/max_partition;
        Combine(angle,A,t,D);

        volume = Volume(N,pt,angle);
        if ( volume < volumeMin )
        {
            imin = i;
            volumeMin = volume;
        }
    }

    if ( imin != -1 )
    {
        t = tmin+imin*tran/max_partition;
    }
    else
    {
        t = 0.0;

        // interval in which t=0 lies
        imin = int(-tmin*max_partition/tran+0.5);
    }
    volume = volumeMin;

    double t0 = tmin+(imin-1)*tran/max_partition;
    Combine(angle,A,t0,D);
    double volume0 = Volume(N,pt,angle);

    double t1 = tmin+(imin+1)*tran/max_partition;
    Combine(angle,A,t1,D);
    double volume1 = Volume(N,pt,angle);
    
    // use inverse parabolic interpolation to find the minimum
    const int inv_interp = 64;
    for (i = 0; i <= inv_interp; i++)
    {
        double tMid, volumeMid;

        // test for convergence (do not change these parameters)
        const double epsilon = 1e-08, tol = 1e-04;
//        const double omtol = 1.0-tol;
        if ( fabs(t1-t0) <= 2*tol*fabs(t)+epsilon )
            break;

        // compute vertex of interpolating parabola
        double dt0 = t0-t, dt1 = t1-t;
        double dvolume0 = volume0-volume, dvolume1 = volume1-volume;
        double temp0 = dt0*dvolume1, temp1 = dt1*dvolume0;
        double delta = temp1-temp0;
        if ( fabs(delta) < epsilon )
           break;

        tMid = t+0.5*(dt1*temp1-dt0*temp0)/(temp1-temp0);

        // update bracket
        if ( tMid < t )
        {
            Combine(angle,A,tMid,D);
            volumeMid = Volume(N,pt,angle);
            if ( volumeMid <= volume )
            {
                t1 = t;
                volume1 = volume;
                t = tMid;
                volume = volumeMid;
            }
            else
            {
                t0 = tMid;
                volume0 = volumeMid;
            }
        }
        else if ( tMid > t )
        {
            Combine(angle,A,tMid,D);
            volumeMid = Volume(N,pt,angle);
            if ( volumeMid <= volume )
            {
                t0 = t;
                volume0 = volume;
                t = tMid;
                volume = volumeMid;
            }
            else
            {
                t1 = tMid;
                volume1 = volumeMid;
            }
        }
        else
        {
            // bracket middle already vertex of parabola
            break;
        }
    }

    Combine(A,A,t,D);
    return volume;
}
//---------------------------------------------------------------------------
double AnalyticGeometryTool::MinimizeOnLattice (int N, Point3* pt, double A[3], int layers,
                                 double thickness)
{
    int xmin = 0, ymin = 0, zmin = 0;
    double volume = Volume(N,pt,A);

    double angle[3];
    for (int z = -layers; z <= layers; z++)
    {
        angle[2] = A[2]+thickness*z/layers;
        for (int y = -layers; y <= layers; y++)
        {
            angle[1] = A[1]+thickness*y/layers;
            for (int x = -layers; x <= layers; x++)
            {
                angle[0] = A[0]+thickness*x/layers;

                double v = Volume(N,pt,angle);
                if ( v < volume )
                {
                    xmin = x;
                    ymin = y;
                    zmin = z;
                    volume = v;
                }
            }
        }
    }

    A[0] += thickness*xmin/layers;
    A[1] += thickness*ymin/layers;
    A[2] += thickness*zmin/layers;

    return volume;
}
//---------------------------------------------------------------------------
void AnalyticGeometryTool::InitialGuess (int N, Point3* pt, double angle[3])
{
    int i;

    // compute mean of points
    double xsum = 0.0f, ysum = 0.0f, zsum = 0.0f;;
    for (i = 0; i < N; i++)
    {
        xsum += pt[i].x;
        ysum += pt[i].y;
        zsum += pt[i].z;
    }
    double xmean = xsum/N;
    double ymean = ysum/N;
    double zmean = zsum/N;

    // compute covariances of points
    double xxsum = 0.0f, xysum = 0.0f, xzsum = 0.0f;
    double yysum = 0.0f, yzsum = 0.0f, zzsum = 0.0f;
    for (i = 0; i < N; i++)
    {
        double dx = pt[i].x - xmean;
        double dy = pt[i].y - ymean;
        double dz = pt[i].z - zmean;
        xxsum += dx*dx;
        xysum += dx*dy;
        xzsum += dx*dz;
        yysum += dy*dy;
        yzsum += dy*dz;
        zzsum += dz*dz;
    }
    double xxcov = xxsum/N;
    double xycov = xysum/N;
    double xzcov = xzsum/N;
    double yycov = yysum/N;
    double yzcov = yzsum/N;
    double zzcov = zzsum/N;

    // compute eigenvectors for covariance matrix
    mgcEigenD eig(3);
    eig.Matrix(0,0) = xxcov;
    eig.Matrix(0,1) = xycov;
    eig.Matrix(0,2) = xzcov;
    eig.Matrix(1,0) = xycov;
    eig.Matrix(1,1) = yycov;
    eig.Matrix(1,2) = yzcov;
    eig.Matrix(2,0) = xzcov;
    eig.Matrix(2,1) = yzcov;
    eig.Matrix(2,2) = zzcov;
    eig.EigenStuff3();

    // Use eigenvectors as the box axes.  Eigenmatrix must not have a
    // reflection component, thus the check for negative determinant.
    const double epsilon = 1e-06;
    double** R = (double**)eig.Eigenvector();
    double det =
        +R[0][0]*R[1][1]*R[2][2]
        +R[0][1]*R[1][2]*R[2][0]
        +R[0][2]*R[1][0]*R[2][1]
        -R[0][2]*R[1][1]*R[2][0]
        -R[0][1]*R[1][0]*R[2][2]
        -R[0][0]*R[1][2]*R[2][1];
    if ( det < 0.0 )
    {
        R[0][0] = -R[0][0];
        R[1][0] = -R[1][0];
        R[2][0] = -R[2][0];
    }

    // extract angles from rotation axis = (cos(u)sin(v),sin(u)sin(v),cos(v))
    double axis[3];
    MatrixToAngleAxis(R,angle[2],axis);
    if ( -1+epsilon < axis[2] )
    {
        if ( axis[2] < 1-epsilon )
        {
            angle[0] = atan2(axis[1],axis[0]);
            angle[1] = acos(axis[2]);
        }
        else
        {
            angle[0] = 0;
            angle[1] = 0;
        }
    }
    else
    {
        angle[0] = 0;
        angle[1] = AGT_PI;
    }
}
//---------------------------------------------------------------------------
OBBox3 AnalyticGeometryTool::MinimalBox3 (int N, Point3* pt)
{
    // compute a good initial guess for an oriented bounding box
    double angle[3];
    InitialGuess(N,pt,angle);
    double oldVolume = Volume(N,pt,angle);
    double saveAngle[3] = { angle[0], angle[1], angle[2] };

    // Powell's direction set method
    double U[3][3], volume;
    const int maxiters = 3*32;
    for (int iter = 0; iter < maxiters; iter++)
    {
        // reset directions to avoid linear dependence degeneration
        if ( iter % 3 == 0 )
        {
            U[0][0] = 1.0;  U[0][1] = 0.0;  U[0][2] = 0.0;
            U[1][0] = 0.0;  U[1][1] = 1.0;  U[1][2] = 0.0;
            U[2][0] = 0.0;  U[2][1] = 0.0;  U[2][2] = 1.0;
        }

        // find minima in specified directions
        for (int d = 0; d < 3; d++)
            volume = MinimizeOnInterval(N,pt,angle,U[d]);

        // estimate a conjugate direction
        double conj[3] =
        {
            angle[0]-saveAngle[0],
            angle[1]-saveAngle[1],
            angle[2]-saveAngle[2]
        };
        double length = sqrt(conj[0]*conj[0]+conj[1]*conj[1]+conj[2]*conj[2]);
        if ( length >= 1e-06 )
        {
            double invLen = 1.0/length;
            conj[0] *= invLen;
            conj[1] *= invLen;
            conj[2] *= invLen;
            
            // minimize in conjugate direction
            volume = MinimizeOnInterval(N,pt,angle,conj);
        }
        else
        {
            // Possible local, but not global, minimum.  Search nearby for
            // a smaller volume.
            volume = MinimizeOnLattice(N,pt,angle,2,0.0001);
            volume = MinimizeOnLattice(N,pt,angle,2,0.0010);
            volume = MinimizeOnLattice(N,pt,angle,2,0.0100);
            volume = MinimizeOnLattice(N,pt,angle,2,0.1000);
        }

        // test for convergence
        const double epsilon = 1e-04;
        double diff = fabs(volume-oldVolume);
        if ( diff <= epsilon )
        {
            // Possible local, but not global, minimum.  Search nearby for
            // a smaller volume.
            volume = MinimizeOnLattice(N,pt,angle,2,0.0001);
            volume = MinimizeOnLattice(N,pt,angle,2,0.0010);
            volume = MinimizeOnLattice(N,pt,angle,2,0.0100);
            volume = MinimizeOnLattice(N,pt,angle,2,0.1000);
            diff = fabs(volume-oldVolume);
            if ( diff <= epsilon )
                break;
        }

        // cycle the directions and add conjugate direction to list
        U[0][0] = U[1][0];  U[0][1] = U[1][1];  U[0][2] = U[1][2];
        U[1][0] = U[2][0];  U[1][1] = U[2][1];  U[1][2] = U[2][2];
        U[2][0] = conj[0];  U[2][1] = conj[1];  U[2][2] = conj[2];

        // set parameters for next pass
        oldVolume = volume;
        saveAngle[0] = angle[0];
        saveAngle[1] = angle[1];
        saveAngle[2] = angle[2];
    }

    OBBox3 box;
    MinimalBoxForAngles(N,pt,angle,box);
    return box;
}
//---------------------------------------------------------------------------

#ifdef MINBOX3_TEST

#define RAND (rand()/double(RAND_MAX))

void main ()
{
    // build box with axes parallel to coordinate axes
    const int N = 16;
    const double ex = 1.0;
    const double ey = 2.0;
    const double ez = 3.0;
    Point3 pt[N];
    pt[0].x = -ex;  pt[0].y = -ey;  pt[0].z = -ez;
    pt[1].x = -ex;  pt[1].y = +ey;  pt[1].z = -ez;
    pt[2].x = +ex;  pt[2].y = +ey;  pt[2].z = -ez;
    pt[3].x = +ex;  pt[3].y = -ey;  pt[3].z = -ez;
    pt[4].x = -ex;  pt[4].y = -ey;  pt[4].z = +ez;
    pt[5].x = -ex;  pt[5].y = +ey;  pt[5].z = +ez;
    pt[6].x = +ex;  pt[6].y = +ey;  pt[6].z = +ez;
    pt[7].x = +ex;  pt[7].y = -ey;  pt[7].z = +ez;

    for (int k = 8; k < N; k++)
    {
        // generate random points inside box to confound initial Gaussian fit
        pt[k].x = -ex+2.0*ex*RAND;
        pt[k].y = -ey+2.0*ey*RAND;
        pt[k].z = -ez+2.0*ez*RAND;
    }

    double maxNorm = 0.0;
    int iMaxNorm = -1;
    for (int iter = 0; iter < 1024; iter++)
    {
        // build arbitrary rotation matrix
        double angle = RAND;
        double line[3] = { RAND, RAND, RAND };
        double rot[3][3];
        double length = sqrt(line[0]*line[0]+line[1]*line[1]+line[2]*line[2]);
        line[0] /= length;
        line[1] /= length;
        line[2] /= length;
        AngleAxisToMatrix(angle,line,rot);
        
        // rotate box
        Point3 rpt[N];
        for (int i = 0; i < N; i++)
        {
            rpt[i].x = rot[0][0]*pt[i].x+rot[0][1]*pt[i].y+rot[0][2]*pt[i].z;
            rpt[i].y = rot[1][0]*pt[i].x+rot[1][1]*pt[i].y+rot[1][2]*pt[i].z;
            rpt[i].z = rot[2][0]*pt[i].x+rot[2][1]*pt[i].y+rot[2][2]*pt[i].z;
        }
        
        OBBox3 minimal = MinimalBox3(N,rpt);

        int index[3];
        if ( minimal.extent[0] <= minimal.extent[1] )
        {
            if ( minimal.extent[1] <= minimal.extent[2] )
            {
                index[0] = 0;
                index[1] = 1;
                index[2] = 2;
            }
            else
            {
                if ( minimal.extent[0] <= minimal.extent[2] )
                {
                    index[0] = 0;
                    index[1] = 2;
                    index[2] = 1;
                }
                else
                {
                    index[0] = 2;
                    index[1] = 0;
                    index[2] = 1;
                }
            }
        }
        else
        {
            if ( minimal.extent[0] <= minimal.extent[2] )
            {
                index[0] = 1;
                index[1] = 0;
                index[2] = 2;
            }
            else
            {
                if ( minimal.extent[1] <= minimal.extent[2] )
                {
                    index[0] = 1;
                    index[1] = 2;
                    index[2] = 0;
                }
                else
                {
                    index[0] = 2;
                    index[1] = 1;
                    index[2] = 0;
                }
            }
        }

        double dx = ex-minimal.extent[index[0]];
        double dy = ey-minimal.extent[index[1]];
        double dz = ez-minimal.extent[index[2]];
        double norm = sqrt(dx*dx+dy*dy+dz*dz);
        if ( norm > maxNorm )
        {
            maxNorm = norm;
            iMaxNorm = iter;
        }
    }
}

#endif // MINBOX3_TEST

// FILE: eigen.cpp
//===========================================================================
// error handling
int mgcEigen::verbose1 = 0;
unsigned mgcEigen::error = 0;
const unsigned mgcEigen::invalid_size      = 0x00000001;
const unsigned mgcEigen::allocation_failed = 0x00000002;
const unsigned mgcEigen::ql_exceeded       = 0x00000004;
const char* mgcEigen::message[3] = {
        "invalid matrix size",
        "allocation failed",
        "QL algorithm - exceeded maximum iterations"
};
//---------------------------------------------------------------------------
void mgcEigen::
Tridiagonal2 (float** pmat, float* pdiag, float* psubd)
{
        // matrix is already tridiagonal

        pdiag[0] = pmat[0][0];
        pdiag[1] = pmat[1][1];
        psubd[0] = pmat[0][1];
        psubd[1] = 0;
        pmat[0][0] = 1;  pmat[0][1] = 0;
        pmat[1][0] = 0;  pmat[1][1] = 1;
}
//---------------------------------------------------------------------------
void mgcEigen::
Tridiagonal3 (float** pmat, float* pdiag, float* psubd)
{
        float a = pmat[0][0], b = pmat[0][1], c = pmat[0][2],
                                                 d = pmat[1][1], e = pmat[1][2],
                                                                                f = pmat[2][2];

        pdiag[0] = a;
        psubd[2] = 0;
        if ( c != 0 ) {
                float ell = float(sqrt(b*b+c*c));
                b /= ell;
                c /= ell;
                float q = 2*b*e+c*(f-d);
                pdiag[1] = d+c*q;
                pdiag[2] = f-c*q;
                psubd[0] = ell;
                psubd[1] = e-b*q;
                pmat[0][0] = 1; pmat[0][1] = 0; pmat[0][2] = 0;
                pmat[1][0] = 0; pmat[1][1] = b; pmat[1][2] = c;
                pmat[2][0] = 0; pmat[2][1] = c; pmat[2][2] = -b;
        }
        else {
                pdiag[1] = d;
                pdiag[2] = f;
                psubd[0] = b;
                psubd[1] = e;
                pmat[0][0] = 1; pmat[0][1] = 0; pmat[0][2] = 0;
                pmat[1][0] = 0; pmat[1][1] = 1; pmat[1][2] = 0;
                pmat[2][0] = 0; pmat[2][1] = 0; pmat[2][2] = 1;
        }
}
//---------------------------------------------------------------------------
void mgcEigen::
Tridiagonal4 (float** pmat, float* pdiag, float* psubd)
{
        // save pmatrix M
        float
        a = pmat[0][0], b = pmat[0][1], c = pmat[0][2], d = pmat[0][3],
                                   e = pmat[1][1], f = pmat[1][2], g = pmat[1][3],
                                                                  h = pmat[2][2], i = pmat[2][3],
                                                                                                 j = pmat[3][3];

        pdiag[0] = a;
        psubd[3] = 0;

        pmat[0][0] = 1; pmat[0][1] = 0; pmat[0][2] = 0; pmat[0][3] = 0;
        pmat[1][0] = 0;
        pmat[2][0] = 0;
        pmat[3][0] = 0;

        if ( c != 0 || d != 0 ) {
                float q11, q12, q13;
                float q21, q22, q23;
                float q31, q32, q33;

                // build column Q1
                float len = float(sqrt(b*b+c*c+d*d));
                q11 = b/len;
                q21 = c/len;
                q31 = d/len;

                psubd[0] = len;

                // compute S*Q1
                float v0 = e*q11+f*q21+g*q31;
                float v1 = f*q11+h*q21+i*q31;
                float v2 = g*q11+i*q21+j*q31;

                pdiag[1] = q11*v0+q21*v1+q31*v2;

        // build column Q3 = Q1x(S*Q1)
                q13 = q21*v2-q31*v1;
                q23 = q31*v0-q11*v2;
                q33 = q11*v1-q21*v0;
                len = float(sqrt(q13*q13+q23*q23+q33*q33));
                if ( len > 0 ) {
                        q13 /= len;
                        q23 /= len;
                        q33 /= len;

                        // build column Q2 = Q3xQ1
                        q12 = q23*q31-q33*q21; 
                        q22 = q33*q11-q13*q31;
                        q32 = q13*q21-q23*q11;

                        v0 = q12*e+q22*f+q32*g;
                        v1 = q12*f+q22*h+q32*i;
                        v2 = q12*g+q22*i+q32*j;
                        psubd[1] = q11*v0+q21*v1+q31*v2;
                        pdiag[2] = q12*v0+q22*v1+q32*v2;
                        psubd[2] = q13*v0+q23*v1+q33*v2;

                        v0 = q13*e+q23*f+q33*g;
                        v1 = q13*f+q23*h+q33*i;
                        v2 = q13*g+q23*i+q33*j;
                        pdiag[3] = q13*v0+q23*v1+q33*v2;
                }
                else {  // S*Q1 parallel to Q1, choose any valid Q2 and Q3
                        psubd[1] = 0;

                        len = q21*q21+q31*q31;
                        if ( len > 0 ) {
                                float tmp = q11-1;
                                q12 = -q21;
                                q22 = 1+tmp*q21*q21/len;
                q32 = tmp*q21*q31/len;

                                q13 = -q31;
                                q23 = q32;
                                q33 = 1+tmp*q31*q31/len;

                                v0 = q12*e+q22*f+q32*g;
                                v1 = q12*f+q22*h+q32*i;
                                v2 = q12*g+q22*i+q32*j;
                                pdiag[2] = q12*v0+q22*v1+q32*v2;
                                psubd[2] = q13*v0+q23*v1+q33*v2;

                                v0 = q13*e+q23*f+q33*g;
                                v1 = q13*f+q23*h+q33*i;
                                v2 = q13*g+q23*i+q33*j;
                                pdiag[3] = q13*v0+q23*v1+q33*v2;
                        }
                        else {  // Q1 = (+-1,0,0)
                                q12 = 0; q22 = 1; q32 = 0;
                                q13 = 0; q23 = 0; q33 = 1;

                                pdiag[2] = h;
                                pdiag[3] = j;
                                psubd[2] = i;
                        }
                }

                pmat[1][1] = q11; pmat[1][2] = q12; pmat[1][3] = q13;
                pmat[2][1] = q21; pmat[2][2] = q22; pmat[2][3] = q23;
                pmat[3][1] = q31; pmat[3][2] = q32; pmat[3][3] = q33;
        }
        else {
                pdiag[1] = e;
                psubd[0] = b;
                pmat[1][1] = 1;
                pmat[2][1] = 0;
                pmat[3][1] = 0; 

                if ( g != 0 ) {
                        float ell = float(sqrt(f*f+g*g));
                        f /= ell;
                        g /= ell;
                        float Q = 2*f*i+g*(j-h);

                        pdiag[2] = h+g*Q;
                        pdiag[3] = j-g*Q;
                        psubd[1] = ell;
                        psubd[2] = i-f*Q;
                        pmat[1][2] = 0;  pmat[1][3] = 0;
                        pmat[2][2] = f;  pmat[2][3] = g;
                        pmat[3][2] = g;  pmat[3][3] = -f;
                }
                else {
                        pdiag[2] = h;
                        pdiag[3] = j;
                        psubd[1] = f;
                        psubd[2] = i;
                        pmat[1][2] = 0;  pmat[1][3] = 0;
                        pmat[2][2] = 1;  pmat[2][3] = 0;
                        pmat[3][2] = 0;  pmat[3][3] = 1;
                }
        }
}
//---------------------------------------------------------------------------
void mgcEigen::
TridiagonalN (int n, float** pmat, float* pdiag, float* psubd)
{
        int i, j, k, ell;

        for (i = n-1, ell = n-2; i >= 1; i--, ell--) {
                float h = 0, scale = 0;

                if ( ell > 0 ) {
                        for (k = 0; k <= ell; k++)
                                scale += float(fabs(pmat[i][k]));
                        if ( scale == 0 )
                                psubd[i] = pmat[i][ell];
                        else {
                                for (k = 0; k <= ell; k++) {
                                        pmat[i][k] /= scale;
                                        h += pmat[i][k]*pmat[i][k];
                                }
                                float f = pmat[i][ell];
                                float g = ( f > 0 ? -float(sqrt(h)) : float(sqrt(h)) );
                                psubd[i] = scale*g;
                                h -= f*g;
                                pmat[i][ell] = f-g;
                                f = 0;
                                for (j = 0; j <= ell; j++) {
                                        pmat[j][i] = pmat[i][j]/h;
                                        g = 0;
                                        for (k = 0; k <= j; k++)
                                                g += pmat[j][k]*pmat[i][k];
                                        for (k = j+1; k <= ell; k++)
                                                g += pmat[k][j]*pmat[i][k];
                                        psubd[j] = g/h;
                                        f += psubd[j]*pmat[i][j];
                                }
                                float hh = f/(h+h);
                                for (j = 0; j <= ell; j++) {
                                        f = pmat[i][j];
                                        psubd[j] = g = psubd[j] - hh*f;
                                        for (k = 0; k <= j; k++)
                                                pmat[j][k] -= f*psubd[k]+g*pmat[i][k];
                                }
            }
                }
                else
                        psubd[i] = pmat[i][ell];

                pdiag[i] = h;
        }

        pdiag[0] = psubd[0] = 0;
        for (i = 0, ell = -1; i <= n-1; i++, ell++) {
                if ( pdiag[i] ) {
                        for (j = 0; j <= ell; j++) {
                                float sum = 0;
                                for (k = 0; k <= ell; k++)
                                        sum += pmat[i][k]*pmat[k][j];
                                for (k = 0; k <= ell; k++)
                                        pmat[k][j] -= sum*pmat[k][i];
                        }
                }
                pdiag[i] = pmat[i][i];
                pmat[i][i] = 1;
                for (j = 0; j <= ell; j++)
                        pmat[j][i] = pmat[i][j] = 0;
        }

        // re-ordering if mgcEigen::QLAlgorithm is used subsequently
        for (i = 1, ell = 0; i < n; i++, ell++)
                psubd[ell] = psubd[i];
        psubd[n-1] = 0;
}
//---------------------------------------------------------------------------
void mgcEigen::
QLAlgorithm (int n, float* pdiag, float* psubd, float** pmat)
{
        const int eigen_maxiter = 30;

        for (int ell = 0; ell < n; ell++) {
                int iter;
                for (iter = 0; iter < eigen_maxiter; iter++) {
                        int m;
                        for (m = ell; m <= n-2; m++) {
                                float dd = float(fabs(pdiag[m])+fabs(pdiag[m+1]));
                                if ( (float)(fabs(psubd[m])+dd) == dd )
                                        break;
                        }
                        if ( m == ell )
                                break;

                        float g = (pdiag[ell+1]-pdiag[ell])/(2*psubd[ell]);
                        float r = float(sqrt(g*g+1));
                        if ( g < 0 )
                                g = pdiag[m]-pdiag[ell]+psubd[ell]/(g-r);
                        else
                                g = pdiag[m]-pdiag[ell]+psubd[ell]/(g+r);
                        float s = 1, c = 1, p = 0;
                        for (int i = m-1; i >= ell; i--) {
                                float f = s*psubd[i], b = c*psubd[i];
                                if ( fabs(f) >= fabs(g) ) {
                                        c = g/f;
                                        r = float(sqrt(c*c+1));
                                        psubd[i+1] = f*r;
                                        c *= (s = 1/r);
                                }
                                else {
                                        s = f/g;
                                        r = float(sqrt(s*s+1));
                                        psubd[i+1] = g*r;
                                        s *= (c = 1/r);
                                }
                                g = pdiag[i+1]-p;
                                r = (pdiag[i]-g)*s+2*b*c;
                                p = s*r;
                                pdiag[i+1] = g+p;
                                g = c*r-b;

                                for (int k = 0; k < n; k++) {
                                        f = pmat[k][i+1];
                                        pmat[k][i+1] = s*pmat[k][i]+c*f;
                                        pmat[k][i] = c*pmat[k][i]-s*f;
                                }
                        }
                        pdiag[ell] -= p;
                        psubd[ell] = g;
                        psubd[m] = 0;
                }
                if ( iter == eigen_maxiter ) {
                        Report(ql_exceeded);
                        return;
                }
        }
}
//---------------------------------------------------------------------------
void mgcEigen::
DecreasingSort (int n, float* eigval, float** eigvec)
{
        // sort eigenvalues in decreasing order, e[0] >= ... >= e[n-1]
        for (int i = 0, k; i <= n-2; i++) {
                // locate maximum eigenvalue
                float max = eigval[k=i];
                int j;
                for (j = i+1; j < n; j++)
                        if ( eigval[j] > max )
                                max = eigval[k=j];

                if ( k != i ) {
                        // swap eigenvalues
                        eigval[k] = eigval[i];
                        eigval[i] = max;

                        // swap eigenvectors
                        for (j = 0; j < n; j++) {
                                float tmp = eigvec[j][i];
                                eigvec[j][i] = eigvec[j][k];
                                eigvec[j][k] = tmp;
                        }
                }
        }
}
//---------------------------------------------------------------------------
void mgcEigen::
IncreasingSort (int n, float* eigval, float** eigvec)
{
        // sort eigenvalues in increasing order, e[0] <= ... <= e[n-1]
        for (int i = 0, k; i <= n-2; i++) {
                // locate minimum eigenvalue
                float min = eigval[k=i];
        int j;
                for (j = i+1; j < n; j++)
                        if ( eigval[j] < min )
                                min = eigval[k=j];

                if ( k != i ) {
                        // swap eigenvalues
                        eigval[k] = eigval[i];
                        eigval[i] = min;

                        // swap eigenvectors
                        for (j = 0; j < n; j++) {
                                float tmp = eigvec[j][i];
                                eigvec[j][i] = eigvec[j][k];
                                eigvec[j][k] = tmp;
                        }
                }
        }
}
//---------------------------------------------------------------------------
int mgcEigen::
Number (unsigned single_error)
{
        int result;
        for (result = -1; single_error; single_error >>= 1)
                result++;
        return result;
}
//---------------------------------------------------------------------------
void mgcEigen::
Report (unsigned single_error)
{
        if ( mgcEigen::verbose1 )
                cout << "mgcEigen: " << message[Number(single_error)] << endl;
        else  {
                ofstream ostr("eigen.err",ios::out|ios::app);
                ostr << "mgcEigen: " << message[Number(single_error)] << endl;
        }
        error |= single_error;
}
//---------------------------------------------------------------------------
void mgcEigen::
Report (ostream& ostr)
{
        for (unsigned single_error = 1; single_error; single_error <<= 1)
                if ( error & single_error )
                        ostr << "mgcEigen: " << message[Number(single_error)] << endl;

        error = 0;
}
//===========================================================================
// error handling
int mgcEigenD::verbose1 = 0;
unsigned mgcEigenD::error = 0;
const unsigned mgcEigenD::invalid_size      = 0x00000001;
const unsigned mgcEigenD::allocation_failed = 0x00000002;
const unsigned mgcEigenD::ql_exceeded       = 0x00000004;
const char* mgcEigenD::message[3] = {
        "invalid matrix size",
        "allocation failed",
        "QL algorithm - exceeded maximum iterations"
};
//---------------------------------------------------------------------------
mgcEigenD::
mgcEigenD (int _size)
{
        if ( (size = _size) <= 1 ) {
                Report(invalid_size);
                return;
        }
        if ( (mat = new double*[size]) == 0 ) {
                Report(allocation_failed);
                return;
        }
        for (int d = 0; d < size; d++)
                if ( (mat[d] = new double[size]) == 0 ) {
                        Report(allocation_failed);
                        return;
                }
        if ( (diag = new double[size]) == 0 ) {
                Report(allocation_failed);
                return;
        }
        if ( (subd = new double[size]) == 0 ) {
                Report(allocation_failed);
                return;
        }
}
//---------------------------------------------------------------------------
mgcEigenD::
~mgcEigenD ()
{
        delete[] subd;
        delete[] diag;
        for (int d = 0; d < size; d++)
                delete[] mat[d];
        delete[] mat;
}
//---------------------------------------------------------------------------
void mgcEigenD::
Tridiagonal2 (double** pmat, double* pdiag, double* psubd)
{
        // matrix is already tridiagonal

        pdiag[0] = pmat[0][0];
        pdiag[1] = pmat[1][1];
        psubd[0] = pmat[0][1];
        psubd[1] = 0;
        pmat[0][0] = 1;  pmat[0][1] = 0;
        pmat[1][0] = 0;  pmat[1][1] = 1;
}
//---------------------------------------------------------------------------
void mgcEigenD::
Tridiagonal3 (double** pmat, double* pdiag, double* psubd)
{
        double a = pmat[0][0], b = pmat[0][1], c = pmat[0][2],
        d = pmat[1][1], e = pmat[1][2], f = pmat[2][2];

        pdiag[0] = a;
        psubd[2] = 0;
        if ( c != 0 ) {
                double ell = sqrt(b*b+c*c);
                b /= ell;
                c /= ell;
                double q = 2*b*e+c*(f-d);
                pdiag[1] = d+c*q;
                pdiag[2] = f-c*q;
                psubd[0] = ell;
                psubd[1] = e-b*q;
                pmat[0][0] = 1; pmat[0][1] = 0; pmat[0][2] = 0;
                pmat[1][0] = 0; pmat[1][1] = b; pmat[1][2] = c;
                pmat[2][0] = 0; pmat[2][1] = c; pmat[2][2] = -b;
        }
        else {
                pdiag[1] = d;
                pdiag[2] = f;
                psubd[0] = b;
                psubd[1] = e;
                pmat[0][0] = 1; pmat[0][1] = 0; pmat[0][2] = 0;
                pmat[1][0] = 0; pmat[1][1] = 1; pmat[1][2] = 0;
                pmat[2][0] = 0; pmat[2][1] = 0; pmat[2][2] = 1;
        }
}
//---------------------------------------------------------------------------
void mgcEigenD::
Tridiagonal4 (double** pmat, double* pdiag, double* psubd)
{
        // save pmatrix M
        double
        a = pmat[0][0], b = pmat[0][1], c = pmat[0][2], d = pmat[0][3],
        e = pmat[1][1], f = pmat[1][2], g = pmat[1][3],
        h = pmat[2][2], i = pmat[2][3], j = pmat[3][3];

        pdiag[0] = a;
        psubd[3] = 0;

        pmat[0][0] = 1; pmat[0][1] = 0; pmat[0][2] = 0; pmat[0][3] = 0;
        pmat[1][0] = 0;
        pmat[2][0] = 0;
        pmat[3][0] = 0;

        if ( c != 0 || d != 0 ) {
                double q11, q12, q13;
                double q21, q22, q23;
                double q31, q32, q33;

                // build column Q1
                double len = sqrt(b*b+c*c+d*d);
                q11 = b/len;
                q21 = c/len;
                q31 = d/len;

                psubd[0] = len;

                // compute S*Q1
                double v0 = e*q11+f*q21+g*q31;
                double v1 = f*q11+h*q21+i*q31;
                double v2 = g*q11+i*q21+j*q31;

                pdiag[1] = q11*v0+q21*v1+q31*v2;

        // build column Q3 = Q1x(S*Q1)
                q13 = q21*v2-q31*v1;
                q23 = q31*v0-q11*v2;
                q33 = q11*v1-q21*v0;
                len = sqrt(q13*q13+q23*q23+q33*q33);
                if ( len > 0 ) {
                        q13 /= len;
                        q23 /= len;
                        q33 /= len;

                        // build column Q2 = Q3xQ1
                        q12 = q23*q31-q33*q21;
                        q22 = q33*q11-q13*q31;
                        q32 = q13*q21-q23*q11;

                        v0 = q12*e+q22*f+q32*g;
                        v1 = q12*f+q22*h+q32*i;
                        v2 = q12*g+q22*i+q32*j;
                        psubd[1] = q11*v0+q21*v1+q31*v2;
                        pdiag[2] = q12*v0+q22*v1+q32*v2;
                        psubd[2] = q13*v0+q23*v1+q33*v2;

                        v0 = q13*e+q23*f+q33*g;
                        v1 = q13*f+q23*h+q33*i;
                        v2 = q13*g+q23*i+q33*j;
                        pdiag[3] = q13*v0+q23*v1+q33*v2;
                }
                else {  // S*Q1 parallel to Q1, choose any valid Q2 and Q3
                        psubd[1] = 0;

                        len = q21*q21+q31*q31;
                        if ( len > 0 ) {
                                double tmp = q11-1;
                                q12 = -q21;
                                q22 = 1+tmp*q21*q21/len;
                q32 = tmp*q21*q31/len;

                                q13 = -q31;
                                q23 = q32;
                                q33 = 1+tmp*q31*q31/len;

                                v0 = q12*e+q22*f+q32*g;
                                v1 = q12*f+q22*h+q32*i;
                                v2 = q12*g+q22*i+q32*j;
                                pdiag[2] = q12*v0+q22*v1+q32*v2;
                                psubd[2] = q13*v0+q23*v1+q33*v2;

                                v0 = q13*e+q23*f+q33*g;
                                v1 = q13*f+q23*h+q33*i;
                                v2 = q13*g+q23*i+q33*j;
                                pdiag[3] = q13*v0+q23*v1+q33*v2;
                        }
                        else {  // Q1 = (+-1,0,0)
                                q12 = 0; q22 = 1; q32 = 0;
                                q13 = 0; q23 = 0; q33 = 1;

                                pdiag[2] = h;
                                pdiag[3] = j;
                                psubd[2] = i;
                        }
                }

                pmat[1][1] = q11; pmat[1][2] = q12; pmat[1][3] = q13;
                pmat[2][1] = q21; pmat[2][2] = q22; pmat[2][3] = q23;
                pmat[3][1] = q31; pmat[3][2] = q32; pmat[3][3] = q33;
        }
        else {
                pdiag[1] = e;
                psubd[0] = b;
                pmat[1][1] = 1;
                pmat[2][1] = 0;
                pmat[3][1] = 0; 

                if ( g != 0 ) {
                        double ell = sqrt(f*f+g*g);
                        f /= ell;
                        g /= ell;
                        double Q = 2*f*i+g*(j-h);

                        pdiag[2] = h+g*Q;
                        pdiag[3] = j-g*Q;
                        psubd[1] = ell;
                        psubd[2] = i-f*Q;
                        pmat[1][2] = 0;  pmat[1][3] = 0;
                        pmat[2][2] = f;  pmat[2][3] = g;
                        pmat[3][2] = g;  pmat[3][3] = -f;
                }
                else {
                        pdiag[2] = h;
                        pdiag[3] = j;
                        psubd[1] = f;
                        psubd[2] = i;
                        pmat[1][2] = 0;  pmat[1][3] = 0;
                        pmat[2][2] = 1;  pmat[2][3] = 0;
                        pmat[3][2] = 0;  pmat[3][3] = 1;
                }
        }
}
//---------------------------------------------------------------------------
void mgcEigenD::
TridiagonalN (int n, double** pmat, double* pdiag, double* psubd)
{
        int i, j, k, ell;

        for (i = n-1, ell = n-2; i >= 1; i--, ell--) {
                double h = 0, scale = 0;

                if ( ell > 0 ) {
                        for (k = 0; k <= ell; k++)
                                scale += fabs(pmat[i][k]);
                        if ( scale == 0 )
                                psubd[i] = pmat[i][ell];
                        else {
                                for (k = 0; k <= ell; k++) {
                                        pmat[i][k] /= scale;
                                        h += pmat[i][k]*pmat[i][k];
                                }
                                double f = pmat[i][ell];
                                double g = ( f > 0 ? -sqrt(h) : sqrt(h) );
                                psubd[i] = scale*g;
                                h -= f*g;
                                pmat[i][ell] = f-g;
                                f = 0;
                                for (j = 0; j <= ell; j++) {
                                        pmat[j][i] = pmat[i][j]/h;
                                        g = 0;
                                        for (k = 0; k <= j; k++)
                                                g += pmat[j][k]*pmat[i][k];
                                        for (k = j+1; k <= ell; k++)
                                                g += pmat[k][j]*pmat[i][k];
                                        psubd[j] = g/h;
                                        f += psubd[j]*pmat[i][j];
                                }
                                double hh = f/(h+h);
                                for (j = 0; j <= ell; j++) {
                                        f = pmat[i][j];
                                        psubd[j] = g = psubd[j] - hh*f;
                                        for (k = 0; k <= j; k++)
                                                pmat[j][k] -= f*psubd[k]+g*pmat[i][k];
                                }
            }
                }
                else
                        psubd[i] = pmat[i][ell];

                pdiag[i] = h;
        }

        pdiag[0] = psubd[0] = 0;
        for (i = 0, ell = -1; i <= n-1; i++, ell++) {
                if ( pdiag[i] ) {
                        for (j = 0; j <= ell; j++) {
                                double sum = 0;
                                for (k = 0; k <= ell; k++)
                                        sum += pmat[i][k]*pmat[k][j];
                                for (k = 0; k <= ell; k++)
                                        pmat[k][j] -= sum*pmat[k][i];
                        }
                }
                pdiag[i] = pmat[i][i];
                pmat[i][i] = 1;
                for (j = 0; j <= ell; j++)
                        pmat[j][i] = pmat[i][j] = 0;
        }

        // re-ordering if mgcEigenD::QLAlgorithm is used subsequently
        for (i = 1, ell = 0; i < n; i++, ell++)
                psubd[ell] = psubd[i];
        psubd[n-1] = 0;
}
//---------------------------------------------------------------------------
void mgcEigenD::
QLAlgorithm (int n, double* pdiag, double* psubd, double** pmat)
{
        const int eigen_maxiter = 30;

        for (int ell = 0; ell < n; ell++) {
                int iter;
                for (iter = 0; iter < eigen_maxiter; iter++) {
                        int m;
                        for (m = ell; m <= n-2; m++) {
                                double dd = fabs(pdiag[m])+fabs(pdiag[m+1]);
                                if ( (double)(fabs(psubd[m])+dd) == dd )
                                        break;
                        }
                        if ( m == ell )
                                break;

                        double g = (pdiag[ell+1]-pdiag[ell])/(2*psubd[ell]);
                        double r = sqrt(g*g+1);
                        if ( g < 0 )
                                g = pdiag[m]-pdiag[ell]+psubd[ell]/(g-r);
                        else
                                g = pdiag[m]-pdiag[ell]+psubd[ell]/(g+r);
                        double s = 1, c = 1, p = 0;
                        for (int i = m-1; i >= ell; i--) {
                                double f = s*psubd[i], b = c*psubd[i];
                                if ( fabs(f) >= fabs(g) ) {
                                        c = g/f;
                                        r = sqrt(c*c+1);
                                        psubd[i+1] = f*r;
                                        c *= (s = 1/r);
                                }
                                else {
                                        s = f/g;
                                        r = sqrt(s*s+1);
                                        psubd[i+1] = g*r;
                                        s *= (c = 1/r);
                                }
                                g = pdiag[i+1]-p;
                                r = (pdiag[i]-g)*s+2*b*c;
                                p = s*r;
                                pdiag[i+1] = g+p;
                                g = c*r-b;

                                for (int k = 0; k < n; k++) {
                                        f = pmat[k][i+1];
                                        pmat[k][i+1] = s*pmat[k][i]+c*f;
                                        pmat[k][i] = c*pmat[k][i]-s*f;
                                }
                        }
                        pdiag[ell] -= p;
                        psubd[ell] = g;
                        psubd[m] = 0;
                }
                if ( iter == eigen_maxiter ) {
                        Report(ql_exceeded);
                        return;
                }
        }
}
//---------------------------------------------------------------------------
void mgcEigenD::
DecreasingSort (int n, double* eigval, double** eigvec)
{
        // sort eigenvalues in decreasing order, e[0] >= ... >= e[n-1]
        for (int i = 0, k; i <= n-2; i++) {
                // locate maximum eigenvalue
                double max = eigval[k=i];
                int j;
                for (j = i+1; j < n; j++)
                        if ( eigval[j] > max )
                                max = eigval[k=j];

                if ( k != i ) {
                        // swap eigenvalues
                        eigval[k] = eigval[i];
                        eigval[i] = max;

                        // swap eigenvectors
                        for (j = 0; j < n; j++) {
                                double tmp = eigvec[j][i];
                                eigvec[j][i] = eigvec[j][k];
                                eigvec[j][k] = tmp;
                        }
                }
        }
}
//---------------------------------------------------------------------------
void mgcEigenD::
IncreasingSort (int n, double* eigval, double** eigvec)
{
        // sort eigenvalues in increasing order, e[0] <= ... <= e[n-1]
        for (int i = 0, k; i <= n-2; i++) {
                // locate minimum eigenvalue
                double min = eigval[k=i];
        int j;
                for (j = i+1; j < n; j++)
                        if ( eigval[j] < min )
                                min = eigval[k=j];

                if ( k != i ) {
                        // swap eigenvalues
                        eigval[k] = eigval[i];
                        eigval[i] = min;

                        // swap eigenvectors
                        for (j = 0; j < n; j++) {
                                double tmp = eigvec[j][i];
                                eigvec[j][i] = eigvec[j][k];
                                eigvec[j][k] = tmp;
                        }
                }
        }
}
//---------------------------------------------------------------------------
mgcEigenD& mgcEigenD::
Matrix (double** inmat)
{
        for (int row = 0; row < size; row++)
                for (int col = 0; col < size; col++)
                        mat[row][col] = inmat[row][col];
        return *this;
}
//---------------------------------------------------------------------------
void mgcEigenD::
EigenStuff3 ()
{
        Tridiagonal3(mat,diag,subd);
        QLAlgorithm(size,diag,subd,mat);
}
//---------------------------------------------------------------------------
int mgcEigenD::
Number (unsigned single_error)
{
        int result;
        for (result = -1; single_error; single_error >>= 1)
                result++;
        return result;
}
//---------------------------------------------------------------------------
void mgcEigenD::
Report (unsigned single_error)
{
        if ( mgcEigenD::verbose1 )
                cout << "mgcEigenD: " << message[Number(single_error)] << endl;
        else {
          ofstream ostr("eigen.err",ios::out|ios::app);
          ostr << "mgcEigenD: " << message[Number(single_error)] << endl;
        }
        error |= single_error;
}
//---------------------------------------------------------------------------
void mgcEigenD::
Report (ostream& ostr)
{
        for (unsigned single_error = 1; single_error; single_error <<= 1)
                if ( error & single_error )
                        ostr << "mgcEigenD: " << message[Number(single_error)] << endl;

        error = 0;
}
//===========================================================================

#ifdef EIGEN_TEST

int main ()
{
        mgcEigenD eig(3);

        eig.Matrix(0,0) = 2;  eig.Matrix(0,1) = 1;  eig.Matrix(0,2) = 1;
        eig.Matrix(1,0) = 1;  eig.Matrix(1,1) = 2;  eig.Matrix(1,2) = 1;
        eig.Matrix(2,0) = 1;  eig.Matrix(2,1) = 1;  eig.Matrix(2,2) = 2;

        eig.IncrSortEigenStuff3();

        cout.setf(ios::fixed);

        cout << "eigenvalues = " << endl;
        for (int row = 0; row < 3; row++)
                cout << eig.Eigenvalue(row) << ' ';
        cout << endl;

        cout << "eigenvectors = " << endl;
        for (row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++)
                        cout << eig.Eigenvector(row,col) << ' ';
                cout << endl;
        }

        // eigenvalues =
        //    1.000000 1.000000 4.000000
        // eigenvectors =
        //    0.411953  0.704955 0.577350
        //    0.404533 -0.709239 0.577350
        //   -0.816485  0.004284 0.577350

        return 0;
}
#endif
#endif

// FILE: tri3tri3.cpp
//---------------------------------------------------------------------------
// MinTriangleTriangle
//
// The quadratic form representing the squared distance between two planes
// never has an isolated global minimum.  The generic case is for it to
// have an entire line of zeros (the line of intersection of the two
// planes).  Therefore, it is sufficient to compare edges of each triangle
// with the entire other triangle, looking for the minimum distance.
//---------------------------------------------------------------------------

// NOTE:  This code is not fully optimized as is the code in pt3lin3,
// pt3tri3, and lin3lin3.

//---------------------------------------------------------------------------
double 
AnalyticGeometryTool::MinTriangleTriangle (Triangle3& tri0, 
                                           Triangle3& tri1,
                                           double& s, double& t, double& u, 
                                           double& v)
{
    double s0, t0, u0, v0, min, min0;
    Line3 line;

    // compare edges of tri0 against all of tri1
    line.b = tri0.b;
    line.m = tri0.e0;
    min = MinLineSegmentTriangle(line,tri1,s,u,v);
    t = 0;

    line.m = tri0.e1;
    min0 = MinLineSegmentTriangle(line,tri1,t0,u0,v0);
    s0 = 0;
    if ( min0 < min )
    {
        min = min0;
        s = s0;
        t = t0;
        u = u0;
        v = v0;
    }

    line.b = line.b + tri0.e0;
    line.m = line.m - tri0.e0;
    min0 = MinLineSegmentTriangle(line,tri1,t0,u0,v0);
    s0 = 1-t0;
    if ( min0 < min )
    {
        min = min0;
        s = s0;
        t = t0;
        u = u0;
        v = v0;
    }

    // compare edges of tri1 against all of tri0
    line.b = tri1.b;
    line.m = tri1.e0;
    min0 = MinLineSegmentTriangle(line,tri0,u0,s0,t0);
    v0 = 0;
    if ( min0 < min )
    {
        min = min0;
        s = s0;
        t = t0;
        u = u0;
        v = v0;
    }

    line.m = tri1.e1;
    min0 = MinLineSegmentTriangle(line,tri0,v0,s0,t0);
    u0 = 0;
    if ( min0 < min )
    {
        min = min0;
        s = s0;
        t = t0;
        u = u0;
        v = v0;
    }

    line.b = line.b + tri1.e0;
    line.m = line.m - tri1.e0;
    min0 = MinLineSegmentTriangle(line,tri0,v0,s0,t0);
    u0 = 1-v0;
    if ( min0 < min )
    {
        min = min0;
        s = s0;
        t = t0;
        u = u0;
        v = v0;
    }

    return min;
}

//---------------------------------------------------------------------------

#ifdef TRI3TRI3_TEST

//#define RAND (2.0f*rand()/double(RAND_MAX)-1)

//ofstream ostr("data.txt");

void main ()
{
    Triangle3 tri0, tri1;
    Point3 p, q, diff;
    double u, v, s, t, u0, v0, s0, t0, min0, min1, dist;

    double maxdiff = 0;

    for (int i = 0; i < 128; i++)
    {
        tri0.b.x = RAND;
        tri0.b.y = RAND;
        tri0.b.z = RAND;
        tri0.e0.x = RAND;
        tri0.e0.y = RAND;
        tri0.e0.z = RAND;
        tri0.e1.x = RAND;
        tri0.e1.y = RAND;
        tri0.e1.z = RAND;

        tri1.b.x = RAND;
        tri1.b.y = RAND;
        tri1.b.z = RAND;
        tri1.e0.x = RAND;
        tri1.e0.y = RAND;
        tri1.e0.z = RAND;
        tri1.e1.x = RAND;
        tri1.e1.y = RAND;
        tri1.e1.z = RAND;

        min0 = FLT_MAX;
        int max = 32;
        for (int w = 0; w <= max; w++)
        {
            t0 = w/double(max);
            for (int z = 0; w+z <= max; z++)
            {
                s0 = z/double(max);
                p = tri0.b+s0*tri0.e0+t0*tri0.e1;
                for (int y = 0; y <= max; y++)
                {
                    v0 = y/double(max);
                    for (int x = 0; x+y <= max; x++)
                    {
                        u0 = x/double(max);
                        q = tri1.b+u0*tri1.e0+v0*tri1.e1;
                        
                        diff = p-q;
                        dist = Length(diff);
                        if ( dist < min0 )
                        {
                            min0 = dist;
                            s = s0;
                            t = t0;
                            u = u0;
                            v = v0;
                        }
                    }
                }
            }
        }
        ostr << "i = " << i << endl;
        ostr << "sampled = " << s << ' ' << t << ' ' << u << ' '
             << v << ' ' << min0 << endl;

        min1 = MinTriangleTriangle(tri0,tri1,s,t,u,v);
        ostr << "analytic = " << s << ' ' << t << ' ' << u << ' '
             << v << ' ' << min1 << endl;

        ostr << "diff = " << min1-min0 << endl;

        if ( min1-min0 > maxdiff )
            maxdiff = min1-min0;

        ostr << endl;
    }

    ostr << "max diff = " << maxdiff << endl;
}
#endif

//---------------------------------------------------------------------------
double 
AnalyticGeometryTool::MinLineSegmentLineSegment (const Line3& seg0, const Line3& seg1,
                                                 double& s, double& t)
{
    Point3 diff = seg0.b - seg1.b;
    double A = Dot(seg0.m,seg0.m);
    double B = -Dot(seg0.m,seg1.m);
    double C = Dot(seg1.m,seg1.m);
    double D = Dot(seg0.m,diff);
    double E;  // -Dot(seg1.m,diff), defer until needed
    double F = Dot(diff,diff);
    double det = Abs(A*C-B*B);  // A*C-B*B = |Cross(M0,M1)|^2 >= 0

    double tmp;

    if ( det >= par_tolerance )
    {
        // line segments are not parallel
        E = -Dot(seg1.m,diff);
        s = B*E-C*D;
        t = B*D-A*E;
        
        if ( s >= 0 )
        {
            if ( s <= det )
            {
                if ( t >= 0 )
                {
                    if ( t <= det )  // region 0 (interior)
                    {
                        // minimum at two interior points of 3D lines
                        double invDet = 1.0f/det;
                        s *= invDet;
                        t *= invDet;
                        return DIST(s*(A*s+B*t+2*D)+t*(B*s+C*t+2*E)+F);
                    }
                    else  // region 3 (side)
                    {
                        t = 1;
                        tmp = B+D;
                        if ( tmp >= 0 )
                        {
                            s = 0;
                            return DIST(C+2*E+F);
                        }
                        else if ( -tmp >= A )
                        {
                            s = 1;
                            return DIST(A+C+F+2*(E+tmp));
                        }
                        else
                        {
                            s = -tmp/A;
                            return DIST(tmp*s+C+2*E+F);
                        }
                    }
                }
                else  // region 7 (side)
                {
                    t = 0;
                    if ( D >= 0 )
                    {
                        s = 0;
                        return DIST(F);
                    }
                    else if ( -D >= A )
                    {
                        s = 1;
                        return DIST(A+2*D+F);
                    }
                    else
                    {
                        s = -D/A;
                        return DIST(D*s+F);
                    }
                }
            }
            else
            {
                if ( t >= 0 )
                {
                    if ( t <= det )  // region 1 (side)
                    {
                        s = 1;
                        tmp = B+E;
                        if ( tmp >= 0 )
                        {
                            t = 0;
                            return DIST(A+2*D+F);
                        }
                        else if ( -tmp >= C )
                        {
                            t = 1;
                            return DIST(A+C+F+2*(D+tmp));
                        }
                        else
                        {
                            t = -tmp/C;
                            return DIST(tmp*t+A+2*D+F);
                        }
                    }
                    else  // region 2 (corner)
                    {
                        tmp = B+D;
                        if ( -tmp <= A )
                        {
                            t = 1;
                            if ( tmp >= 0 )
                            {
                                s = 0;
                                return DIST(C+2*E+F);
                            }
                            else
                            {
                                 s = -tmp/A;
                                 return DIST(tmp*s+C+2*E+F);
                            }
                        }
                        else
                        {
                            s = 1;
                            tmp = B+E;
                            if ( tmp >= 0 )
                            {
                                t = 0;
                                return DIST(A+2*D+F);
                            }
                            else if ( -tmp >= C )
                            {
                                t = 1;
                                return DIST(A+C+F+2*(D+tmp));
                            }
                            else
                            {
                                t = -tmp/C;
                                return DIST(tmp*t+A+2*D+F);
                            }
                        }
                    }
                }
                else  // region 8 (corner)
                {
                    if ( -D < A )
                    {
                        t = 0;
                        if ( D >= 0 )
                        {
                            s = 0;
                            return DIST(F);
                        }
                        else
                        {
                            s = -D/A;
                            return DIST(D*s+F);
                        }
                    }
                    else
                    {
                        s = 1;
                        tmp = B+E;
                        if ( tmp >= 0 )
                        {
                            t = 0;
                            return DIST(A+2*D+F);
                        }
                        else if ( -tmp >= C )
                        {
                            t = 1;
                            return DIST(A+C+F+2*(D+tmp));
                        }
                        else
                        {
                            t = -tmp/C;
                            return DIST(tmp*t+A+2*D+F);
                        }
                    }
                }
            }
        }
        else 
        {
            if ( t >= 0 )
            {
                if ( t <= det )  // region 5 (side)
                {
                    s = 0;
                    if ( E >= 0 )
                    {
                        t = 0;
                        return DIST(F);
                    }
                    else if ( -E >= C )
                    {
                        t = 1;
                        return DIST(C+2*E+F);
                    }
                    else
                    {
                        t = -E/C;
                        return DIST(E*t+F);
                    }
                }
                else  // region 4 (corner)
                {
                    tmp = B+D;
                    if ( tmp < 0 )
                    {
                        t = 1;
                        if ( -tmp >= A )
                        {
                            s = 1;
                            return DIST(A+C+F+2*(E+tmp));
                        }
                        else
                        {
                            s = -tmp/A;
                            return DIST(tmp*s+C+2*E+F);
                        }
                    }
                    else
                    {
                        s = 0;
                        if ( E >= 0 )
                        {
                            t = 0;
                            return DIST(F);
                        }
                        else if ( -E >= C )
                        {
                            t = 1;
                            return DIST(C+2*E+F);
                        }
                        else
                        {
                            t = -E/C;
                            return DIST(E*t+F);
                        }
                    }
                }
            }
            else   // region 6 (corner)
            {
                if ( D < 0 )
                {
                    t = 0;
                    if ( -D >= A )
                    {
                        s = 1;
                        return DIST(A+2*D+F);
                    }
                    else
                    {
                        s = -D/A;
                        return DIST(D*s+F);
                    }
                }
                else
                {
                    s = 0;
                    if ( E >= 0 )
                    {
                        t = 0;
                        return DIST(F);
                    }
                    else if ( -E >= C )
                    {
                        t = 1;
                        return DIST(C+2*E+F);
                    }
                    else
                    {
                        t = -E/C;
                        return DIST(E*t+F);
                    }
                }
            }
        }
    }
    else
    {
        // line segments are parallel
        if ( B > 0 )
        {
            // direction vectors form an obtuse angle
            if ( D >= 0 )
            {
                s = 0;
                t = 0;
                return DIST(F);
            }
            else if ( -D <= A )
            {
                s = -D/A;
                t = 0;
                return DIST(D*s+F);
            }
            else
            {
                E = -Dot(seg1.m,diff);
                s = 1;
                tmp = A+D;
                if ( -tmp >= B )
                {
                    t = 1;
                    return DIST(A+C+F+2*(B+D+E));
                }
                else
                {
                    t = -tmp/B;
                    return DIST(A+2*D+F+t*(C*t+2*(B+E)));
                }
            }
        }
        else
        {
            // direction vectors form an acute angle
            if ( -D >= A )
            {
                s = 1;
                t = 0;
                return DIST(A+2*D+F);
            }
            else if ( D <= 0 )
            {
                s = -D/A;
                t = 0;
                return DIST(D*s+F);
            }
            else
            {
                E = -Dot(seg1.m,diff);
                s = 0;
                if ( D >= -B )
                {
                    t = 1;
                    return DIST(C+2*E+F);
                }
                else
                {
                    t = -D/B;
                    return DIST(F+t*(2*E+C*t));
                }
            }
        }
    }
}
//---------------------------------------------------------------------------

#ifdef LIN3LIN3_TEST

//#include <stdlib.h>
//#include <fstream.h>

//#define RAND (2.0f*rand()/double(RAND_MAX)-1)

//ofstream ostr("data.txt");

void TestSegSeg ()
{
    Line3 seg0, seg1;
    Point3 p0, p1, diff;
    double s, t, s0, t0, min0, min1, dist;
    double maxDiff = 0.0f;

    for (int i = 0; i < 128; i++)
    {
        seg0.b.x = RAND;
        seg0.b.y = RAND;
        seg0.b.z = RAND;

        seg0.m.x = RAND;
        seg0.m.y = RAND;
        seg0.m.z = RAND;

        seg1.b.x = RAND;
        seg1.b.y = RAND;
        seg1.b.z = RAND;

        if ( i % 2 )
        {
            // non-parallel line segments
            seg1.m.x = RAND;
            seg1.m.y = RAND;
            seg1.m.z = RAND;
        }
        else
        {
            // parallel line segments
            double scale = RAND;
            seg1.m = scale*seg0.m;
        }

        min0 = FLT_MAX;
        int ymax = 128, xmax = 128;
        for (int y = 0; y < ymax; y++)
        {
            s0 = y/double(ymax-1);
            p0 = seg0.b+s0*seg0.m;
            for (int x = 0; x < xmax; x++)
            {
                t0 = x/double(xmax-1);
                p1 = seg1.b+t0*seg1.m;

                diff = p1-p0;
                dist = Length(diff);
                if ( dist < min0 )
                {
                    min0 = dist;
                    s = s0;
                    t = t0;
                }
            }
        }
        ostr << "sampled = " << s << ' ' << t << ' ' << min0 << endl;

        min1 = MinLineSegmentLineSegment(seg0,seg1,s,t);
        ostr << "analytic = " << s << ' ' << t << ' ' << min1 << endl;

        double compDiff = min1-min0;
        ostr << "diff = " << compDiff << endl;
        if ( compDiff > maxDiff )
            maxDiff = compDiff;

        ostr << endl;
    }
    ostr << "max diff = " << maxDiff << endl;
}
#endif

// FILE: pt3tri3.cpp
//---------------------------------------------------------------------------
double 
AnalyticGeometryTool::MinPointTriangle (const Point3& p, const Triangle3& tri,
                                        double& s, double& t)
{
    Point3 diff = tri.b - p;
    double A = Dot(tri.e0,tri.e0);
    double B = Dot(tri.e0,tri.e1);
    double C = Dot(tri.e1,tri.e1);
    double D = Dot(tri.e0,diff);
    double E = Dot(tri.e1,diff);
    double F = Dot(diff,diff);
    double det = Abs(A*C-B*B);  // A*C-B*B = |Cross(e0,e1)|^2 >= 0
    // A-2*B+C = Dot(e0,e0)-2*Dot(e0,e1)+Dot(e1,e1) = |e0-e1|^2 > 0

    s = B*E-C*D;
    t = B*D-A*E;

    if ( s+t <= det )
    {
        if ( s < 0 )
        {
            if ( t < 0 )  // region 4
            {
                if ( D < 0 )
                {
                    t = 0;
                    if ( -D >= A )
                    {
                        s = 1;
                        return DIST(A+2*D+F);
                    }
                    else
                    {
                        s = -D/A;
                        return DIST(D*s+F);
                    }
                }
                else
                {
                    s = 0;
                    if ( E >= 0 )
                    {
                        t = 0;
                        return DIST(F);
                    }
                    else if ( -E >= C )
                    {
                        t = 1;
                        return DIST(C+2*E+F);
                    }
                    else
                    {
                        t = -E/C;
                        return DIST(E*t+F);
                    }
                }
            }
            else  // region 3
            {
                s = 0;
                if ( E >= 0 )
                {
                    t = 0;
                    return DIST(F);
                }
                else if ( -E >= C )
                {
                    t = 1;
                    return DIST(C+2*E+F);
                }
                else
                {
                    t = -E/C;
                    return DIST(E*t+F);
                }
            }
        }
        else if ( t < 0 )  // region 5
        {
            t = 0;
            if ( D >= 0 )
            {
                s = 0;
                return DIST(F);
            }
            else if ( -D >= A )
            {
                s = 1;
                return DIST(A+2*D+F);
            }
            else
            {
                s = -D/A;
                return DIST(D*s+F);
            }
        }
        else  // region 0
        {
            // minimum at interior point
          if( det == 0.0 )
          {
            //PRINT_WARNING( "Found zero determinant\n" );
            return CUBIT_DBL_MAX;
          }

            double invDet = 1.0f/det;
            s *= invDet;
            t *= invDet;
            return DIST(s*(A*s+B*t+2*D)+t*(B*s+C*t+2*E)+F);
        }
    }
    else
    {
        double tmp0, tmp1, numer, denom;

        if ( s < 0 )  // region 2
        {
            tmp0 = B+D;
            tmp1 = C+E;
            if ( tmp1 > tmp0 )
            {
                numer = tmp1 - tmp0;
                denom = A-2*B+C;
                if ( numer >= denom )
                {
                    s = 1;
                    t = 0;
                    return DIST(A+2*D+F);
                }
                else
                {
                    s = numer/denom;
                    t = 1-s;
                    return DIST(s*(A*s+B*t+2*D)+t*(B*s+C*t+2*E)+F);
                }
            }
            else
            {
                s = 0;
                if ( tmp1 <= 0 )
                {
                    t = 1;
                    return DIST(C+2*E+F);
                }
                else if ( E >= 0 )
                {
                    t = 0;
                    return DIST(F);
                }
                else
                {
                    t = -E/C;
                    return DIST(E*t+F);
                }
            }
        }
        else if ( t < 0 )  // region 6
        {
            tmp0 = B+E;
            tmp1 = A+D;
            if ( tmp1 > tmp0 )
            {
                numer = tmp1 - tmp0;
                denom = A-2*B+C;
                if ( numer >= denom )
                {
                    t = 1;
                    s = 0;
                    return DIST(C+2*E+F);
                }
                else
                {
                    t = numer/denom;
                    s = 1-t;
                    return DIST(s*(A*s+B*t+2*D)+t*(B*s+C*t+2*E)+F);
                }
            }
            else
            {
                t = 0;
                if ( tmp1 <= 0 )
                {
                    s = 1;
                    return DIST(A+2*D+F);
                }
                else if ( D >= 0 )
                {
                    s = 0;
                    return DIST(F);
                }
                else
                {
                    s = -D/A;
                    return DIST(D*s+F);
                }
            }
        }
        else  // region 1
        {
            numer = C+E-B-D;
            if ( numer <= 0 )
            {
                s = 0;
                t = 1;
                return DIST(C+2*E+F);
            }
            else
            {
                denom = A-2*B+C;
                if ( numer >= denom )
                {
                    s = 1;
                    t = 0;
                    return DIST(A+2*D+F);
                }
                else
                {
                    s = numer/denom;
                    t = 1-s;
                    return DIST(s*(A*s+B*t+2*D)+t*(B*s+C*t+2*E)+F);
                }
            }
        }
    }
}
//---------------------------------------------------------------------------

#ifdef PT3TRI3_TEST

//#include <float.h>
//#include <fstream.h>
//#include <stdlib.h>

//ofstream ostr("data.txt");

//#define RAND (2.0f*rand()/double(RAND_MAX)-1)

void main ()
{
    Triangle3 tri;
    Point3 p, q, diff;
    double s, t, s0, t0, min0, min1, dist;

    double maxdiff = 0;

    for (int i = 0; i < 128; i++)
    {
        tri.b.x = RAND;
        tri.b.y = RAND;
        tri.b.z = RAND;
        tri.e0.x = RAND;
        tri.e0.y = RAND;
        tri.e0.z = RAND;
        tri.e1.x = RAND;
        tri.e1.y = RAND;
        tri.e1.z = RAND;

        p.x = RAND;
        p.y = RAND;
        p.z = RAND;

        min0 = FLT_MAX;
        int max = 128;
        for (int y = 0; y <= max; y++)
        {
            s0 = y/double(max);
            for (int x = 0; x+y <= max; x++)
            {
                t0 = x/double(max);
                q = tri.b+s0*tri.e0+t0*tri.e1;

                diff = p-q;
                dist = Length(diff);
                if ( dist < min0 )
                {
                    min0 = dist;
                    s = s0;
                    t = t0;
                }
            }
        }
        ostr << "sampled = " << s << ' ' << t << ' ' << min0 << endl;

        min1 = MinPointTriangle(p,tri,s,t);
        ostr << "analytic = " << s << ' ' << t << ' ' << min1 << endl;

        ostr << "diff = " << min1-min0 << endl;

        if ( min1-min0 > maxdiff )
            maxdiff = min1-min0;

        ostr << endl;
    }
    ostr << "max diff = " << maxdiff << endl;
}
#endif

// FILE: lin3tri3.cpp
//---------------------------------------------------------------------------
// This code computes the closest points of a line L(r) = a0+r*a1 and a
// triangle Tri(s,t) = b0+s*b1+t*b2, where 0 <= r <= 1 and 0 <= s <= 1,
// 0 <= t <= 1, and 0 <= s+t <= 1.
//
// In calculus terms, the goal is to minimize the squared-distance function
// Q(r,s,t) = Dot(L(r)-Tri(s,t),L(r)-Tri(s,t)) over the prism domain
// 0 <= r <= 1, 0 <= s <= 1, 0 <= t <= 1, 0 <= s+t <= 1.  This function is
//
//   Q(r,s,t) = [r s t] A [r s t] + 2 B [r s t] + C
//
// where A is a 3x3 symmetric matrix, B is a 3x1 column vector, and C is
// a constant.  The entries of A and B are various dot products derived
// from the vectors a0, a1, b0, b1, and b2.
//
// The analysis is similar to that of MinPointTriangle.  The (s,t) domain
// is partitioned similarly, but the prism is obtained by extruding the
// domain in the r-direction.  The three cases for r are: r < 0, 0 <= r <= 1,
// and r > 1.  For each case the partitioning in (s,t) is identical to that
// of MinPointTriangle.  The region numbers for r < 0 are appended with an
// 'm' and the region numbers for r > 1 are appended with a 'p'.
//---------------------------------------------------------------------------

// NOTE:  This code is not fully optimized as is the code in pt3lin3,
// pt3tri3, and lin3lin3.

//---------------------------------------------------------------------------
double 
AnalyticGeometryTool::MinLineSegmentTriangle (const Line3& seg, const Triangle3& tri,
                                              double& r, double& s, double& t)
{
    Point3 diff = tri.b - seg.b;
    double A00 = Dot(seg.m,seg.m);
    double A01 = -Dot(seg.m,tri.e0);
    double A02 = -Dot(seg.m,tri.e1);
    double A11 = Dot(tri.e0,tri.e0);
    double A12 = Dot(tri.e0,tri.e1);
    double A22 = Dot(tri.e1,tri.e1);
    double B0  = -Dot(diff,seg.m);
    double B1  = Dot(diff,tri.e0);
    double B2  = Dot(diff,tri.e1);
    double cof00 = A11*A22-A12*A12;
    double cof01 = A02*A12-A01*A22;
    double cof02 = A01*A12-A02*A11;
    double det = A00*cof00+A01*cof01+A02*cof02;

    Line3 triseg;
    Point3 pt;
    double min, min0, r0, s0, t0;

    if ( Abs(det) >= par_tolerance )
    {
        double cof11 = A00*A22-A02*A02;
        double cof12 = A02*A01-A00*A12;
        double cof22 = A00*A11-A01*A01;
        double invDet = 1.0f/det;
        double rhs0 = -B0*invDet;
        double rhs1 = -B1*invDet;
        double rhs2 = -B2*invDet;

        r = cof00*rhs0+cof01*rhs1+cof02*rhs2;
        s = cof01*rhs0+cof11*rhs1+cof12*rhs2;
        t = cof02*rhs0+cof12*rhs1+cof22*rhs2;

        if ( r < 0 )
        {
            if ( s+t <= 1 )
            {
                if ( s < 0 )
                {
                    if ( t < 0 )  // region 4m
                    {
                        // min on face s=0 or t=0 or r=0
                        triseg.b = tri.b;
                        triseg.m = tri.e1;
                        min = MinLineSegmentLineSegment(seg,triseg,r,t);
                        s = 0;
                        triseg.b = tri.b;
                        triseg.m = tri.e0;
                        min0 = MinLineSegmentLineSegment(seg,triseg,r0,s0);
                        t0 = 0;
                        if ( min0 < min )
                        {
                            min = min0;
                            r = r0;
                            s = s0;
                            t = t0;
                        }
                        min0 = MinPointTriangle(seg.b,tri,s0,t0);
                        r0 = 0;
                        if ( min0 < min )
                        {
                            min = min0;
                            r = r0;
                            s = s0;
                            t = t0;
                        }
                    }
                    else  // region 3m
                    {
                        // min on face s=0 or r=0
                        triseg.b = tri.b;
                        triseg.m = tri.e1;
                        min = MinLineSegmentLineSegment(seg,triseg,r,t);
                        s = 0;
                        min0 = MinPointTriangle(seg.b,tri,s0,t0);
                        r0 = 0;
                        if ( min0 < min )
                        {
                            min = min0;
                            r = r0;
                            s = s0;
                            t = t0;
                        }
                    }
                }
                else if ( t < 0 )  // region 5m
                {
                    // min on face t=0 or r=0
                    triseg.b = tri.b;
                    triseg.m = tri.e0;
                    min = MinLineSegmentLineSegment(seg,triseg,r,s);
                    t = 0;
                    min0 = MinPointTriangle(seg.b,tri,s0,t0);
                    r0 = 0;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                }
                else  // region 0m
                {
                    // min face on r=0
                    min = MinPointTriangle(seg.b,tri,s,t);
                    r = 0;
                }
            }
            else
            {
                if ( s < 0 )  // region 2m
                {
                    // min on face s=0 or s+t=1 or r=0
                    triseg.b = tri.b;
                    triseg.m = tri.e1;
                    min = MinLineSegmentLineSegment(seg,triseg,r,t);
                    s = 0;
                    triseg.b = tri.b+tri.e0;
                    triseg.m = tri.e1-tri.e0;
                    min0 = MinLineSegmentLineSegment(seg,triseg,r0,t0);
                    s0 = 1-t0;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                    min0 = MinPointTriangle(seg.b,tri,s0,t0);
                    r0 = 0;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                }
                else if ( t < 0 )  // region 6m
                {
                    // min on face t=0 or s+t=1 or r=0
                    triseg.b = tri.b;
                    triseg.m = tri.e0;
                    min = MinLineSegmentLineSegment(seg,triseg,r,s);
                    t = 0;
                    triseg.b = tri.b+tri.e0;
                    triseg.m = tri.e1-tri.e0;
                    min0 = MinLineSegmentLineSegment(seg,triseg,r0,t0);
                    s0 = 1-t0;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                    min0 = MinPointTriangle(seg.b,tri,s0,t0);
                    r0 = 0;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                }
                else  // region 1m
                {
                    // min on face s+t=1 or r=0
                    triseg.b = tri.b+tri.e0;
                    triseg.m = tri.e1-tri.e0;
                    min = MinLineSegmentLineSegment(seg,triseg,r,t);
                    s = 1-t;
                    min0 = MinPointTriangle(seg.b,tri,s0,t0);
                    r0 = 0;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                }
            }
        }
        else if ( r <= 1 )
        {
            if ( s+t <= 1 )
            {
                if ( s < 0 )
                {
                    if ( t < 0 )  // region 4
                    {
                        // min on face s=0 or t=0
                        triseg.b = tri.b;
                        triseg.m = tri.e1;
                        min = MinLineSegmentLineSegment(seg,triseg,r,t);
                        s = 0;
                        triseg.b = tri.b;
                        triseg.m = tri.e0;
                        min0 = MinLineSegmentLineSegment(seg,triseg,r0,s0);
                        t0 = 0;
                        if ( min0 < min )
                        {
                            min = min0;
                            r = r0;
                            s = s0;
                            t = t0;
                        }
                    }
                    else  // region 3
                    {
                        // min on face s=0
                        triseg.b = tri.b;
                        triseg.m = tri.e1;
                        min = MinLineSegmentLineSegment(seg,triseg,r,t);
                        s = 0;
                    }
                }
                else if ( t < 0 )  // region 5
                {
                    // min on face t=0
                    triseg.b = tri.b;
                    triseg.m = tri.e0;
                    min = MinLineSegmentLineSegment(seg,triseg,r,s);
                    t = 0;
                }
                else  // region 0
                {
                    // global minimum is interior, done
                    min = Sqrt(Abs(r*(A00*r+A01*s+A02*t+2.0f*B0)
                          +s*(A01*r+A11*s+A12*t+2.0f*B1)
                          +t*(A02*r+A12*s+A22*t+2.0f*B2)
                          +Dot(diff,diff)));
                }
            }
            else
            {
                if ( s < 0 )  // region 2
                {
                    // min on face s=0 or s+t=1
                    triseg.b = tri.b;
                    triseg.m = tri.e1;
                    min = MinLineSegmentLineSegment(seg,triseg,r,t);
                    s = 0;
                    triseg.b = tri.b+tri.e0;
                    triseg.m = tri.e1-tri.e0;
                    min0 = MinLineSegmentLineSegment(seg,triseg,r0,t0);
                    s0 = 1-t0;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                }
                else if ( t < 0 )  // region 6
                {
                    // min on face t=0 or s+t=1
                    triseg.b = tri.b;
                    triseg.m = tri.e0;
                    min = MinLineSegmentLineSegment(seg,triseg,r,s);
                    t = 0;
                    triseg.b = tri.b+tri.e0;
                    triseg.m = tri.e1-tri.e0;
                    min0 = MinLineSegmentLineSegment(seg,triseg,r0,t0);
                    s0 = 1-t0;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                }
                else  // region 1
                {
                    // min on face s+t=1
                    triseg.b = tri.b+tri.e0;
                    triseg.m = tri.e1-tri.e0;
                    min = MinLineSegmentLineSegment(seg,triseg,r,t);
                    s = 1-t;
                }
            }
        }
        else  // r > 1
        {
            if ( s+t <= 1 )
            {
                if ( s < 0 )
                {
                    if ( t < 0 )  // region 4p
                    {
                        // min on face s=0 or t=0 or r=1
                        triseg.b = tri.b;
                        triseg.m = tri.e1;
                        min = MinLineSegmentLineSegment(seg,triseg,r,t);
                        s = 0;
                        triseg.b = tri.b;
                        triseg.m = tri.e0;
                        min0 = MinLineSegmentLineSegment(seg,triseg,r0,s0);
                        t0 = 0;
                        if ( min0 < min )
                        {
                            min = min0;
                            r = r0;
                            s = s0;
                            t = t0;
                        }
                        pt = seg.b+seg.m;
                        min0 = MinPointTriangle(pt,tri,s0,t0);
                        r0 = 1;
                        if ( min0 < min )
                        {
                            min = min0;
                            r = r0;
                            s = s0;
                            t = t0;
                        }
                    }
                    else  // region 3p
                    {
                        // min on face s=0 or r=1
                        triseg.b = tri.b;
                        triseg.m = tri.e1;
                        min = MinLineSegmentLineSegment(seg,triseg,r,t);
                        s = 0;
                        pt = seg.b+seg.m;
                        min0 = MinPointTriangle(pt,tri,s0,t0);
                        r0 = 1;
                        if ( min0 < min )
                        {
                            min = min0;
                            r = r0;
                            s = s0;
                            t = t0;
                        }
                    }
                }
                else if ( t < 0 )  // region 5p
                {
                    // min on face t=0 or r=1
                    triseg.b = tri.b;
                    triseg.m = tri.e0;
                    min = MinLineSegmentLineSegment(seg,triseg,r,s);
                    t = 0;
                    pt = seg.b+seg.m;
                    min0 = MinPointTriangle(pt,tri,s0,t0);
                    r0 = 1;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                }
                else  // region 0p
                {
                    // min face on r=1
                    pt = seg.b+seg.m;
                    min = MinPointTriangle(pt,tri,s,t);
                    r = 1;
                }
            }
            else
            {
                if ( s < 0 )  // region 2p
                {
                    // min on face s=0 or s+t=1 or r=1
                    triseg.b = tri.b;
                    triseg.m = tri.e1;
                    min = MinLineSegmentLineSegment(seg,triseg,r,t);
                    s = 0;
                    triseg.b = tri.b+tri.e0;
                    triseg.m = tri.e1-tri.e0;
                    min0 = MinLineSegmentLineSegment(seg,triseg,r0,t0);
                    s0 = 1-t0;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                    pt = seg.b+seg.m;
                    min0 = MinPointTriangle(pt,tri,s0,t0);
                    r0 = 1;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                }
                else if ( t < 0 )  // region 6p
                {
                    // min on face t=0 or s+t=1 or r=1
                    triseg.b = tri.b;
                    triseg.m = tri.e0;
                    min = MinLineSegmentLineSegment(seg,triseg,r,s);
                    t = 0;
                    triseg.b = tri.b+tri.e0;
                    triseg.m = tri.e1-tri.e0;
                    min0 = MinLineSegmentLineSegment(seg,triseg,r0,t0);
                    s0 = 1-t0;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                    pt = seg.b+seg.m;
                    min0 = MinPointTriangle(pt,tri,s0,t0);
                    r0 = 1;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                }
                else  // region 1p
                {
                    // min on face s+t=1 or r=1
                    triseg.b = tri.b+tri.e0;
                    triseg.m = tri.e1-tri.e0;
                    min = MinLineSegmentLineSegment(seg,triseg,r,t);
                    s = 1-t;
                    pt = seg.b+seg.m;
                    min0 = MinPointTriangle(pt,tri,s0,t0);
                    r0 = 1;
                    if ( min0 < min )
                    {
                        min = min0;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                }
            }
        }
    }
    else
    {
        // line and triangle are parallel
        triseg.b = tri.b;
        triseg.m = tri.e0;
        min = MinLineSegmentLineSegment(seg,triseg,r,s);
        t = 0;

        triseg.m = tri.e1;
        min0 = MinLineSegmentLineSegment(seg,triseg,r0,t0);
        s0 = 0;
        if ( min0 < min )
        {
            min = min0;
            r = r0;
            s = s0;
            t = t0;
        }

        triseg.b = triseg.b + tri.e0;
        triseg.m = triseg.m - tri.e0;
        min0 = MinLineSegmentLineSegment(seg,triseg,r0,t0);
        s0 = 1-t0;
        if ( min0 < min )
        {
            min = min0;
            r = r0;
            s = s0;
            t = t0;
        }

        min0 = MinPointTriangle(seg.b,tri,s0,t0);
        r0 = 0;
        if ( min0 < min )
        {
            min = min0;
            r = r0;
            s = s0;
            t = t0;
        }

        pt = seg.b+seg.m;
        min0 = MinPointTriangle(pt,tri,s0,t0);
        r0 = 1;
        if ( min0 < min )
        {
            min = min0;
            r = r0;
            s = s0;
            t = t0;
        }
    }

    return min;
}
//---------------------------------------------------------------------------

#ifdef LIN3TRI3_TEST

//#define RAND (2.0f*rand()/float(RAND_MAX)-1)

//ofstream ostr("data.txt");

void main ()
{
    Triangle3 tri;
    Line3 line;
    Point3 p, q, diff;
    double r, s, t, r0, s0, t0, min0, min1, dist;

    double maxdiff = 0;

    for (int i = 0; i < 128; i++)
    {
        tri.b.x = RAND;
        tri.b.y = RAND;
        tri.b.z = RAND;
        tri.e0.x = RAND;
        tri.e0.y = RAND;
        tri.e0.z = RAND;
        tri.e1.x = RAND;
        tri.e1.y = RAND;
        tri.e1.z = RAND;

        line.b.x = RAND;
        line.b.y = RAND;
        line.b.z = RAND;
        if ( i % 2 )
        {
            // non-parallel line and triangle
            line.m.x = RAND;
            line.m.y = RAND;
            line.m.z = RAND;
        }
        else
        {
            // line is parallel to triangle
            double c0 = RAND;
            double c1 = RAND;
            line.m = c0*tri.e0+c1*tri.e1;
        }

        min0 = FLT_MAX;
        int max = 32;
        for (int z = 0; z <= max; z++)
        {
            r0 = z/double(max);
            p = line.b+r0*line.m;
            for (int y = 0; y <= max; y++)
            {
                s0 = y/double(max);
                for (int x = 0; x+y <= max; x++)
                {
                    t0 = x/double(max);
                    q = tri.b+s0*tri.e0+t0*tri.e1;
                    
                    diff = p-q;
                    dist = Length(diff);
                    if ( dist < min0 )
                    {
                        min0 = dist;
                        r = r0;
                        s = s0;
                        t = t0;
                    }
                }
            }
        }
        ostr << "i = " << i << endl;
        ostr << "sampled = " << r << ' ' << s << ' ' << t << ' '
             << min0 << endl;

        min1 = MinLineSegmentTriangle(line,tri,r,s,t);
        ostr << "analytic = " << r << ' ' << s << ' ' << t << ' '
             << min1 << endl;

        ostr << "diff = " << min1-min0 << endl;

        if ( min1-min0 > maxdiff )
            maxdiff = min1-min0;

        ostr << endl;
    }

    ostr << "max diff = " << maxdiff << endl;
}
#endif

//FILE: triasect.cpp
//---------------------------------------------------------------------------
AgtLine 
AnalyticGeometryTool::EdgeToLine (Point2* v0, Point2* v1)
{
    // assert:  v0 and v1 are distinct

    Point2 edge = { v1->x - v0->x, v1->y - v0->y };

    AgtLine line;
    line.N.x = edge.y;
    line.N.y = -edge.x;
    line.c = line.N.x*v0->x + line.N.y*v0->y;

    return line;
}
//---------------------------------------------------------------------------
void 
AnalyticGeometryTool::TriangleLines (Triangle* T, AgtLine line[3])
{
    line[0] = EdgeToLine(&T->v[0],&T->v[1]);
    line[1] = EdgeToLine(&T->v[0],&T->v[2]);
    line[2] = EdgeToLine(&T->v[1],&T->v[2]);

    // make sure normals point outwards
    Point2 avr = { 0.0, 0.0 };
    int i;
    for (i = 0; i < 3; i++)
    {
        avr.x += T->v[i].x;
        avr.y += T->v[i].y;
    }
    static double oneThird = 1.0/3.0;
    avr.x *= oneThird;
    avr.y *= oneThird;

    for (i = 0; i < 3; i++)
    {
        double dot = line[i].N.x*avr.x+line[i].N.y*avr.y-line[i].c;
        if ( dot > 0.0 )
        {
            line[i].N.x = -line[i].N.x;
            line[i].N.y = -line[i].N.y;
            line[i].c = -line[i].c;
        }
    }
}
//---------------------------------------------------------------------------
AgtTriList* 
AnalyticGeometryTool::SplitAndDecompose (Triangle* T, AgtLine* line)
{
    double c[3];
    int positive = 0, ip[3];
    int negative = 0, in[3];
    int i;
    for (i = 0; i < 3; i++)
    {
        c[i] = line->N.x*T->v[i].x + line->N.y*T->v[i].y - line->c;
        if ( c[i] > 0.0 )
        {
            ip[positive++] = i;
        }
        else if ( c[i] < 0.0 )
        {
            in[negative++] = i;
        }
    }

    // For a split to occur, one of the c_i must be positive and one must
    // be negative.
    AgtTriList* list = NULL;

    if ( negative == 0 ) // T is completely on the positive side of line
        return 0;

    if ( positive == 0 )
    {
        // T is completely on the negative side of line
        list = new AgtTriList;
        list->tri = T;
        list->next = 0;
        return list;
    }

    // T is split by line.  Determine how it is split and how to decompose
    // the negative-side portion into triangles (3 cases).
    double w0, w1, cdiff;
    Point2 intrp[2];

    if ( positive == 2 )
    {
        // ++-
        for (i = 0; i < 2 /* = positive */; i++)
        {
            cdiff = c[ip[i]]-c[in[0]];
            w0 = -c[in[0]]/cdiff;
            w1 = c[ip[i]]/cdiff;
            T->v[ip[i]].x = w0*T->v[ip[i]].x+w1*T->v[in[0]].x;
            T->v[ip[i]].y = w0*T->v[ip[i]].y+w1*T->v[in[0]].y;
        }
        list = new AgtTriList;
        list->tri = T;
        list->next = 0;
    }
    else if ( positive == 1 )
    {
        if ( negative == 2 )
        {
            // +--
            for (i = 0; i < 2 /* = negative */; i++)
            {
                cdiff = c[ip[0]]-c[in[i]];
                w0 = -c[in[i]]/cdiff;
                w1 = c[ip[0]]/cdiff;
                intrp[i].x = w0*T->v[ip[0]].x+w1*T->v[in[i]].x;
                intrp[i].y = w0*T->v[ip[0]].y+w1*T->v[in[i]].y;
            }

            T->v[ip[0]] = intrp[0];
            list = new AgtTriList;
            list->tri = T;

            Triangle* T1 = new Triangle;
            T1->v[0] = intrp[0];
            T1->v[1] = T->v[in[1]];
            T1->v[2] = intrp[1];
            list->next = new AgtTriList;
            list->next->tri = T1;
            list->next->next = 0;
        }
        else
        {
            // +-0
            cdiff = c[ip[0]]-c[in[0]];
            w0 = -c[in[0]]/cdiff;
            w1 = c[ip[0]]/cdiff;
            T->v[ip[0]].x = w0*T->v[ip[0]].x+w1*T->v[in[0]].x;
            T->v[ip[0]].y = w0*T->v[ip[0]].y+w1*T->v[in[0]].y;
            list = new AgtTriList;
            list->tri = T;
            list->next = 0;
        }
    }

    return list;
}
//---------------------------------------------------------------------------
AgtTriList* 
AnalyticGeometryTool::Intersection (Triangle* T0, Triangle* T1)
{
    // build edges of T0
    AgtLine line[3];
    TriangleLines(T0,line);

    // initial list is copy of T1 (since triangle may be deleted)
    AgtTriList* list = new AgtTriList;
    list->tri = new Triangle;
    memcpy(list->tri,T1,sizeof(Triangle));
    list->next = 0;

    // process subtriangles of T1 against lines of T0
    for (int i = 0; i < 3; i++)
    {
        AgtTriList* save = 0;
        while ( list )
        {
            // get head of list
            Triangle* T = list->tri;
            
            AgtTriList* sad = SplitAndDecompose(T,&line[i]);
            if ( sad )
            {
                // search for end of list
                AgtTriList* end = sad;
                while ( end->next )
                    end = end->next;

                // attach decomposition to front of save-list
                end->next = save;
                save = sad;
            }

            // remove head of list
            AgtTriList* tmp = list;
            list = list->next;

            if( sad == NULL )
              delete tmp->tri;
            delete tmp;
        }
        list = save;
    }

    return list;
}
//---------------------------------------------------------------------------
double 
AnalyticGeometryTool::AreaTriangle (Triangle* T)
{
    Point2 e1 = { T->v[1].x - T->v[0].x, T->v[1].y - T->v[0].y };
    Point2 e2 = { T->v[2].x - T->v[0].x, T->v[2].y - T->v[0].y };

    return fabs(0.5*(e1.x*e2.y-e1.y*e2.x));
}
//---------------------------------------------------------------------------
double 
AnalyticGeometryTool::Area (AgtTriList* list)
{
    double area = 0.0;

    while ( list )
    {
        area += AreaTriangle(list->tri);
        list = list->next;
    }

    return area;
}
//---------------------------------------------------------------------------

double
AnalyticGeometryTool::AreaIntersection( Triangle &tri1, Triangle &tri2 )
{
  AgtTriList* list = Intersection(&tri1,&tri2);
  double area = Area(list);
  while ( list )
  {
    AgtTriList* save = list;
    list = list->next;
    delete save->tri;
    delete save;
  }
  delete list;
  return area;
}

#ifdef TRIASECT_TEST

void main ()
{
    Triangle T0;
    T0.v[0].x = 0.0;  T0.v[0].y = 0.0;
    T0.v[1].x = 1.0;  T0.v[1].y = 0.0;
    T0.v[2].x = 0.0;  T0.v[2].y = 1.0;

    Triangle T1;
    const double eps = 0.001;
    T1.v[0].x = 0.5+eps;  T1.v[0].y = 0.5+eps;
    T1.v[1].x = -eps;  T1.v[1].y = 0.5;
    T1.v[2].x = 0.5;  T1.v[2].y = -eps;

    AgtTriList* list = Intersection(&T0,&T1);
    double area = Area(list);

    while ( list )
    {
        AgtTriList* save = list;
        list = list->next;
        delete save->tri;
        delete save;
    }
}
#endif

double
AnalyticGeometryTool::Angle( Triangle3& tri1, Triangle3& tri2 )
{
  double norm1[3];
  Normal( tri1, norm1 );

  double norm2[3];
  Normal( tri2, norm2 );

  return angle_vec_vec( norm1, norm2 );
}


void
AnalyticGeometryTool::Normal( Triangle3& tri, double normal[3] )
{
  double vec1[3];
  vec1[0]=tri.e0.x; vec1[1]=tri.e0.y; vec1[2]=tri.e0.z;

  double vec2[3];
  vec2[0]=tri.e1.x; vec2[1]=tri.e1.y; vec2[2]=tri.e1.z;

  cross_vec_unit( vec1, vec2, normal );
}

double
AnalyticGeometryTool::ProjectedOverlap( Triangle3& tri1, Triangle3& tri2 )
{
  // Transform both triangles into local coordinate system of tri1
  // to eliminate z-coordinate. 

  // Use normal as z-axis
  double z[3];
  Normal( tri1, z );

  // x-axis goes from b to e0
  double x[3];
  x[0] = tri1.e0.x;
  x[1] = tri1.e0.y;
  x[2] = tri1.e0.z;

  // Cross z with x to get y-axis
  double y[3];
  cross_vec( z, x, y );

  // Unitize them all
  unit_vec( x, x );
  unit_vec( y, y );
  unit_vec( z, z );

  double origin[3];
  origin[0] = tri1.b.x;
  origin[1] = tri1.b.y;
  origin[2] = tri1.b.z;

  double mtxGlobal2TriOne[4][4];
  vecs_to_mtx( x, y, z, origin, mtxGlobal2TriOne ); // Really mtxTriOne2Global
  inv_trans_mtx( mtxGlobal2TriOne, mtxGlobal2TriOne );
  
  double v0[3], v1[3], v2[3]; // Vertices of triangle
  v0[0]=tri1.b.x;   v0[1]=tri1.b.y;  v0[2]=tri1.b.z;
  v1[0]=tri1.e0.x+tri1.b.x; v1[1]=tri1.e0.y+tri1.b.y; v1[2]=tri1.e0.z+tri1.b.z;
  v2[0]=tri1.e1.x+tri1.b.x; v2[1]=tri1.e1.y+tri1.b.y; v2[2]=tri1.e1.z+tri1.b.z;
  transform_pnt( mtxGlobal2TriOne, v0, v0 );
  transform_pnt( mtxGlobal2TriOne, v1, v1 );
  transform_pnt( mtxGlobal2TriOne, v2, v2 );

  // Load this into a 2D triangle T1 ( z-axis is now zero )
  Triangle T1;
  T1.v[0].x = v0[0]; T1.v[0].y = v0[1];
  T1.v[1].x = v1[0]; T1.v[1].y = v1[1]; 
  T1.v[2].x = v2[0]; T1.v[2].y = v2[1];

  // Setup the next triangle, in this coordinate system
  v0[0]=tri2.b.x;   v0[1]=tri2.b.y;  v0[2]=tri2.b.z;
  v1[0]=tri2.e0.x+tri2.b.x; v1[1]=tri2.e0.y+tri2.b.y; v1[2]=tri2.e0.z+tri2.b.z;
  v2[0]=tri2.e1.x+tri2.b.x; v2[1]=tri2.e1.y+tri2.b.y; v2[2]=tri2.e1.z+tri2.b.z;

  // The following creates coordinates of tri2 in tri1's csys
  transform_pnt( mtxGlobal2TriOne, v0, v0 );
  transform_pnt( mtxGlobal2TriOne, v1, v1 );
  transform_pnt( mtxGlobal2TriOne, v2, v2 );

  // Now we can project tri2 to tri1 simply by dropping the z-coordinate
  Triangle T2;
  T2.v[0].x = v0[0]; T2.v[0].y = v0[1];
  T2.v[1].x = v1[0]; T2.v[1].y = v1[1]; 
  T2.v[2].x = v2[0]; T2.v[2].y = v2[1];

  // Now find area of overlap
  return AreaIntersection( T1, T2 );
}
