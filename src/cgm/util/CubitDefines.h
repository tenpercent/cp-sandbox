/* Class:       CubitDefines 
 * Description: CubitDefines class - all global definitions for CUBIT.
 * Owner:       Tim Tautges
 * Checked by:
 * Version: $Id: 
 */

#ifndef CUBITOBJECT_HPP
#define CUBITOBJECT_HPP

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <sys/types.h>
#include <assert.h>
#include "CubitUtilConfigure.h"

#ifdef NT
#pragma warning ( 4 : 4291 4244 4305 4018 4786)
#endif
/* Adds DBL_MIN, DBL_MAX definitions on some platforms*/
#include <float.h> 

/* sanity check defines (these ARE needed!) */

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef NULL
#define NULL 0
#endif

/* typedef for integers and flags of just a few bits.
 * Usage: "Bit isMarked : 1;" where 1 is bit width of isMarked.
 * If there are just a few small flags in a class, cBit and IttyBit might
 * save some space; the length of int and short are machine dependent,
 * but char is as short as possible on all platforms.
 */
typedef unsigned int Bit;
typedef unsigned short IttyBit;
typedef unsigned char cBit;

static int const INVALID_ENTITY_ID = -1 ;

/* Boolean flags. */
#ifdef __cplusplus
  typedef bool CubitBoolean;
  const bool CUBIT_TRUE = true;
  const bool CUBIT_FALSE = false;
#else
  typedef char CubitBoolean;
  #define CUBIT_TRUE 1;
  #define CUBIT_FALSE 0;
#endif

/* Completion status */
enum CubitStatus { CUBIT_FAILURE = 0, CUBIT_SUCCESS = 1 } ;

/* Flags that indicate the Forward and Reversed senses. */
enum CubitSense {CUBIT_UNKNOWN = -1, CUBIT_FORWARD = 0, CUBIT_REVERSED = 1};

/* Flags to indicate whether a point is outside or off of an entity, inside 
   or on an entity, or on the boundary of an entity.  Defined so that calling 
   code can do true/false checks if it doesn't care if point is inside or on 
   the boundary. */
enum CubitPointContainment {CUBIT_PNT_UNKNOWN = -1, CUBIT_PNT_OUTSIDE = 0, 
                            CUBIT_PNT_OFF = 0, CUBIT_PNT_INSIDE = 1, 
                            CUBIT_PNT_ON = 1, CUBIT_PNT_BOUNDARY = 2};

enum FirmnessType {LIMP, SOFT, FIRM, HARD};
/* Firmness continuum - how bad does the user want something? 
   LIMP means "don't care" or "unset", SOFT means "about right" or
   "automatically set", HARD means "user set" or "must be this way".
*/

/* Relationships between types of objects in the Entity Relationship
 * Diagram */
enum CubitEntityRelation 
{ 
    CUBIT_RELATION_NONE = -1, 
    CUBIT_RELATION_PARENT = 0, 
    CUBIT_RELATION_CHILD = 1 
} ;

/* Flag types for associativity data stored out and read back from
 * ExodusII files */
enum CubitAssocDataType 
{ 
   CUBIT_ASSOC_NULL = 0, 
   CUBIT_ASSOC_GEOM_REFENTITY = 1, 
   CUBIT_ASSOC_BC_BLOCK = 2, 
   CUBIT_ASSOC_BC_NSET = 3, 
   CUBIT_ASSOC_BC_SSET = 4 
} ;

enum ObserverType {
    MODEL_OBSERVER,
      /* big brother observer */

    DRAWING_OBSERVER,
      /* any observer concerned with drawing */

    GUI_OBSERVER
      /* a gui-type observer */
};

/* CSG operations for facetbool */
enum CubitFacetboolOp
{
   CUBIT_FB_UNION,
   CUBIT_FB_INTERSECTION,
   CUBIT_FB_SUBTRACTION
};

/* #define's for the entire CUBIT system */

#ifdef __cplusplus
#ifdef M_PI
const double CUBIT_PI           =              M_PI;
#else
const double CUBIT_PI           =              3.1415926535897932384626;
#endif
#else
#ifdef M_PI
#define CUBIT_PI         M_PI
#else
#define CUBIT_PI         3.1415926535897932384626
#endif
#endif
#define DEGREES_TO_RADIANS(angle) ( (angle) / 180.0 * CUBIT_PI )
#define RADIANS_TO_DEGREES(angle) ( 180.0 * (angle) / CUBIT_PI )
#ifdef __cplusplus
#ifdef INT_MAX
const int CUBIT_INT_MAX = INT_MAX;
#else
const int CUBIT_INT_MAX = 2147483647;
const int INT_MAX = CUBIT_INT_MAX;
#endif
#else
#ifdef INT_MAX
#define CUBIT_INT_MAX INT_MAX
#else
#define CUBIT_INT_MAX 2147483647
#define INT_MAX CUBIT_INT_MAX
#endif
#endif

#ifdef __cplusplus
#ifdef INT_MIN
const int CUBIT_INT_MIN = INT_MIN;
#else
const int CUBIT_INT_MIN = -2147483647;
const int INT_MIN = CUBIT_INT_MIN;
#endif
#else
#ifdef INT_MIN
#define CUBIT_INT_MIN INT_MIN
#else
#define CUBIT_INT_MIN -2147483647
#define INT_MIN CUBIT_INT_MIN
#endif
#endif

#ifdef __cplusplus
#ifdef DBL_MIN
const double CUBIT_DBL_MIN      =              DBL_MIN;
#else
const double CUBIT_DBL_MIN      =              1.0E-30;
#endif
#else
#ifdef DBL_MIN
#define CUBIT_DBL_MIN           DBL_MIN
#else
#define CUBIT_DBL_MIN           1.0E-30
#endif
#endif

#ifdef __cplusplus
#ifdef DBL_MAX
const double CUBIT_DBL_MAX      =              DBL_MAX;
#else
const double CUBIT_DBL_MAX      =              1.0E30;
#endif 
#else
#ifdef DBL_MAX
#define CUBIT_DBL_MAX           DBL_MAX
#else
#define CUBIT_DBL_MAX           1.0E30
#endif 
#endif

#define CUBIT_MIN(a,b)                   ( (a) < (b) ? (a) : (b) )
#define CUBIT_MAX(a,b)                   ( (a) > (b) ? (a) : (b) )

#define CUBIT_MIN_4(a,b,c,d)             (( (a) < (b) ? (a) : (b) ) < \
                                          ( (c) < (d) ? (c) : (d) ) ? \
                                          ( (a) < (b) ? (a) : (b) ) : \
                                          ( (c) < (d) ? (c) : (d) ))

#define CUBIT_MAX_4(a,b,c,d)             (( (a) > (b) ? (a) : (b) ) > \
                                          ( (c) > (d) ? (c) : (d) ) ? \
                                          ( (a) > (b) ? (a) : (b) ) : \
                                          ( (c) > (d) ? (c) : (d) ))


  
/*  Setting this flag (asynchronously) to CUBIT_TRUE 
 *  will cause Cubit/CGM to attempt to abort any current
 *  operations and return.  
 * 
 * NOTE: IT IS THE RESPONSIBILITY OF THE APPLICATION 
 *       USING CGM TO RESET THIS FLAG TO CUBIT_FALSE!!!
 *       For Cubit, this flag is reset in UserInterface.
 *       
 *  CubitApp provides a default signal hander that will
 *  set this flag to CUBIT_TRUE whenever a SIGINT (^C)
 *  is detected.  It is still the responsibility of the
 *  application using CGM to reset the flag to false!
 *  The signal hander provided by CubitApp may be set by 
 *  calling CubitApp::instance()->catch_interrupt(CUBIT_TRUE).
 *
 *  The storage space for this flag is defined in CubitApp.cpp
 *  and initialized in CubitApp::initialize().
 */    
CUBIT_UTIL_EXPORT extern volatile CubitBoolean cubit_intr;

/* Operators on enum types */
#if 0 /* Does not work on broken Microsoft compiler */
#ifdef __cplusplus

  /* Logical Operators for CubitBoolean type */
  
  inline CubitBoolean operator!( const CubitBoolean b )
  { return (CubitBoolean)(!(bool)b); }

  inline CubitBoolean operator&&( const CubitBoolean a, const CubitBoolean b )
  { return (CubitBoolean)((bool)a && (bool)b); }

  inline CubitBoolean operator||( const CubitBoolean a, const CubitBoolean b )
  { return (CubitBoolean)((bool)a || (bool)b); }


  /* Logical Operators for CubitStatus type */

  inline CubitStatus operator!( const CubitStatus b )
  { return (CubitStatus)(!(bool)b); }

  inline CubitStatus operator&&( const CubitStatus a, const CubitStatus b )
  { return (CubitStatus)((bool)a && (bool)b); }

  inline CubitStatus operator||( const CubitStatus a, const CubitStatus b )
  { return (CubitStatus)((bool)a || (bool)b); }


  /* Arithmatic Operators for CubitSense type */

  inline CubitSense operator-( const CubitSense s )
  { return (s == CUBIT_UNKNOWN) ? CUBIT_UNKNOWN : (CubitSense)(1 - s); }

  inline CubitSense operator*( const CubitSense a, const CubitSense b )
  { return (( a == CUBIT_UNKNOWN ) || ( b == CUBIT_UNKNOWN )) ? 
      CUBIT_UNKNOWN : ( a == b ) ? CUBIT_FORWARD : CUBIT_REVERSED; }

  inline CubitSense operator*=( CubitSense& a, const CubitSense b )
  { return a = a * b; }

#endif
#endif

#endif

