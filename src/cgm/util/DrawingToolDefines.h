// define constants that are part of the public DrawingTool interface
//   here in DrawingTool.hpp (versus DrawingToolInstance.hpp)
#ifndef DRAWINGTOOL_DEFINES_H
#define DRAWINGTOOL_DEFINES_H
const int CUBIT_INVALID_DRAW_MODE     = 0;
const int CUBIT_WIRE_FRAME       = 1;
const int CUBIT_HIDDEN_LINE      = 2;
const int CUBIT_POLYGON_FILL     = 3;
const int CUBIT_FLAT_SHADING     = 4;
const int CUBIT_SMOOTH_SHADING   = 5;
const int CUBIT_PAINTERS         = 6;
const int CUBIT_TRUE_HIDDEN_LINE = 7;
const int CUBIT_SMOOTH_TRANSPARENT = 8;

const int CUBIT_SCREEN    = 1;
const int CUBIT_WORLD     = 2;
const int CUBIT_CAMERA    = 3;

const int CUBIT_CAMERA_MOVE_FLY  = 1;
const int CUBIT_CAMERA_MOVE_PAN  = 2;
const int CUBIT_CAMERA_MOVE_SPIN = 3;


enum CubitGraphicsPointStyle { CUBIT_GRAPHICS_POINT_OPEN_CIRCLE,
                               CUBIT_GRAPHICS_POINT_SOLID_CIRCLE,
                               CUBIT_GRAPHICS_POINT_CIRCLED_DOT,
                               CUBIT_GRAPHICS_POINT_X,
                               CUBIT_GRAPHICS_POINT_PLUS         };


// define constants for the labeling type
enum CubitLabelType 
{
   CUBIT_LABEL_NONE=0,
   CUBIT_LABEL_ID,
   CUBIT_LABEL_GENESIS_ID,
   CUBIT_LABEL_NAME,
   CUBIT_LABEL_INTERVAL,
   CUBIT_LABEL_SIZE,
   CUBIT_LABEL_MERGE,
   CUBIT_LABEL_IS_MERGED,
   CUBIT_LABEL_FIRMNESS,
   CUBIT_LABEL_SCHEME,
   CUBIT_LABEL_NAME_ID,
   CUBIT_LABEL_NAME_ONLY
};



// OR these together as the parameter to set_mouse_zoom_direction()
#define CUBIT_ZOOM_UP    1
#define CUBIT_ZOOM_DOWN  2
#define CUBIT_ZOOM_RIGHT 4
#define CUBIT_ZOOM_LEFT  8

#endif

