//- Title: Database Storage Size Constants
//- Filename: database.hpp
//- Description: This file sets initial storage sizes for database entities
//-             storage is dynamic and will not overflow unless free store is
//-             exhausted.  These values affect efficiency, not robustness
//- Owner: Greg Sjaardema
//- Checked By: Mark Whitely, June 13, 1994
//- Version: $Id: 

const int HEX_INCREMENT      =     5000;
const int NODE_INCREMENT     =     HEX_INCREMENT;
const int FACE_USE_INCREMENT = 6 * HEX_INCREMENT;
const int FACE_INCREMENT     = 4 * HEX_INCREMENT;
const int EDGE_USE_INCREMENT = 4 * FACE_INCREMENT;
const int EDGE_INCREMENT     =     FACE_INCREMENT;
const int TRI_INCREMENT      =     FACE_INCREMENT;
const int TET_INCREMENT      =     HEX_INCREMENT;

const int MESH_SIZE_PER_EDGE    = 15;
const int MESH_SIZE_PER_SURFACE = MESH_SIZE_PER_EDGE    * MESH_SIZE_PER_EDGE;
const int MESH_SIZE_PER_VOLUME    = MESH_SIZE_PER_SURFACE * MESH_SIZE_PER_EDGE;
