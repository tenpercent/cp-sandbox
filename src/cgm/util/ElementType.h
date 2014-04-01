#ifndef ELEMENT_TYPE_H
#define ELEMENT_TYPE_H
/*This list and the one in ElementBlock.cc must be
  syncronized correctly.  The invalid type MUST be 0
  if we have an unknown element type in the .cc file.
  DRW 10/20/97                                       */
  
enum ElementType { SPHERE_EXO=0,
                   BAR, BAR2, BAR3,
                   BEAM, BEAM2, BEAM3,
                   TRUSS, TRUSS2, TRUSS3,
                   SPRING,
                   TRI, TRI3, TRI6, TRI7,
                   TRISHELL, TRISHELL3, TRISHELL6, TRISHELL7,
                   SHEL, SHELL4, SHELL8, SHELL9,
                   QUAD, QUAD4, QUAD5, QUAD8, QUAD9,
                   TETRA, TETRA4, TETRA8, TETRA10, TETRA14,
                   PYRAMID, PYRAMID5, PYRAMID8, PYRAMID13, PYRAMID18,
                   HEX, HEX8, HEX9, HEX20, HEX27, HEXSHELL, 
                   INVALID_ELEMENT_TYPE};
/*Invalid element type must be the last type...*/
#endif

