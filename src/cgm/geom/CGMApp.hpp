//-------------------------------------------------------------------------
// Filename      : CGMApp.hpp 
//
// Purpose       : This file represents the Cubit application itself.
//
// Special Notes : 
//
// Creator       : Byron Hanks
//
// Date          : 06/08/98
//
// Owner         : Darryl Melander
//-------------------------------------------------------------------------

//- *********************************************************************
//- Copyright notice:
//-
//- This program was prepared by Sandia National Laboratories under
//- Contract No. DE-AC04-94-AL8500 with the U.S. Department of Energy
//- (DOE).  All rights in the program are reserved by DOE on behalf of the
//- Government and Sandia pursuant to the contract.  You are authorized to
//- use this program for U.S. Government purposes, but it is not to be
//- released or distributed to the public.  Neither the U.S. Government
//- nor Sandia makes any warranty, express or implied, or assumes any
//- liability or responsibility for the use of this software.
//-
//- *********************************************************************
//- More official Copyright text
//-
//- Copyright 2001 Sandia Corporation.  Under the terms of Contract
//- DEAC04-94AL85000, there is a non-exclusive license for use of this
//- work by or on behalf of the U.S. Government.  Export of this program
//- may require a license from the United States Government.
//-
//- *********************************************************************

#ifndef CGM_APP_HPP
#define CGM_APP_HPP

#include "CubitDefines.h"
#include "CubitAttribManager.hpp"
#include "CubitGeomConfigure.h"

class CUBIT_GEOM_EXPORT CGMApp
{
public:
   static CGMApp* instance();
     //- Access to the application object

   static void delete_instance() {if (instance_) {delete instance_; instance_ = NULL;}};
        
   ~CGMApp();

  void startup(int argc, char **argv);
   //-  Contains startup code for cubit

  void shutdown();
   //- Contains shutdown code for cubit

  static void initialize_settings();

  CubitAttribManager* attrib_manager();

  void save_current_attribute_states();
  void restore_previous_attribute_states();
  
private:

   static CGMApp* instance_;
   CubitBoolean mAppStarted;
   CubitAttribManager mAttribManager;

   CGMApp();

   void register_attributes();
   
   CubitBoolean *caUpdateFlgs;
   CubitBoolean *caActuateFlgs;
   CubitBoolean *caWriteFlgs;
   CubitBoolean *caReadFlgs;
};


#endif

