//-------------------------------------------------------------------------
// Filename      : CGMApp.cpp
//
// Purpose       : Init and cleanup for CGM
//
// Special Notes : 
//
// Creator       : Byron Hanks
//
// Date          : 01/09/2002
//
// Owner         : 
//-------------------------------------------------------------------------

#include "CGMApp.hpp"
#include "AppUtil.hpp"

#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "GeometryHealerTool.hpp"
#include "BoundingBoxTool.hpp"
#include "RefEntityName.hpp"
#include "SurfaceOverlapTool.hpp"
#include "CubitSimpleAttrib.hpp"
#include "TDUniqueId.hpp"
#include "MergeTool.hpp"
#include "OldUnmergeCode.hpp"

// TODO - I'm not sure this is the place to include all the attributes
#include "CADefines.hpp"
#include "CAMergePartner.hpp"
#include "CAEntityName.hpp"
#include "CAGroup.hpp"
#include "CAEntityId.hpp"
#include "CAUniqueId.hpp"
#include "CADeferredAttrib.hpp"
#include "CAEntityColor.hpp"
#include "CAMergeStatus.hpp"
#include "CASourceFeature.hpp"

CGMApp* CGMApp::instance_ = NULL;

// Returns the singleton instance of the app object.
CGMApp* CGMApp::instance()
{
   if (instance_ == NULL)
   {
      instance_ = new CGMApp();
      if (!instance_ )
      {
         PRINT_ERROR(" *** Unable to initialize application ***\n");
         exit(1);
      }
   }
   return instance_;
}

CGMApp::CGMApp()
{
    mAppStarted = CUBIT_FALSE;

      // CGMApp depends on there being an AppUtil
    AppUtil::instance();

      // initialize my settings
    initialize_settings();    
  
    caUpdateFlgs = NULL;
    caActuateFlgs = NULL;
    caWriteFlgs = NULL;
    caReadFlgs = NULL;
}

CGMApp::~CGMApp()
{
    instance_ = NULL;
}

void CGMApp::startup(int argc, char **argv)
{
   if (mAppStarted)
       return;

    // make sure apputil has started
  AppUtil::instance()->startup(argc, argv);

     // register attributes
   register_attributes();

   mAppStarted = CUBIT_TRUE;
}

void CGMApp::shutdown()
{
   delete GeometryHealerTool::instance();
   delete GeometryModifyTool::instance();
   delete GeometryQueryTool::instance();

   CGMApp::delete_instance();

   AppUtil::instance()->shutdown();
   AppUtil::delete_instance();
}

void CGMApp::initialize_settings()
{
  GeometryQueryTool::instance()->initialize_settings();
  GeometryModifyTool::instance()->initialize_settings();
  BoundingBoxTool::initialize_settings();  
  CubitSimpleAttrib::initialize_settings();
  RefEntityName::initialize_settings();
  SurfaceOverlapTool::initialize_settings();
  MergeTool::initialize_settings();
  OldUnmergeCode::initialize_settings();
  
}

void CGMApp::register_attributes()
{
  CubitStatus result;
  result = mAttribManager.register_attrib_type(CA_MERGE_PARTNER, "merge", "MERGE_PARTNER", 
                                               CAMergePartner_creator, CUBIT_FALSE,
                                               CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_FALSE, CUBIT_FALSE);
  assert (CUBIT_SUCCESS == result);

  result = mAttribManager.register_attrib_type(CA_ENTITY_NAME, "name", "ENTITY_NAME", 
                                               CAEntityName_creator, CUBIT_TRUE,
                                               CUBIT_TRUE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_TRUE, CUBIT_FALSE);
  assert (CUBIT_SUCCESS == result);

  result = mAttribManager.register_attrib_type(CA_GROUP, "group", "GROUP",
                                               CAGroup_creator, CUBIT_FALSE,
                                               CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_FALSE, CUBIT_FALSE);
  assert (CUBIT_SUCCESS == result);

  result = mAttribManager.register_attrib_type(CA_ENTITY_ID, "id", "ENTITY_ID",
                                               CAEntityId_creator, CUBIT_FALSE,
                                               CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_TRUE, CUBIT_FALSE);
  assert (CUBIT_SUCCESS == result);

  result = mAttribManager.register_attrib_type(CA_UNIQUE_ID, "unique id", "UNIQUE_ID",
                                               CAUniqueId_creator, CUBIT_FALSE,
                                               CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_TRUE, CUBIT_FALSE);
  assert (CUBIT_SUCCESS == result);

  result = mAttribManager.register_attrib_type(CA_DEFERRED_ATTRIB, "deferred attrib", "DEFERRED_ATTRIB",
                                               CADeferredAttrib_creator, CUBIT_FALSE,
                                               CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_FALSE, CUBIT_FALSE);
  assert (CUBIT_SUCCESS == result);

  result = mAttribManager.register_attrib_type(CA_ENTITY_COLOR, "color", "ENTITY_COLOR",
                                               CAEntityColor_creator, CUBIT_FALSE,
                                               CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_FALSE, CUBIT_TRUE);
  assert (CUBIT_SUCCESS == result);

  result = mAttribManager.register_attrib_type(CA_MERGE_STATUS, "merge status", "MERGE_STATUS",
                                               CAMergeStatus_creator, CUBIT_FALSE,
                                               CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_FALSE, CUBIT_FALSE);
  assert (CUBIT_SUCCESS == result);

  result = CGMApp::instance()->attrib_manager()->register_attrib_type(
                                               CA_SOURCE_FEATURE, "source feature", "SOURCE_FEATURE",
                                               CASourceFeature_creator, CUBIT_TRUE,
                                               CUBIT_TRUE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_TRUE, CUBIT_FALSE);
  assert (CUBIT_SUCCESS == result);
}

CubitAttribManager* CGMApp::attrib_manager()
{
  return &mAttribManager;
}


void CGMApp::save_current_attribute_states()
{
  //all pointers must be NULL
  if( !!caUpdateFlgs && !!caActuateFlgs &&
      !!caWriteFlgs && !!caReadFlgs )
  {
    PRINT_ERROR("Problem setting attribute flags\n");
    return;
  }
    
  DLIList<int> attrib_types;
  mAttribManager.get_registered_types(attrib_types);

  caUpdateFlgs = new CubitBoolean[attrib_types.size()];
  caActuateFlgs = new CubitBoolean[attrib_types.size()];
  caWriteFlgs = new CubitBoolean[attrib_types.size()];
  caReadFlgs = new CubitBoolean[attrib_types.size()];

  int num_att = attrib_types.size();
  attrib_types.reset();
  int i;
  for( i=0; i<num_att; i++ )
  {
    int att_type = attrib_types.get();
    caUpdateFlgs[i] = mAttribManager.auto_update_flag( att_type );
    caActuateFlgs[i] = mAttribManager.auto_actuate_flag( att_type );
    caWriteFlgs[i] = mAttribManager.auto_write_flag( att_type );
    caReadFlgs[i] = mAttribManager.auto_read_flag( att_type );
    attrib_types.step();
  }
}

void CGMApp::restore_previous_attribute_states()
{
  DLIList<int> attrib_types;
  mAttribManager.get_registered_types(attrib_types);

  int num_att = attrib_types.size();
  attrib_types.reset();
  int i;
  for( i=0; i<num_att; i++ )
  {
    int att_type = attrib_types.get();
    mAttribManager.set_auto_update_flag( att_type, caUpdateFlgs[i] );
    mAttribManager.set_auto_actuate_flag( att_type, caActuateFlgs[i] );
    mAttribManager.set_auto_write_flag( att_type, caWriteFlgs[i] );
    mAttribManager.set_auto_read_flag( att_type, caReadFlgs[i] );
    attrib_types.step();
  }
  
  delete [] caUpdateFlgs;
  delete [] caActuateFlgs;
  delete [] caWriteFlgs;
  delete [] caReadFlgs;

  caUpdateFlgs = NULL;
  caActuateFlgs = NULL;
  caWriteFlgs = NULL;
  caReadFlgs = NULL;
}

