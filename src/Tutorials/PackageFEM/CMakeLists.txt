# Ani3D/src/Tutorials/PackageFEM

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../aniFEM)

add_executable       (femTetrahedron mainTetrahedron.f)
target_link_libraries(femTetrahedron aniFEM3D aniMBA3D)

add_executable       (femSimple      mainSimple.f)
target_link_libraries(femSimple      aniFEM3D aniMBA3D)

add_executable       (femBC          mainBC.f)
target_link_libraries(femBC          aniFEM3D aniMBA3D)

add_executable       (femTemplate    mainTemplate.f)
target_link_libraries(femTemplate    aniFEM3D aniMBA3D)

if(WIN32)
   set_target_properties(femTetrahedron PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(femSimple      PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(femBC          PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(femTemplate    PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
endif(WIN32)
