# Ani3D/src/Tutorials/PackageVIEW

add_executable       (mainView main.f)
target_link_libraries(mainView aniVIEW3D aniMBA3D)

if(WIN32)
   set_target_properties(mainView PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
endif(WIN32)
