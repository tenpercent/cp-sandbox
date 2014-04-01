# Ani3D/src/Tutorials/PackageLMR

add_executable       (Nodal2MetricRecovery mainNodal2Metric.f)
target_link_libraries(Nodal2MetricRecovery aniLMR3D aniMBA3D aniFEM3D)

add_executable       (Est2MetricRecovery mainEst2Metric.f esterr.f)
target_link_libraries(Est2MetricRecovery aniLMR3D aniMBA3D aniFEM3D)

add_executable       (Func2GradMetric mainFunc2GradMetric.f)
target_link_libraries(Func2GradMetric aniLMR3D aniMBA3D aniFEM3D)

if(WIN32)
   set_target_properties(Nodal2MetricRecovery PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(Est2MetricRecovery PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
   set_target_properties(Func2GradMetric PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
endif(WIN32)
