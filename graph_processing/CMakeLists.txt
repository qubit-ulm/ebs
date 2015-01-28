set(BUILD_SHARED_LIBS false)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${BoostSourceDir})

###################
# graph_processing
###################
set(GP_SRCS
    ${GP_SRCS}
    ${COMMON_SRCS}
    graph_processing_main.cpp
    helpers.h
    cmd_helpers.h
    bk_max_flow.h
    energy.h
    sitesstore.h
    binopt.h
    counting_statistics.h
    runtime_statistics.h)

add_executable(graph_processing ${GP_SRCS})
add_dependencies(graph_processing boost_program_options 
								  boost_log 
								  boost_system 
								  boost_thread
								  boost_filesystem
								  boost_date_time
								  boost_graph)
target_link_libraries(graph_processing
                      ${BoostProgramOptionsLibs} 
                      ${BoostLogLibs} 
                      ${BoostSystemLibs} 
                      ${BoostThreadLibs}
                      ${BoostFilesystemLibs}
                      ${BoostDateTimeLibs}
                      ${BoostChronoLibs}
                      ${BoostGraphLibs})

default_target_compile_options(graph_processing)

##################
# level_generator
##################
set(LG_SRCS
    ${LG_SRCS}
    ${COMMON_SRCS}
    level_generator_main.cpp
    helpers.h
    cmd_helpers.h)

add_executable(level_generator ${LG_SRCS})
add_dependencies(level_generator boost_program_options
                                 boost_log
                                 boost_system
                                 boost_thread
                                 boost_filesystem
                                 boost_date_time)
target_link_libraries(level_generator
                      ${BoostProgramOptionsLibs} 
                      ${BoostLogLibs} 
                      ${BoostSystemLibs} 
                      ${BoostThreadLibs}
                      ${BoostFilesystemLibs}
                      ${BoostDateTimeLibs}
                      ${BoostChronoLibs})

default_target_compile_options(level_generator)


