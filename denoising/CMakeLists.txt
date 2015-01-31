set(BUILD_SHARED_LIBS false)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${BoostSourceDir})

############
# denoising
############
set(DENOISE_SRCS
    ${DENOISE_SRCS}
    ${COMMON_SRCS}
    denoising_main.cpp
    helpers.h
    condat_denoise.h)

add_executable(denoising ${DENOISE_SRCS})
add_dependencies(denoising boost_program_options 
						   boost_log 
						   boost_system 
						   boost_thread
						   boost_filesystem
						   boost_date_time)
target_link_libraries(denoising 
		      ${BoostProgramOptionsLibs} 
		      ${BoostLogLibs} 
		      ${BoostSystemLibs} 
		      ${BoostThreadLibs}
			  ${BoostFilesystemLibs}
			  ${BoostDateTimeLibs}
			  ${BoostChronoLibs})
if(UNIX)
    target_link_libraries(denoising rt)
endif()
			  
default_target_compile_options(denoising)

############
# lambdaopt
############
set(LAMBDAOPT_SRCS
    ${LAMBDAOPT_SRCS}
    ${COMMON_SRCS}
    lambdaopt_main.cpp
    helpers.h
    lambda_max.h
    lambda_opt.h)

add_executable(lambdaopt ${LAMBDAOPT_SRCS})
add_dependencies(lambdaopt boost_program_options 
						   boost_log 
						   boost_system 
						   boost_thread
						   boost_filesystem
						   boost_date_time)
target_link_libraries(lambdaopt
		      ${BoostProgramOptionsLibs} 
		      ${BoostLogLibs} 
		      ${BoostSystemLibs} 
		      ${BoostThreadLibs}
			  ${BoostFilesystemLibs}
			  ${BoostDateTimeLibs}
			  ${BoostChronoLibs})

if(UNIX)
    target_link_libraries(lambdaopt rt)
endif()
			
			  
default_target_compile_options(lambdaopt)
