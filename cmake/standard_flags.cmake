
if(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LibCXX} ${LibCXXAbi} -lpthread")
  if(CcacheExe)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # Workaround for GCC bug https://bugs.launchpad.net/ubuntu/+source/gcc-defaults/+bug/1228201
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed -pthread")
endif()
