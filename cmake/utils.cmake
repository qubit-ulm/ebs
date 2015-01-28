function(underscores_to_camel_case VarIn VarOut)
  string(REPLACE "_" ";" Pieces ${VarIn})
  foreach(Part ${Pieces})
    string(SUBSTRING ${Part} 0 1 Initial)
    string(SUBSTRING ${Part} 1 -1 Part)
    string(TOUPPER ${Initial} Initial)
    set(CamelCase ${CamelCase}${Initial}${Part})
  endforeach()
  set(${VarOut} ${CamelCase} PARENT_SCOPE)
endfunction()

function(get_target_architecture)
  if(TargetArchitecture)
    return()
  endif()
  if(APPLE AND CMAKE_OSX_ARCHITECTURES)
    # On OS X we use CMAKE_OSX_ARCHITECTURES *if* it was set
    # First let's normalize the order of the values

    # Note that it's not possible to compile PowerPC applications if you are using
    # the OS X SDK version 10.6 or later - you'll need 10.4/10.5 for that, so we
    # disable it by default
    # See this page for more information:
    # http://stackoverflow.com/questions/5333490/how-can-we-restore-ppc-ppc64-as-well-as-full-10-4-10-5-sdk-support-to-xcode-4

    # Architecture defaults to i386 or ppc on OS X 10.5 and earlier, depending on the CPU type detected at runtime.
    # On OS X 10.6+ the default is x86_64 if the CPU supports it, i386 otherwise.
    foreach(osx_arch ${CMAKE_OSX_ARCHITECTURES})
      if("${osx_arch}" STREQUAL "ppc" AND ppc_support)
        set(osx_arch_ppc TRUE)
      elseif("${osx_arch}" STREQUAL "i386")
        set(osx_arch_i386 TRUE)
      elseif("${osx_arch}" STREQUAL "x86_64")
        set(osx_arch_x86_64 TRUE)
      elseif("${osx_arch}" STREQUAL "ppc64" AND ppc_support)
        set(osx_arch_ppc64 TRUE)
      else()
        message(FATAL_ERROR "Invalid OS X arch name: ${osx_arch}")
      endif()
    endforeach()

    # Now add all the architectures in our normalized order
    if(osx_arch_ppc)
      list(APPEND ARCH ppc)
    endif()

    if(osx_arch_i386)
      list(APPEND ARCH i386)
    endif()

    if(osx_arch_x86_64)
      list(APPEND ARCH x86_64)
    endif()

    if(osx_arch_ppc64)
      list(APPEND ARCH ppc64)
    endif()
  else()

    # Based on the Qt 5 processor detection code, so should be very accurate
    # https://qt.gitorious.org/qt/qtbase/blobs/master/src/corelib/global/qprocessordetection.h
    # Currently handles arm (v5, v6, v7), x86 (32/64), ia64, and ppc (32/64)
    #
    # Regarding POWER/PowerPC, just as is noted in the Qt source,
    # "There are many more known variants/revisions that we do not handle/detect."
    set(archdetect_c_code "
    #if defined(__arm__) || defined(__TARGET_ARCH_ARM)
    #if defined(__ARM_ARCH_7__) \\
    || defined(__ARM_ARCH_7A__) \\
    || defined(__ARM_ARCH_7R__) \\
    || defined(__ARM_ARCH_7M__) \\
    || (defined(__TARGET_ARCH_ARM) && __TARGET_ARCH_ARM-0 >= 7)
    #error cmake_ARCH armv7
    #elif defined(__ARM_ARCH_6__) \\
    || defined(__ARM_ARCH_6J__) \\
    || defined(__ARM_ARCH_6T2__) \\
    || defined(__ARM_ARCH_6Z__) \\
    || defined(__ARM_ARCH_6K__) \\
    || defined(__ARM_ARCH_6ZK__) \\
    || defined(__ARM_ARCH_6M__) \\
    || (defined(__TARGET_ARCH_ARM) && __TARGET_ARCH_ARM-0 >= 6)
    #error cmake_ARCH armv6
    #elif defined(__ARM_ARCH_5TEJ__) \\
    || (defined(__TARGET_ARCH_ARM) && __TARGET_ARCH_ARM-0 >= 5)
    #error cmake_ARCH armv5
    #else
    #error cmake_ARCH arm
    #endif
    #elif defined(__i386) || defined(__i386__) || defined(_M_IX86)
    #error cmake_ARCH i386
    #elif defined(__x86_64) || defined(__x86_64__) || defined(__amd64) || defined(_M_X64)
    #error cmake_ARCH x86_64
    #elif defined(__ia64) || defined(__ia64__) || defined(_M_IA64)
    #error cmake_ARCH ia64
    #elif defined(__ppc__) || defined(__ppc) || defined(__powerpc__) \\
    || defined(_ARCH_COM) || defined(_ARCH_PWR) || defined(_ARCH_PPC) \\
    || defined(_M_MPPC) || defined(_M_PPC)
    #if defined(__ppc64__) || defined(__powerpc64__) || defined(__64BIT__)
    #error cmake_ARCH ppc64
    #else
    #error cmake_ARCH ppc
    #endif
    #endif
    #error cmake_ARCH unknown
    ")
    file(WRITE "${CMAKE_BINARY_DIR}/arch.c" "${archdetect_c_code}")
    enable_language(C)

    # Detect the architecture in a rather creative way...
    # This compiles a small C program which is a series of ifdefs that selects a
    # particular #error preprocessor directive whose message string contains the
    # target architecture. The program will always fail to compile (both because
    # file is not a valid C program, and obviously because of the presence of the
    # #error preprocessor directives... but by exploiting the preprocessor in this
    # way, we can detect the correct target architecture even when cross-compiling,
    # since the program itself never needs to be run (only the compiler/preprocessor)
    try_run(run_result_unused
            compile_result_unused
            "${CMAKE_BINARY_DIR}"
            "${CMAKE_BINARY_DIR}/arch.c"
            COMPILE_OUTPUT_VARIABLE ARCH
            CMAKE_FLAGS CMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES})

    # Parse the architecture name from the compiler output
    string(REGEX MATCH "cmake_ARCH ([a-zA-Z0-9_]+)" ARCH "${ARCH}")

    # Get rid of the value marker leaving just the architecture name
    string(REPLACE "cmake_ARCH " "" ARCH "${ARCH}")

    # If we are compiling with an unknown architecture this variable should
    # already be set to "unknown" but in the case that it's empty (i.e. due
    # to a typo in the code), then set it to unknown
    if (NOT ARCH)
      set(ARCH unknown)
    endif()
  endif()

  set(TargetArchitecture "${ARCH}" CACHE INTERNAL "")
endfunction()

# Gets and caches the target platform name
function(get_target_platform)
  if(WIN32)
    # See http://en.wikipedia.org/wiki/Comparison_of_Windows_versions
    if(CMAKE_SYSTEM_VERSION VERSION_EQUAL "6.2")
      # Windows 8
      set(Platform Win8)
    elseif(CMAKE_SYSTEM_VERSION VERSION_EQUAL "6.1")
      # Windows 7, Windows Server 2008 R2, Windows Home Server 2011
      set(Platform Win7)
    elseif(CMAKE_SYSTEM_VERSION VERSION_EQUAL "6.0")
      # Windows Server 2008
      set(Platform Vista)
    else()
      set(Platform Unsupported)
    endif()
  elseif(UNIX)
    if(APPLE)
      # See http://en.wikipedia.org/wiki/Darwin_%28operating_system%29
      if(CMAKE_SYSTEM_VERSION VERSION_LESS "12")
        set(Platform Unsupported)
      elseif(CMAKE_SYSTEM_VERSION VERSION_LESS "13")
        # OS X v10.8 "Mountain Lion"
        set(Platform OSX10.8)
      elseif(CMAKE_SYSTEM_VERSION VERSION_LESS "14")
        # OS X v10.9 "Mavericks"
        set(Platform OSX10.9)
      elseif(CMAKE_SYSTEM_VERSION VERSION_LESS "15")
        # OS X v10.10 "Yosemite"
        set(Platform OSX10.10)
      else()
        set(Platform Unsupported)
      endif()
    else()
      set(Platform Linux)
    endif()
  endif()
  set(TargetPlatform "${Platform}" CACHE INTERNAL "")
endfunction()

# Workaround for the Xcode's missing ability to pass -isystem to the compiler.
function(target_include_system_dirs Target)
  if(XCODE OR (UNIX AND NOT CMAKE_VERSION VERSION_LESS "3.0"))
    foreach(Arg ${ARGN})
      string(REGEX MATCH "\\$<" IsGeneratorExpression "${Arg}")
      if(Arg STREQUAL "PRIVATE" OR Arg STREQUAL "PUBLIC" OR Arg STREQUAL "INTERFACE")
        set(Scope ${Arg})
      elseif(NOT IsGeneratorExpression STREQUAL "")
        message(AUTHOR_WARNING "This function doesn't handle generator expressions; skipping ${Arg}")
      else()
        target_compile_options(${Target} ${Scope} -isystem${Arg})
      endif()
    endforeach()
  else()
    target_include_directories(${Target} SYSTEM ${Scope} ${ARGN})
  endif()
endfunction()

function(default_target_compile_options Target)
	target_compile_options(${Target}
	  PUBLIC
		$<$<BOOL:${MSVC}>:
			/W4      # Set warning level 4.
			#/WX      # Treat warnings as errors.
			/MP7     # Enable multi-processor compilation (max 7).
			/EHsc    # Catches C++ exceptions only and tells the compiler to assume that extern C functions never throw a C++ exception.
			/TP      # Treat sources as C++.
			/wd4351  # Disable C4351 'new behavior: elements of array 'array' will be default initialized'.
					 # Unneeded for new code (only applies to code previously compiled with VS 2005).
			/wd4503  # Disable C4503 'decorated name length exceeded' caused by boost multi-index and signals2.
					 # Disabled as per advice at https://svn.boost.org/trac/boost/wiki/Guidelines/WarningsGuidelines.
			/wd4512  # Disable C4512 'assignment operator could not be generated' caused by boost signals2.
					 # Disabled as per advice at http://lists.boost.org/boost-users/2009/01/44368.php.
			/wd4519  # Disable C4519 'default template arguments are only allowed on a class template'
			/wd4913  # Disable C4913 'default built-in binary operator ',' used' caused by inclusion of boost/utility/result_of.hpp.
					 # Disabled due to boost bug https://svn.boost.org/trac/boost/ticket/7663.
			/wd4996  # Disable C4996 'Function call with parameters that may be unsafe' caused by boost signals2.
					 # Disabled as per advice at https://svn.boost.org/trac/boost/wiki/Guidelines/WarningsGuidelines.
			/wd4714  # Disables warning generated by boost::expected library on function inlining.
			$<$<CONFIG:Release>:
				/O2  # Optimise code for maximum speed.  Implies the following:
					 #      Og (global optimisations)
					 #      Oi (replace some function calls with intrinsic functions),
					 #      Ot (favour fast code),
					 #      Oy (suppress creation of frame pointers on the call stack),
					 #      Ob2 (auto inline),
					 #      Gs (control stack probes),
					 #      GF (eliminate duplicate strings),
					 #      Gy (allows the compiler to package individual functions in the form of
					 #          packaged functions)
				/GL  # Whole program optimisation.
				/MD  # Use the multithread, dynamic version of the C run-time library.
			>
			$<$<CONFIG:Debug>:
				/Zi   # Produce a program database (.pdb) that contains type information and symbolic debugging information.
				/Od   # No optimizations in the program (speeds compilation).
				/RTC1 # Enables stack frame run-time error checking and checking for unintialised variables.
				/MDd  # Use the debug multithread, dynamic version of the C run-time library.
			>
			$<$<CONFIG:MinSizeRel>:
				/MD  # Use the multithread, dynamic version of the C run-time library.
			>
			$<$<CONFIG:RelWithDebInfo>:
				/O2  # Optimise code for maximum speed.
				/GL  # Whole program optimisation.
				/MD  # Use the multithread, dynamic version of the C run-time library.
				/Zi  # Produce a program database (.pdb) that contains type information and symbolic debugging information.
			>
		>
		$<$<BOOL:${UNIX}>:
			-std=c++11
			-pthread
			-W
			-fPIC
			-fstrict-overflow
			-pedantic
			-pedantic-errors
			$<$<CONFIG:Debug>:
				-O0
				-fno-inline
				-fno-eliminate-unused-debug-types
				-g3
				-ggdb
				${CoverageFlags}
			>
			$<$<CONFIG:Release>:-O2>
			$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
				${LibCXX}
				$<$<CONFIG:Debug>:
					-fdiagnostics-format=clang
					-fdiagnostics-show-option
					-fdiagnostics-fixit-info
					-Wno-unused-command-line-argument
				>
			>
			$<$<CXX_COMPILER_ID:GNU>:-static-libstdc++ -Wstrict-null-sentinel>
		>
	)
endfunction()
