# setup verbosity level for compile flags
# set(Hammer_CompileOptions "")
set(Hammer_VerboseOptions "")
set(Hammer_SanitizeOptions "")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -O0 -fno-optimize-sibling-calls -fno-omit-frame-pointer")
    set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -Wall -Wextra -Weffc++ -Wconversion -Wsign-conversion -Wpedantic -Wno-c++98-compat")
    set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -Wno-c++98-compat-pedantic -Wno-padded -Wno-documentation-deprecated-sync")
    set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -Wno-unknown-pragmas -Wno-weak-vtables")
    set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -fstack-protector-strong -Wstack-protector")
    set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -Wno-documentation -Wno-documentation-unknown-command -Wno-poison-system-directories")
    set(Hammer_SanitizeOptions "${Hammer_SanitizeOptions} -fsanitize=address -fsanitize=undefined -fno-sanitize=vptr")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "Detected GNU compiler.")

    set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -O0 -fno-optimize-sibling-calls -fno-omit-frame-pointer")
    set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -fno-sanitize-recover -fstack-protector")
    set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -pedantic -Wall -Wextra -Weffc++ -Wshadow -Wformat=2 -Wfloat-equal -Wconversion -Wsign-conversion -Wlogical-op")
    set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -Wcast-qual -Wcast-align -Wno-unknown-pragmas -Wno-effc++")

    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6.0.0)
        set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -Wshift-overflow=2 -Wduplicated-cond")
    endif()

    set(Hammer_SanitizeOptions "${Hammer_SanitizeOptions} -fsanitize=address -fsanitize=undefined")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
# using Intel C++
else()
    set(Hammer_VerboseOptions "${Hammer_VerboseOptions} -O0")
    set(Hammer_SanitizeOptions "${Hammer_SanitizeOptions}")
endif()

# separate_arguments(Hammer_CompileOptions)
separate_arguments(Hammer_VerboseOptions)
separate_arguments(Hammer_SanitizeOptions)

if(SANITIZE)
    add_library(HammerSanitizers INTERFACE)
    target_compile_options(HammerSanitizers INTERFACE ${Hammer_SanitizeOptions})
    target_link_options(HammerSanitizers INTERFACE ${Hammer_SanitizeOptions})

    install(TARGETS HammerSanitizers
            EXPORT HammerTargets
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            COMPONENT Hammer_Runtime
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
            COMPONENT Hammer_Runtime
            NAMELINK_COMPONENT Hammer_Development
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
            COMPONENT Hammer_Development
        )
endif()

# setting fPIC
message(STATUS "${CMAKE_SYSTEM_PROCESSOR}")

IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" OR CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
    set(Hammer_PIC TRUE)
ELSE()
    set(Hammer_PIC ${CMAKE_POSITION_INDEPENDENT_CODE})
ENDIF()

message(STATUS "Process independent code flag is ${Hammer_PIC}.")

IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
ENDIF()

if(BUILD_SHARED_LIBS)
    # setting correct rpaths
    if(APPLE)
        SET(CMAKE_MACOSX_RPATH TRUE)
    endif()

    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    if(NOT APPLE)
        set(CMAKE_INSTALL_RPATH $ORIGIN)
    else()
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}/Hammer")
    endif()

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# platform information
set(L2_Cache_Width 64)

if(UNIX)
    if(APPLE)
        execute_process(COMMAND sysctl -a machdep.cpu.cache.linesize
            COMMAND sed "s/machdep.cpu.cache.linesize: //"
            OUTPUT_VARIABLE L2_Cache_Width)
    else()
        execute_process(COMMAND getconf LEVEL1_DCACHE_LINESIZE
            OUTPUT_VARIABLE L2_Cache_Width)
    endif()
endif(UNIX)

if(CMAKE_VERSION VERSION_LESS 3.20.0)
    include(TestBigEndianCXX)
    TEST_BIG_ENDIAN_CXX(IS_BIG_ENDIAN)
else()
    if(CMAKE_CXX_BYTE_ORDER STREQUAL "BIG_ENDIAN")
        set(IS_BIG_ENDIAN)
    endif()
endif()
