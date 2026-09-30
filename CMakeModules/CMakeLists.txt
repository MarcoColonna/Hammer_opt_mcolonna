install(SCRIPT "PostInstallMessages.cmake")

set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Hammer)

# set(INCLUDE_DIRS ${CMAKE_INSTALL_FULL_INCLUDEDIR})
# set(HAMMER_LIBDIR ${CMAKE_INSTALL_FULL_LIBDIR})
# get_property(LIB_INSTALL_YAML TARGET yaml-cpp::yaml-cpp PROPERTY INTERFACE_LINK_LIBRARIES)

# if(HAVE_ROOT)
# get_property(INCLUDE_DIRS_ROOT TARGET ROOT::Core PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
# get_property(LIB_ROOT TARGET ROOT::Core PROPERTY INTERFACE_LINK_LIBRARIES)
# endif()
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    HammerConfigVersion.cmake
    VERSION ${HAMMER_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(EXPORT HammerTargets
    FILE HammerTargets.cmake
    NAMESPACE Hammer::
    DESTINATION ${INSTALL_CONFIGDIR}
)

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
set(DATA_INSTALL_DIR ${CMAKE_INSTALL_DATADIR}/Hammer)

configure_package_config_file(HammerConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/HammerConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
    PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR DATA_INSTALL_DIR
)

# Build substitution blocks for HammerConfigDependencies.cmake
if(SANITIZE)
    set(HAMMER_FIND_SANITIZERS "find_dependency(Sanitizers)")
else()
    set(HAMMER_FIND_SANITIZERS "")
endif()

if(HAVE_ROOT)
    if(ROOT_Found_HOW STREQUAL "New")
        string(CONCAT HAMMER_FIND_ROOT
            "find_dependency(ROOT CONFIG COMPONENTS ${ROOT_COMPONENTS_LIST})\n"
            [=[include("${CMAKE_CURRENT_LIST_DIR}/FixROOT.cmake")]=])
    else()
        set(HAMMER_FIND_ROOT [=[list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_dependency(ROOT)
include("${CMAKE_CURRENT_LIST_DIR}/FixROOT.cmake")]=])
    endif()
else()
    set(HAMMER_FIND_ROOT "")
endif()

configure_file(HammerConfigDependencies.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/HammerConfigDependencies.cmake @ONLY)

# Install the config, configversion and find-module helpers
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/HammerConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/HammerConfigVersion.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/HammerConfigDependencies.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

if(HAVE_ROOT)
    install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/FixROOT.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
    )
endif()

# pkg-config file
set(HAMMER_PC_BOOST_CFLAGS "-I${Boost_INCLUDE_DIRS}")
set(_blibs "${Boost_THREAD_LIBRARY}")
string(REPLACE ";" " " HAMMER_PC_BOOST_LIBS "${_blibs}")
unset(_blibs)

if(HAVE_ROOT)
    set(_root_libs ${ROOT_COMPONENTS_LIST})
    # Add the transitive core ROOT libraries that the above components pull in
    list(APPEND _root_libs Core MathCore Matrix GenVector RIO Thread Imt MultiProc Net)
    list(REMOVE_DUPLICATES _root_libs)
    set(_root_lflags)
    foreach(_lib IN LISTS _root_libs)
        list(APPEND _root_lflags "-l${_lib}")
    endforeach()
    string(JOIN " " HAMMER_PC_ROOT_LIBS "-L${ROOT_LIBRARY_DIR}" ${_root_lflags})
    unset(_root_libs)
    unset(_root_lflags)
else()
    set(HAMMER_PC_ROOT_LIBS "")
endif()

configure_file(Hammer.pc.in ${CMAKE_CURRENT_BINARY_DIR}/Hammer.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Hammer.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
    COMPONENT Hammer_Development
)

# Install the legacy ROOT finder only when it was used at build time (ROOT < 6.16).
# Downstream find_package(Hammer) will invoke it as FindROOT.cmake via MODULE mode.
if(HAVE_ROOT AND NOT "${ROOT_Found_HOW}" STREQUAL "New")
    install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/FindROOT_old.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        RENAME FindROOT.cmake
    )
endif()
