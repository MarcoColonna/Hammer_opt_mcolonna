# if(${ROOT_CXX_STANDARD} LESS ${MAX_CXX_STANDARD})
# if(${ROOT_CXX_STANDARD} LESS 14)
# MESSAGE(STATUS "ROOT seems to be compiled with c++11. Hammer will be compiled with c++14.")
# set(MAX_CXX_STD 14)
# else()
# set(MAX_CXX_STD ${ROOT_CXX_STANDARD})
# endif()
# endif()

if(${ROOT_VERSION} VERSION_LESS 6.16)
    # This is required for ROOT < 6.16
    string(REPLACE "-L " "-L" ROOT_EXE_LINKER_FLAGS "${ROOT_EXE_LINKER_FLAGS}")
endif()

# This is required on if there is more than one flag (like on macOS)
separate_arguments(ROOT_EXE_LINKER_FLAGS)

# -bind_at_load and -undefined,dynamic_lookup are Apple ld flags for Cling/cppyy
if(APPLE)
    set(_LLXX "-Wl,-w;-Wl,-bind_at_load;-Wl,-undefined,dynamic_lookup")
else()
    set(_LLXX "")
endif()
set(_CSTDXX "cxx_std_${ROOT_CXX_STANDARD}")

# Create imported target ROOT::Core
if(NOT TARGET ROOT::Core)
    add_library(ROOT::Core SHARED IMPORTED)

    set_target_properties(ROOT::Core PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libCore.so"
        IMPORTED_SONAME "@rpath/libCore.so"
    )
endif()

if(NOT TARGET ROOT::Cling)
    # Create imported target ROOT::Cling
    add_library(ROOT::Cling SHARED IMPORTED)

    set_target_properties(ROOT::Cling PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${_LLXX}"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libCling.so"
        IMPORTED_SONAME "@rpath/libCling.so"
    )
endif()

# Create imported target ROOT::Physics
if(NOT TARGET ROOT::Physics)
    add_library(ROOT::Physics SHARED IMPORTED)

    set_target_properties(ROOT::Physics PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::Matrix;ROOT::MathCore;ROOT::GenVector"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libPhysics.so"
        IMPORTED_SONAME "@rpath/libPhysics.so"
    )
endif()

# Create imported target ROOT::Matrix
if(NOT TARGET ROOT::Matrix)
    add_library(ROOT::Matrix SHARED IMPORTED)

    set_target_properties(ROOT::Matrix PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::MathCore"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libMatrix.so"
        IMPORTED_SONAME "@rpath/libMatrix.so"
    )
endif()

# Create imported target ROOT::MathCore
if(NOT TARGET ROOT::MathCore)
    add_library(ROOT::MathCore SHARED IMPORTED)

    set_target_properties(ROOT::MathCore PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::Core;ROOT::Imt"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libMathCore.so"
        IMPORTED_SONAME "@rpath/libMathCore.so"
    )
endif()

# Create imported target ROOT::Imt
if(NOT TARGET ROOT::Imt)
    add_library(ROOT::Imt SHARED IMPORTED)

    set_target_properties(ROOT::Imt PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::MultiProc;ROOT::Core"
        IMPORTED_LINK_DEPENDENT_LIBRARIES "ROOT::Thread"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libImt.so"
        IMPORTED_SONAME "@rpath/libImt.so"
    )
endif()

# Create imported target ROOT::GenVector
if(NOT TARGET ROOT::GenVector)
    add_library(ROOT::GenVector SHARED IMPORTED)

    set_target_properties(ROOT::GenVector PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::Core;ROOT::MathCore"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libGenVector.so"
        IMPORTED_SONAME "@rpath/libGenVector.so"
    )
endif()

# Create imported target ROOT::Hist
if(NOT TARGET ROOT::Hist)
    add_library(ROOT::Hist SHARED IMPORTED)

    set_target_properties(ROOT::Hist PROPERTIES
        INTERFACE_COMPILE_DEFINITIONS "ROOT_SUPPORT_CLAD"
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::MathCore;ROOT::Matrix;ROOT::RIO"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libHist.so"
        IMPORTED_SONAME "@rpath/libHist.so"
    )
endif()

# Create imported target ROOT::RIO
if(NOT TARGET ROOT::RIO)
    add_library(ROOT::RIO SHARED IMPORTED)

    set_target_properties(ROOT::RIO PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::Core;ROOT::Thread"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libRIO.so"
        IMPORTED_SONAME "@rpath/libRIO.so"
    )
endif()

# Create imported target ROOT::Thread
if(NOT TARGET ROOT::Thread)
    add_library(ROOT::Thread SHARED IMPORTED)

    set_target_properties(ROOT::Thread PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::Core"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libThread.so"
        IMPORTED_SONAME "@rpath/libThread.so"
    )
endif()

# Create imported target ROOT::MultiProc
if(NOT TARGET ROOT::MultiProc)
    add_library(ROOT::MultiProc SHARED IMPORTED)

    set_target_properties(ROOT::MultiProc PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::Core;ROOT::Net"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libMultiProc.so"
        IMPORTED_SONAME "@rpath/libMultiProc.so"
    )
endif()

# Create imported target ROOT::Net
if(NOT TARGET ROOT::Net)
    add_library(ROOT::Net SHARED IMPORTED)

    set_target_properties(ROOT::Net PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::RIO"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libNet.so"
        IMPORTED_SONAME "@rpath/libNet.so"
    )
endif()

# Create imported target ROOT::rootcling
if(NOT TARGET ROOT::rootcling)
    if(DEFINED ROOT_rootcling_CMD AND ROOT_rootcling_CMD)
        set(_rootcling_exe "${ROOT_rootcling_CMD}")
    elseif(DEFINED ROOTSYS)
        set(_rootcling_exe "${ROOTSYS}/bin/rootcling")
    else()
        find_program(_rootcling_exe NAMES rootcling HINTS "${ROOT_LIBRARY_DIR}/../bin")
    endif()

    add_executable(ROOT::rootcling IMPORTED)
    set_property(TARGET ROOT::rootcling PROPERTY ENABLE_EXPORTS 1)
    set_target_properties(ROOT::rootcling PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_LINK_LIBRARIES "ROOT::RIO;ROOT::Cling;ROOT::Core;ROOT::Rint"
        IMPORTED_LOCATION "${_rootcling_exe}"
    )
    unset(_rootcling_exe)
endif()

# Create imported target ROOT::cppyy_backend
if(NOT TARGET ROOT::cppyy_backend)
    add_library(ROOT::cppyy_backend SHARED IMPORTED)

    set_target_properties(ROOT::cppyy_backend PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_LINK_LIBRARIES "ROOT::Core"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libcppyy_backend.so"
        IMPORTED_SONAME "@rpath/libcppyy_backend.so"
    )
endif()

# Create imported target ROOT::cppyy
if(NOT TARGET ROOT::cppyy)
    add_library(ROOT::cppyy SHARED IMPORTED)

    set_target_properties(ROOT::cppyy PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_LINK_LIBRARIES "${_LLXX};ROOT::cppyy_backend"
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${Python3_INCLUDE_DIRS}"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libcppyy.so"
        IMPORTED_SONAME "@rpath/libcppyy.so"
    )
endif()

# Create imported target ROOT::ROOTPythonizations
if(NOT TARGET ROOT::ROOTPythonizations)
    add_library(ROOT::ROOTPythonizations SHARED IMPORTED)

    set_target_properties(ROOT::ROOTPythonizations PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_LINK_LIBRARIES "${_LLXX};ROOT::Core;ROOT::Tree;ROOT::cppyy"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libROOTPythonizations.so"
        IMPORTED_SONAME "@rpath/libROOTPythonizations.so"
    )
endif()

# Create imported target ROOT::PyROOT
if(NOT TARGET ROOT::PyROOT)
    add_library(ROOT::PyROOT INTERFACE IMPORTED)

    set_target_properties(ROOT::PyROOT PROPERTIES
        INTERFACE_LINK_LIBRARIES "ROOT::cppyy_backend;ROOT::cppyy;ROOT::ROOTPythonizations"
    )
endif()

# Create imported target ROOT::ROOTTPython
if(NOT TARGET ROOT::ROOTTPython)
    add_library(ROOT::ROOTTPython SHARED IMPORTED)

    set_target_properties(ROOT::ROOTTPython PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::Core"
        IMPORTED_LINK_DEPENDENT_LIBRARIES "ROOT::cppyy;Python3::Python"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libROOTTPython.so"
        IMPORTED_SONAME "@rpath/libROOTTPython.so"
    )
endif()

# Create imported target ROOT::Rint
if(NOT TARGET ROOT::Rint)
    add_library(ROOT::Rint SHARED IMPORTED)

    set_target_properties(ROOT::Rint PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::Core"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libRint.so"
        IMPORTED_SONAME "@rpath/libRint.so"
    )
endif()

# Create imported target ROOT::Tree
if(NOT TARGET ROOT::Tree)
    add_library(ROOT::Tree SHARED IMPORTED)

    set_target_properties(ROOT::Tree PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::Imt;ROOT::Net;ROOT::RIO;ROOT::MathCore"
        IMPORTED_LOCATION "${ROOT_LIBRARY_DIR}/libTree.so"
        IMPORTED_SONAME "@rpath/libTree.so"
    )
endif()

# Create imported target ROOT::Gpad
if(NOT TARGET ROOT::Gpad)
    add_library(ROOT::Gpad SHARED IMPORTED)

    set_target_properties(ROOT::Gpad PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::Graf;ROOT::Hist"
        IMPORTED_LOCATION_RELEASE "${ROOT_LIBRARY_DIR}/libGpad.so"
        IMPORTED_SONAME_RELEASE "@rpath/libGpad.so"
    )
endif()

# Create imported target ROOT::Graf
if(NOT TARGET ROOT::Graf)
    add_library(ROOT::Graf SHARED IMPORTED)

    set_target_properties(ROOT::Graf PROPERTIES
        INTERFACE_COMPILE_FEATURES ${_CSTDXX}
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "ROOT::Hist;ROOT::Matrix;ROOT::MathCore;ROOT::RIO"
        IMPORTED_LOCATION_RELEASE "${ROOT_LIBRARY_DIR}/libGraf.so"
        IMPORTED_SONAME_RELEASE "@rpath/libGraf.so"
    )
endif()

# Mark ROOT include directories as SYSTEM to suppress compiler warnings from ROOT headers.
# This applies to targets already exported by ROOT's CMake config (CONFIG path).
# Manually created targets above already use INTERFACE_SYSTEM_INCLUDE_DIRECTORIES directly.
set(_hammer_root_patch_targets
    ROOT::Core
    ROOT::Hist
    ROOT::Physics
    ROOT::MathCore
    ROOT::Matrix
    ROOT::GenVector
    ROOT::RIO
    ROOT::Thread
    ROOT::Imt
    ROOT::MultiProc
    ROOT::Net
    ROOT::Tree
    ROOT::Rint
    ROOT::Gpad
    ROOT::Graf
    ROOT::Cling
    ROOT::cppyy
    ROOT::cppyy_backend
    ROOT::ROOTPythonizations
    ROOT::ROOTTPython
    ROOT::PyROOT
    ROOT::rootcling
)
foreach(_hammer_tgt IN LISTS _hammer_root_patch_targets)
    if(TARGET "${_hammer_tgt}")
        get_target_property(_hammer_tgt_incs "${_hammer_tgt}" INTERFACE_INCLUDE_DIRECTORIES)
        if(_hammer_tgt_incs AND NOT "${_hammer_tgt_incs}" MATCHES "NOTFOUND")
            target_include_directories("${_hammer_tgt}" SYSTEM INTERFACE "${_hammer_tgt_incs}")
        endif()
    endif()
endforeach()
unset(_hammer_root_patch_targets)
unset(_hammer_tgt)
unset(_hammer_tgt_incs)
