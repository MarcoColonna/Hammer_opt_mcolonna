message(STATUS "Preparing Include Files.")

SET(FIXROOTINCL "")

if(APPLE AND WITH_ROOT)
    get_property(dirs TARGET HammerCore PROPERTY INCLUDE_DIRECTORIES)
    get_property(tmpdirs1 TARGET yaml-cpp::yaml-cpp PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    list(APPEND dirs ${tmpdirs1})
    get_property(tmpdirs2 TARGET ROOT::Hist PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    list(APPEND dirs ${tmpdirs2})

    if(WITH_EXAMPLES)
        get_property(tmpdirs3 TARGET HepMC3::All PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
        list(APPEND dirs ${tmpdirs3})
    endif()

    foreach(inclpath ${dirs})
        if(EXISTS "${inclpath}/root/ROOT/RConfig.h")
            MESSAGE(STATUS "There may be a ROOT include problem -- trying to fix it.")
            SET(FIXROOTINCL "root/ROOT/")
            break()
        endif()
    endforeach()
endif()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Hammer/Config/HammerConfig.hh.in
    ${PROJECT_BINARY_DIR}/include/Hammer/Config/HammerConfig.hh
    ESCAPE_QUOTES)

file(GLOB AVAIL_AMPLS LIST_DIRECTORIES FALSE RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/Hammer/Amplitudes" "${CMAKE_CURRENT_SOURCE_DIR}/Hammer/Amplitudes/*.hh")
file(GLOB_RECURSE AVAIL_FFS LIST_DIRECTORIES FALSE RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/Hammer/FormFactors" "${CMAKE_CURRENT_SOURCE_DIR}/Hammer/FormFactors/**/*.hh")
file(GLOB AVAIL_RATES LIST_DIRECTORIES FALSE RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/Hammer/Rates" "${CMAKE_CURRENT_SOURCE_DIR}/Hammer/Rates/*.hh")

set(HAMMER_AMPLITUDES_INCLUDES "")
set(AMPL_BODY "")

foreach(FNAME ${AVAIL_AMPLS})
    if(NOT(FNAME MATCHES ".*Base\.hh"))
        get_filename_component(FNAME_NOEXT ${FNAME} NAME_WE)
        string(APPEND HAMMER_AMPLITUDES_INCLUDES "#include \"Hammer/Amplitudes/" "${FNAME}" "\"\n")
        string(APPEND AMPL_BODY "        ADD(" "${FNAME_NOEXT}" ", CONT, AmplitudeBase); \\\n")
    endif()
endforeach()

set(HAMMER_FORMFACTORS_INCLUDES "")
set(FF_BODY "")

foreach(FNAME ${AVAIL_FFS})
    if((NOT(FNAME MATCHES ".*Base\.hh")) AND(FNAME MATCHES "FF.*\.hh"))
        get_filename_component(FNAME_NOEXT ${FNAME} NAME_WE)
        string(APPEND HAMMER_FORMFACTORS_INCLUDES "#include \"Hammer/FormFactors/" "${FNAME}" "\"\n")
        string(APPEND FF_BODY "        ADDV(" "${FNAME_NOEXT}" ", CONT, FormFactorBase); \\\n")
    endif()
endforeach()

set(HAMMER_RATES_INCLUDES "")
set(RATE_BODY "")

foreach(FNAME ${AVAIL_RATES})
    if(NOT(FNAME MATCHES ".*Base\.hh"))
        get_filename_component(FNAME_NOEXT ${FNAME} NAME_WE)
        string(APPEND HAMMER_RATES_INCLUDES "#include \"Hammer/Rates/" "${FNAME}" "\"\n")
        string(APPEND RATE_BODY "        ADD(" "${FNAME_NOEXT}" ", CONT, RateBase); \\\n")
    endif()
endforeach()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Hammer/Config/HammerModules.hh.in
    ${PROJECT_BINARY_DIR}/include/Hammer/Config/HammerModules.hh
)

add_library(HammerConfig INTERFACE)
target_sources(HammerConfig
    PRIVATE
    ${PROJECT_BINARY_DIR}/include/Hammer/Config/HammerModules.hh
    ${PROJECT_BINARY_DIR}/include/Hammer/Config/HammerConfig.hh
)
set_target_properties(HammerConfig PROPERTIES
    PUBLIC_HEADER
    "${PROJECT_BINARY_DIR}/include/Hammer/Config/HammerModules.hh;${PROJECT_BINARY_DIR}/include/Hammer/Config/HammerConfig.hh"
)
target_include_directories(HammerConfig INTERFACE
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)

if(HAVE_ROOT)
    target_link_libraries(HammerCore PUBLIC ROOT::Core)
endif()

install(TARGETS HammerConfig
    EXPORT HammerTargets
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Hammer/Config
    COMPONENT Hammer_Development
)

# Tier 2: additional headers needed for the ROOT/Python reflection workflow
if(WITH_ROOT OR INSTALL_HEADERS)
    install(FILES
        ${PROJECT_SOURCE_DIR}/include/Hammer/Tools/HammerRoot.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Hammer/Tools
        COMPONENT Hammer_Development
    )
endif()

if(WITH_PYTHON OR WITH_ROOT OR INSTALL_HEADERS)
    install(FILES
        ${PROJECT_SOURCE_DIR}/include/Hammer/Tools/HammerSerial.hh
        ${PROJECT_SOURCE_DIR}/include/Hammer/Tools/HammerSerial_generated.h
        ${PROJECT_SOURCE_DIR}/include/Hammer/Tools/ParticleData.hh
        ${PROJECT_SOURCE_DIR}/include/Hammer/Tools/ParticleUtils.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Hammer/Tools
        COMPONENT Hammer_Development
    )
    install(DIRECTORY
        ${PROJECT_SOURCE_DIR}/include/flatbuffers
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT Hammer_Development
    )
endif()

# Tier 3: all development headers — every header under include/Hammer/ in the
# source tree (Amplitudes, FormFactors, Rates, Math/MultiDim, and all
# internal headers) plus the flatbuffers bundled dependency.
if(INSTALL_HEADERS)
    install(DIRECTORY
        ${PROJECT_SOURCE_DIR}/include/Hammer
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT Hammer_Development
        FILES_MATCHING
            PATTERN "*.hh"
            PATTERN "*.fhh"
            PATTERN "*.h"
    )
    # flatbuffers already covered by Tier 2 when INSTALL_HEADERS=ON, but
    # listed here explicitly in case Tier 2 conditions are not met.
    install(DIRECTORY
        ${PROJECT_SOURCE_DIR}/include/flatbuffers
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT Hammer_Development
    )
endif()
