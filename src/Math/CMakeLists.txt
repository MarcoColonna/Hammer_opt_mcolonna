message(STATUS "Preparing Hammer Math library.")

add_subdirectory(MultiDim)

add_library(HammerMath)
target_sources(HammerMath
    PRIVATE
    Tensor.cc
    Integrator.cc
    Histogram.cc
    HistogramDefinition.cc
    HistogramSet.cc
    FourMomentum.cc
    PhaseSpace.cc
    Units.cc
    Utils.cc
)

set(HammerMath_INCLUDE_DIR "../../include/Hammer/Math")

target_sources(HammerMath
    PRIVATE
    ${HammerMath_INCLUDE_DIR}/Constants.hh
    ${HammerMath_INCLUDE_DIR}/Differentiator.hh
    ${HammerMath_INCLUDE_DIR}/FourMomentum.hh
    ${HammerMath_INCLUDE_DIR}/FourMomentumDefs.hh
    ${HammerMath_INCLUDE_DIR}/Histogram.hh
    ${HammerMath_INCLUDE_DIR}/HistogramDefinition.hh
    ${HammerMath_INCLUDE_DIR}/HistogramDefs.hh
    ${HammerMath_INCLUDE_DIR}/HistogramSet.hh
    ${HammerMath_INCLUDE_DIR}/Integrator.fhh
    ${HammerMath_INCLUDE_DIR}/Integrator.hh
    ${HammerMath_INCLUDE_DIR}/PhaseSpace.hh
    ${HammerMath_INCLUDE_DIR}/Tensor.hh
    ${HammerMath_INCLUDE_DIR}/TensorDefs.hh
    ${HammerMath_INCLUDE_DIR}/Units.hh
    ${HammerMath_INCLUDE_DIR}/Utils.hh
)

# pre-Cmake-3.23 lengthy commands
set(Math_Headers_Inst
    ${HammerMath_INCLUDE_DIR}/FourMomentum.hh
    ${HammerMath_INCLUDE_DIR}/FourMomentumDefs.hh
    ${HammerMath_INCLUDE_DIR}/Utils.hh
)
separate_arguments(Math_Headers_Inst)
set_target_properties(HammerMath PROPERTIES
    PUBLIC_HEADER
    "${Math_Headers_Inst}"
)

target_include_directories(HammerMath PRIVATE ${PROJECT_SOURCE_DIR}/include)
set_target_properties(HammerMath PROPERTIES
    POSITION_INDEPENDENT_CODE ${Hammer_PIC}
    VERSION ${Hammer_VERSION}
    SOVERSION ${Hammer_VERSION_MAJOR})
target_compile_features(HammerMath PUBLIC cxx_std_17)

if(VERBOSE_DEBUG)
    target_compile_options(HammerMath PRIVATE ${Hammer_VerboseOptions})
endif()

set(MATH_TARGETS HammerConfig HammerTools Boost::headers HammerTensor)

if(HAVE_ROOT)
    set(MATH_TARGETS "${MATH_TARGETS} ROOT::Physics")
endif()

separate_arguments(MATH_TARGETS)
target_link_libraries(HammerMath PUBLIC ${MATH_TARGETS})

if(SANITIZE)
    target_link_libraries(HammerMath PRIVATE HammerSanitizers)
endif()

install(TARGETS HammerMath
    EXPORT HammerTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT Hammer_Runtime
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
    COMPONENT Hammer_Runtime
    NAMELINK_COMPONENT Hammer_Development
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
    COMPONENT Hammer_Development
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Hammer/Math
    COMPONENT Hammer_Development
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

