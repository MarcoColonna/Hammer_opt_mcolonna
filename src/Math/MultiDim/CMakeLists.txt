message(STATUS "Preparing Hammer Tensor library.")

add_library(HammerTensor)
target_sources(HammerTensor
    PRIVATE
    SequentialIndexing.cc
    AlignedIndexing.cc
    BlockIndexing.cc
    BruteForceIterator.cc
    VectorContainer.cc
    SparseContainer.cc
    ScalarContainer.cc
    OuterContainer.cc
    Operations.cc
    ./Ops/Sum.cc
    ./Ops/Multiply.cc
    ./Ops/Divide.cc
    ./Ops/OuterSquare.cc
    ./Ops/Trace.cc
    ./Ops/Dot.cc
    ./Ops/DotOuterOptimizer.cc
    ./Ops/AddAt.cc
    ./Ops/Optimize.cc
    ./Ops/Convert.cc
    ./Ops/Compare.cc
)

set(HammerTensor_INCLUDE_DIR "../../../include/Hammer/Math/MultiDim")

target_sources(HammerTensor
    PRIVATE
    ${HammerTensor_INCLUDE_DIR}/../MultiDimensional.fhh
    ${HammerTensor_INCLUDE_DIR}/AlignedIndexing.hh
    ${HammerTensor_INCLUDE_DIR}/BinnedIndexing.hh
    ${HammerTensor_INCLUDE_DIR}/BinnedIndexingDefs.hh
    ${HammerTensor_INCLUDE_DIR}/BlockIndexing.hh
    ${HammerTensor_INCLUDE_DIR}/IContainer.hh
    ${HammerTensor_INCLUDE_DIR}/ISingleContainer.hh
    ${HammerTensor_INCLUDE_DIR}/LabeledIndexing.hh
    ${HammerTensor_INCLUDE_DIR}/LabeledIndexingDefs.hh
    ${HammerTensor_INCLUDE_DIR}/OperationDefs.hh
    ${HammerTensor_INCLUDE_DIR}/Operations.fhh
    ${HammerTensor_INCLUDE_DIR}/Operations.hh
    ${HammerTensor_INCLUDE_DIR}/OuterContainer.hh
    ${HammerTensor_INCLUDE_DIR}/ScalarContainer.hh
    ${HammerTensor_INCLUDE_DIR}/SequentialIndexing.hh
    ${HammerTensor_INCLUDE_DIR}/SparseContainer.hh
    ${HammerTensor_INCLUDE_DIR}/VectorContainer.hh
    ${HammerTensor_INCLUDE_DIR}/BruteForceIterator.hh
    ${HammerTensor_INCLUDE_DIR}/Ops/AddAt.hh
    ${HammerTensor_INCLUDE_DIR}/Ops/Compare.hh
    ${HammerTensor_INCLUDE_DIR}/Ops/Convert.hh
    ${HammerTensor_INCLUDE_DIR}/Ops/Divide.hh
    ${HammerTensor_INCLUDE_DIR}/Ops/Dot.hh
    ${HammerTensor_INCLUDE_DIR}/Ops/DotOuterOptimizer.hh
    ${HammerTensor_INCLUDE_DIR}/Ops/Multiply.hh
    ${HammerTensor_INCLUDE_DIR}/Ops/Optimize.hh
    ${HammerTensor_INCLUDE_DIR}/Ops/OuterSquare.hh
    ${HammerTensor_INCLUDE_DIR}/Ops/Sum.hh
    ${HammerTensor_INCLUDE_DIR}/Ops/Trace.hh
)

target_include_directories(HammerTensor PRIVATE ${PROJECT_SOURCE_DIR}/include)
set_target_properties(HammerTensor PROPERTIES
    POSITION_INDEPENDENT_CODE ${Hammer_PIC}
    VERSION ${Hammer_VERSION}
    SOVERSION ${Hammer_VERSION_MAJOR})
target_compile_features(HammerTensor PUBLIC cxx_std_17)

if(VERBOSE_DEBUG)
    target_compile_options(HammerTensor PRIVATE ${Hammer_VerboseOptions})
endif()

target_link_libraries(HammerTensor PUBLIC HammerTools PRIVATE Boost::headers)

if(SANITIZE)
    target_link_libraries(HammerTensor PRIVATE HammerSanitizers)
endif()

install(TARGETS HammerTensor
    EXPORT HammerTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT Hammer_Runtime
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
    COMPONENT Hammer_Runtime
    NAMELINK_COMPONENT Hammer_Development
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
    COMPONENT Hammer_Development
)
