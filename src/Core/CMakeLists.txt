message(STATUS "Preparing Hammer Base library.")

add_library(HammerBase)
target_sources(HammerBase
	PRIVATE
	Particle.cc
	AmplitudeBase.cc
	FormFactorBase.cc
	RateBase.cc
)

set(HammerBase_INCLUDE_DIR "../../include/Hammer")

target_sources(HammerBase
	PRIVATE
	${HammerBase_INCLUDE_DIR}/Particle.fhh
	${HammerBase_INCLUDE_DIR}/Particle.hh
	${HammerBase_INCLUDE_DIR}/ParticleDefs.hh
	${HammerBase_INCLUDE_DIR}/AmplitudeBase.hh
	${HammerBase_INCLUDE_DIR}/FormFactorBase.hh
	${HammerBase_INCLUDE_DIR}/RateBase.hh
)

# pre-Cmake-3.23 lengthy commnands
set(Base_Headers_Inst
	${HammerBase_INCLUDE_DIR}/Particle.fhh
	${HammerBase_INCLUDE_DIR}/Particle.hh
	${HammerBase_INCLUDE_DIR}/ParticleDefs.hh
)
separate_arguments(Base_Headers_Inst)
set_target_properties(HammerBase PROPERTIES
	PUBLIC_HEADER
	"${Base_Headers_Inst}"
)

target_include_directories(HammerBase PRIVATE ${PROJECT_SOURCE_DIR}/include)

set_target_properties(HammerBase PROPERTIES
	POSITION_INDEPENDENT_CODE ${Hammer_PIC}
	VERSION ${Hammer_VERSION}
	SOVERSION ${Hammer_VERSION_MAJOR})
target_compile_features(HammerBase PUBLIC cxx_std_17)

if(VERBOSE_DEBUG)
	target_compile_options(HammerBase PRIVATE ${Hammer_VerboseOptions})
endif()

target_link_libraries(HammerBase PUBLIC HammerTools HammerMath Boost::headers yaml-cpp::yaml-cpp)

if(SANITIZE)
	target_link_libraries(HammerBase PRIVATE HammerSanitizers)
endif()

install(TARGETS HammerBase
	EXPORT HammerTargets
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	COMPONENT Hammer_Runtime
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
	COMPONENT Hammer_Runtime
	NAMELINK_COMPONENT Hammer_Development
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
	COMPONENT Hammer_Development
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Hammer
	COMPONENT Hammer_Development
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

message(STATUS "Preparing Hammer Core library.")

add_library(HammerCore)
target_sources(HammerCore
	PRIVATE
	Process.cc
	ProcResults.cc
	ProcGraph.cc
	DictionaryManager.cc
	PurePhaseSpaceDefs.cc
	SchemeDefinitions.cc
	ProcessDefinitions.cc
	ExternalData.cc
	WCSpecialization.cc
	SpecializationDefinitions.cc
	ProcRates.cc
	ProvidersRepo.cc
	ProcRequirements.cc
	ProcManager.cc
	Histos.cc
	Event.cc
	FFEval.cc
	Hammer.cc
)

set(HammerCore_INCLUDE_DIR "../../include/Hammer")

target_sources(HammerCore
	PRIVATE
	${HammerCore_INCLUDE_DIR}/Hammer.hh
	${HammerCore_INCLUDE_DIR}/Process.hh
	${HammerCore_INCLUDE_DIR}/DictionaryManager.hh
	${HammerCore_INCLUDE_DIR}/Event.hh
	${HammerCore_INCLUDE_DIR}/ExternalData.hh
	${HammerCore_INCLUDE_DIR}/FFEval.hh
	${HammerCore_INCLUDE_DIR}/HammerDoc.hh
	${HammerCore_INCLUDE_DIR}/Histos.fhh
	${HammerCore_INCLUDE_DIR}/Histos.hh
	${HammerCore_INCLUDE_DIR}/ProcessDefinitions.hh
	${HammerCore_INCLUDE_DIR}/ProcGraph.fhh
	${HammerCore_INCLUDE_DIR}/ProcGraph.hh
	${HammerCore_INCLUDE_DIR}/ProcManager.hh
	${HammerCore_INCLUDE_DIR}/ProcRates.hh
	${HammerCore_INCLUDE_DIR}/ProcRequirements.hh
	${HammerCore_INCLUDE_DIR}/ProcResults.hh
	${HammerCore_INCLUDE_DIR}/ProvidersRepo.hh
	${HammerCore_INCLUDE_DIR}/PurePhaseSpaceDefs.hh
	${HammerCore_INCLUDE_DIR}/SchemeDefinitions.hh
	${HammerCore_INCLUDE_DIR}/SpecializationDefinitions.hh
	${HammerCore_INCLUDE_DIR}/WCSpecialization.hh
)

# pre-Cmake-3.23 lengthy commnands
set(Core_Headers_Inst
	${HammerCore_INCLUDE_DIR}/Process.hh
	${HammerCore_INCLUDE_DIR}/Hammer.hh
)
separate_arguments(Core_Headers_Inst)
set_target_properties(HammerCore PROPERTIES
	PUBLIC_HEADER
	"${Core_Headers_Inst}"
)

target_include_directories(HammerCore PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)
set_target_properties(HammerCore PROPERTIES
	POSITION_INDEPENDENT_CODE ${Hammer_PIC}
	VERSION ${Hammer_VERSION}
	SOVERSION ${Hammer_VERSION_MAJOR})
target_compile_features(HammerCore PUBLIC cxx_std_17)

if(VERBOSE_DEBUG)
	target_compile_options(HammerCore PRIVATE ${Hammer_VerboseOptions})
endif()

set(CORE_TARGETS "HammerConfig HammerTools HammerBase HammerMath HammerTensor Amplitudes FormFactors Rates Boost::thread")

if(HAVE_ROOT)
	set(CORE_TARGETS "${CORE_TARGETS} ROOT::Physics ROOT::Hist")
endif()

separate_arguments(CORE_TARGETS)

target_link_libraries(HammerCore PUBLIC ${CORE_TARGETS})

if(SANITIZE)
	target_link_libraries(HammerCore PRIVATE HammerSanitizers)
endif()

install(TARGETS HammerCore
	EXPORT HammerTargets
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	COMPONENT Hammer_Runtime
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
	COMPONENT Hammer_Runtime
	NAMELINK_COMPONENT Hammer_Development
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
	COMPONENT Hammer_Development
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Hammer
	COMPONENT Hammer_Development
)
