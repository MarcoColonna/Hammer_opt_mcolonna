message(STATUS "Preparing Hammer Tools library.")

add_library(HammerTools)

target_sources(HammerTools
	PRIVATE
	Setting.cc
	SettingsConsumer.cc
	Logging.cc
	ParticleData.cc
	ParticleUtils.cc
	Pdg.cc
	IOTypes.cc
	Errors.cc
	HammerSerial.cc
	Utils.cc
	../Core/SettingsHandler.cc
	../Core/IndexLabels.cc
)

set(HammerTools_INCLUDE_DIR "../../include/Hammer/Tools")

target_sources(HammerTools
	PRIVATE
	${HammerTools_INCLUDE_DIR}/HammerRoot.fhh
	${HammerTools_INCLUDE_DIR}/HammerRoot.hh
	${HammerTools_INCLUDE_DIR}/HammerSerial_generated.h
	${HammerTools_INCLUDE_DIR}/HammerSerial.fhh
	${HammerTools_INCLUDE_DIR}/HammerSerial.hh
	${HammerTools_INCLUDE_DIR}/HammerYaml.hh
	${HammerTools_INCLUDE_DIR}/IOTypes.hh
	${HammerTools_INCLUDE_DIR}/Iterators.hh
	${HammerTools_INCLUDE_DIR}/Logging.hh
	${HammerTools_INCLUDE_DIR}/Loki.hh
	${HammerTools_INCLUDE_DIR}/ParticleData.hh
	${HammerTools_INCLUDE_DIR}/ParticleUtils.hh
	${HammerTools_INCLUDE_DIR}/Pdg.fhh
	${HammerTools_INCLUDE_DIR}/Pdg.hh
	${HammerTools_INCLUDE_DIR}/Setting.hh
	${HammerTools_INCLUDE_DIR}/SettingDefs.hh
	${HammerTools_INCLUDE_DIR}/SettingsConsumer.hh
	${HammerTools_INCLUDE_DIR}/SettingsConsumerDefs.hh
	${HammerTools_INCLUDE_DIR}/Utils.hh
	${HammerTools_INCLUDE_DIR}/../IndexLabels.hh
	${HammerTools_INCLUDE_DIR}/../IndexTypes.hh
    ${HammerTools_INCLUDE_DIR}/../Exceptions.hh
	${HammerTools_INCLUDE_DIR}/../SettingsHandler.hh
	${HammerTools_INCLUDE_DIR}/../SettingsHandlerDefs.hh
	${HammerTools_INCLUDE_DIR}/SettingVisitors.hh
)

# pre-Cmake-3.23 lengthy commands
# Tier 1: minimal set needed to compile against Hammer (transitive closure of Hammer.hh)
set(Tools_Headers_Inst
	${HammerTools_INCLUDE_DIR}/HammerRoot.fhh
	${HammerTools_INCLUDE_DIR}/HammerSerial.fhh
	${HammerTools_INCLUDE_DIR}/IOTypes.hh
	${HammerTools_INCLUDE_DIR}/Logging.hh
	${HammerTools_INCLUDE_DIR}/Pdg.fhh
	${HammerTools_INCLUDE_DIR}/Pdg.hh
	${HammerTools_INCLUDE_DIR}/Setting.hh
	${HammerTools_INCLUDE_DIR}/SettingDefs.hh
	${HammerTools_INCLUDE_DIR}/SettingsConsumer.hh
	${HammerTools_INCLUDE_DIR}/SettingsConsumerDefs.hh
	${HammerTools_INCLUDE_DIR}/Utils.hh
)
separate_arguments(Tools_Headers_Inst)
set_target_properties(HammerTools PROPERTIES
	PUBLIC_HEADER
	"${Tools_Headers_Inst}"
)

target_include_directories(HammerTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
set_target_properties(HammerTools PROPERTIES
	POSITION_INDEPENDENT_CODE ${Hammer_PIC}
	VERSION ${Hammer_VERSION}
	SOVERSION ${Hammer_VERSION_MAJOR})
target_compile_features(HammerTools PUBLIC cxx_std_17)

if(VERBOSE_DEBUG)
	target_compile_options(HammerTools PRIVATE ${Hammer_VerboseOptions})
endif()

set(TOOLS_TARGETS "HammerConfig Boost::headers yaml-cpp::yaml-cpp")

if(HAVE_ROOT)
	set(TOOLS_TARGETS "${TOOLS_TARGETS} ROOT::Core")
endif()

separate_arguments(TOOLS_TARGETS)
target_link_libraries(HammerTools PUBLIC ${TOOLS_TARGETS})

if(SANITIZE)
	target_link_libraries(HammerTools PRIVATE HammerSanitizers)
endif()

install(TARGETS HammerTools
	EXPORT HammerTargets
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	COMPONENT Hammer_Runtime
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
	COMPONENT Hammer_Runtime
	NAMELINK_COMPONENT Hammer_Development
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
	COMPONENT Hammer_Development
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Hammer/Tools
	COMPONENT Hammer_Development
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# pre-Cmake-3.23 compatible to add core headers (instead of FILE_SET HEADERS)
install(FILES
	${HammerTools_INCLUDE_DIR}/../IndexTypes.hh
	${HammerTools_INCLUDE_DIR}/../Exceptions.hh
	${HammerTools_INCLUDE_DIR}/../SettingsHandler.hh
	${HammerTools_INCLUDE_DIR}/../SettingsHandlerDefs.hh
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Hammer
	COMPONENT Hammer_Development
)


