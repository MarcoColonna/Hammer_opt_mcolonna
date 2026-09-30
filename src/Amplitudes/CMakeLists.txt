message(STATUS "Preparing Hammer Amplitudes library.")

add_library(Amplitudes)
target_sources(Amplitudes
	PRIVATE
	AmplBToQLepNuBase.cc
	AmplBDLepNu.cc
	AmplBDstarLepNu.cc
	AmplBDstarDPiLepNu.cc
	AmplBDstarDGamLepNu.cc
	AmplBD0starLepNu.cc
	AmplBD1starLepNu.cc
	AmplBD1LepNu.cc
	AmplBD2starLepNu.cc
	AmplBcJpsiEllEllLepNu.cc
	AmplLbLcLepNu.cc
	AmplLbLcstar12LepNu.cc
	AmplLbLcstar32LepNu.cc
	AmplBRhoPiPiLepNu.cc
	AmplBOmega3PiLepNu.cc
	AmplD1starDstarDPiPi.cc
	AmplD1starDstarDGamPi.cc
	AmplD1DstarDPiPi.cc
	AmplD1DstarDGamPi.cc
	AmplD2starDstarDPiPi.cc
	AmplD2starDstarDGamPi.cc
	AmplD2starDPi.cc
	AmplTauEllNuNu.cc
	AmplTauPiNu.cc
	AmplTau3PiNu.cc
)

set(Amplitudes_INCLUDE_DIR "../../include/Hammer/Amplitudes")

target_sources(Amplitudes
	PRIVATE
	${Amplitudes_INCLUDE_DIR}/AmplBToQLepNuBase.hh
	${Amplitudes_INCLUDE_DIR}/AmplBDLepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplBDstarLepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplBDstarDPiLepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplBDstarDGamLepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplBD0starLepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplBD1starLepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplBD1LepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplBD2starLepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplBcJpsiEllEllLepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplLbLcLepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplLbLcstar12LepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplLbLcstar32LepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplBRhoPiPiLepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplBOmega3PiLepNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplD1starDstarDPiPi.hh
	${Amplitudes_INCLUDE_DIR}/AmplD1starDstarDGamPi.hh
	${Amplitudes_INCLUDE_DIR}/AmplD1DstarDPiPi.hh
	${Amplitudes_INCLUDE_DIR}/AmplD1DstarDGamPi.hh
	${Amplitudes_INCLUDE_DIR}/AmplD2starDstarDPiPi.hh
	${Amplitudes_INCLUDE_DIR}/AmplD2starDstarDGamPi.hh
	${Amplitudes_INCLUDE_DIR}/AmplD2starDPi.hh
	${Amplitudes_INCLUDE_DIR}/AmplTauEllNuNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplTauPiNu.hh
	${Amplitudes_INCLUDE_DIR}/AmplTau3PiNu.hh
)

target_include_directories(Amplitudes PRIVATE ${PROJECT_SOURCE_DIR}/include)
set_target_properties(Amplitudes PROPERTIES
	POSITION_INDEPENDENT_CODE ${Hammer_PIC}
	VERSION ${Hammer_VERSION}
	SOVERSION ${Hammer_VERSION_MAJOR})
target_compile_features(Amplitudes PUBLIC cxx_std_17)

if(VERBOSE_DEBUG)
	target_compile_options(Amplitudes PRIVATE ${Hammer_VerboseOptions})
endif()

target_link_libraries(Amplitudes PUBLIC HammerTools HammerBase HammerMath HammerTensor)

if(SANITIZE)
	target_link_libraries(Amplitudes PRIVATE HammerSanitizers)
endif()

install(TARGETS Amplitudes
	EXPORT HammerTargets
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	COMPONENT Hammer_Runtime
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
	COMPONENT Hammer_Runtime
	NAMELINK_COMPONENT Hammer_Development
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
	COMPONENT Hammer_Development
)
