message(STATUS "Preparing Hammer Rates library.")

add_library(Rates)
target_sources(Rates
	PRIVATE
	RateBDLepNu.cc
	RateBDstarLepNu.cc
	RateBD0starLepNu.cc
	RateBD1starLepNu.cc
	RateBD1LepNu.cc
	RateBD2starLepNu.cc
	RateLbLcLepNu.cc
	RateLbLcstar12LepNu.cc
	RateLbLcstar32LepNu.cc
	RateBRhoOmegaLepNu.cc
	RateTau3PiNu.cc
	RateD1starDstarPi.cc
	RateD1DstarPi.cc
	RateD2starDstarPi.cc
	RateD2starDPi.cc
)

set(Rates_INCLUDE_DIR "../../include/Hammer/Rates")

target_sources(Rates
	PRIVATE
	${Rates_INCLUDE_DIR}/RateBDLepNu.hh
	${Rates_INCLUDE_DIR}/RateBDstarLepNu.hh
	${Rates_INCLUDE_DIR}/RateBD0starLepNu.hh
	${Rates_INCLUDE_DIR}/RateBD1starLepNu.hh
	${Rates_INCLUDE_DIR}/RateBD1LepNu.hh
	${Rates_INCLUDE_DIR}/RateBD2starLepNu.hh
	${Rates_INCLUDE_DIR}/RateLbLcLepNu.hh
	${Rates_INCLUDE_DIR}/RateLbLcstar12LepNu.hh
	${Rates_INCLUDE_DIR}/RateLbLcstar32LepNu.hh
	${Rates_INCLUDE_DIR}/RateBRhoOmegaLepNu.hh
	${Rates_INCLUDE_DIR}/RateTau3PiNu.hh
	${Rates_INCLUDE_DIR}/RateD1starDstarPi.hh
	${Rates_INCLUDE_DIR}/RateD1DstarPi.hh
	${Rates_INCLUDE_DIR}/RateD2starDstarPi.hh
	${Rates_INCLUDE_DIR}/RateD2starDPi.hh
)

target_include_directories(Rates PRIVATE ${PROJECT_SOURCE_DIR}/include)
set_target_properties(Rates PROPERTIES
	POSITION_INDEPENDENT_CODE ${Hammer_PIC}
	VERSION ${Hammer_VERSION}
	SOVERSION ${Hammer_VERSION_MAJOR})
target_compile_features(Rates PUBLIC cxx_std_17)

if(VERBOSE_DEBUG)
	target_compile_options(Rates PRIVATE ${Hammer_VerboseOptions})
endif()

target_link_libraries(Rates PUBLIC HammerTools HammerBase HammerMath HammerTensor)

if(SANITIZE)
	target_link_libraries(Rates PRIVATE HammerSanitizers)
endif()

install(TARGETS Rates
	EXPORT HammerTargets
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	COMPONENT Hammer_Runtime
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
	COMPONENT Hammer_Runtime
	NAMELINK_COMPONENT Hammer_Development
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
	COMPONENT Hammer_Development
)