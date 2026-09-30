# set minimum required cmake version
cmake_minimum_required(VERSION 3.21...4.3)

# make sure we are not going to build in the main source directory
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -------------------------------------------------------------------------------
# Project metadata
# -------------------------------------------------------------------------------
set(HAMMER_VERSION_MAJOR "2")
set(HAMMER_VERSION_MINOR "0")
set(HAMMER_VERSION_PATCH "0")
set(HAMMER_VERSION "${HAMMER_VERSION_MAJOR}.${HAMMER_VERSION_MINOR}.${HAMMER_VERSION_PATCH}")

project(Hammer VERSION "${HAMMER_VERSION}"
	DESCRIPTION "Helicity Amplitude Module for Matrix Element Reweighting"
	HOMEPAGE_URL "https://hammer.physics.lbl.gov"
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	set(CMAKE_CXX_EXTENSIONS OFF)
	set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules/")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/lib)

# -------------------------------------------------------------------------------
# Options
# -------------------------------------------------------------------------------
option(ENABLE_TESTS "Enable Hammer unit tests" OFF)
option(WITH_PYTHON "Install the Hammer python bindings" ON)
option(PYTHON_USE_CPPYY "Use cppyy to generate the python bindings. (OFF = use cython)" ON)
option(WITH_EXAMPLES "Install and build Hammer examples" OFF)
option(WITH_ROOT "Turn on ROOT histogram support" OFF)
option(INSTALL_EXTERNAL_DEPENDENCIES "Enable automatic download/install of external dependencies" OFF)
option(BUILD_DOCUMENTATION "Build the doxygen documentation" OFF)
option(VERBOSE_DEBUG "Turn on verbose debugging" OFF)
option(SANITIZE "Turn on sanitizers" OFF)
option(INSTALL_HEADERS "Install all Hammer development headers" OFF)
option(COVERAGE "Enable coverage reporting" OFF)
option(FORCE_YAMLCPP_INSTALL "Force YAMLCPP installation even if present" OFF)
option(FORCE_HEPMC_INSTALL "Force HepMC installation even if present" OFF)

# -------------------------------------------------------------------------------
# Build type adjustments
# -------------------------------------------------------------------------------
if(SANITIZE AND NOT VERBOSE_DEBUG)
	set(VERBOSE_DEBUG ON CACHE BOOL "" FORCE)
endif()

if(VERBOSE_DEBUG AND NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
	message(STATUS "Verbose debugging turned on: setting build type to 'Debug'")
	set(CMAKE_BUILD_TYPE "Debug")
endif()

if(COVERAGE)
	message("Building with coverage instrumentation...")
	if(NOT ${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
		message(STATUS "Coverage turned on: setting build type to 'RelWithDebInfo'")
		set(CMAKE_BUILD_TYPE "RelWithDebInfo")
	endif()
	add_compile_options(-fprofile-instr-generate -fcoverage-mapping)
	add_link_options(-fprofile-instr-generate)
endif()

# -------------------------------------------------------------------------------
# Dependencies
# -------------------------------------------------------------------------------

# --- ROOT ---
set(HAVE_ROOT FALSE)
if(WITH_ROOT)
	set(_ROOT_MIN_VERSION "6.25")
	set(ROOT_COMPONENTS_LIST Hist Physics)
	if(WITH_PYTHON)
		list(APPEND ROOT_COMPONENTS_LIST ROOTTPython)
	endif()
	if(WITH_EXAMPLES)
		list(APPEND ROOT_COMPONENTS_LIST Tree Gpad)
	endif()

	find_package(ROOT ${_ROOT_MIN_VERSION} CONFIG COMPONENTS ${ROOT_COMPONENTS_LIST})

	if(NOT ROOT_FOUND)
		include(FindROOT_old)
		if(NOT ROOT_FOUND)
			message(FATAL_ERROR
				"WITH_ROOT=ON but ROOT >= ${_ROOT_MIN_VERSION} was not found. "
				"Install ROOT (>= ${_ROOT_MIN_VERSION}) or set WITH_ROOT=OFF.")
		endif()
		if(DEFINED ROOT_VERSION AND ROOT_VERSION VERSION_LESS _ROOT_MIN_VERSION)
			message(FATAL_ERROR
				"Found ROOT ${ROOT_VERSION} but >= ${_ROOT_MIN_VERSION} is required.")
		endif()
		set(ROOT_Found_HOW "Old")
	else()
		set(ROOT_Found_HOW "New")
	endif()

	set(HAVE_ROOT TRUE)
	include(FixROOT)
	unset(_ROOT_MIN_VERSION)
endif()

# --- Compiler checks ---
include(CompilerChecks)

# --- Boost ---
set(BOOST_MINREQ_VERSION "1.75.0")
find_package(Boost ${BOOST_MINREQ_VERSION} REQUIRED COMPONENTS thread)
message(STATUS "Found Boost version ${Boost_VERSION}.")

# --- yaml-cpp ---
set(YamlCpp_MIN_VERSION "0.6.0")
find_package(yaml-cpp ${YamlCpp_MIN_VERSION} CONFIG QUIET)
if((NOT yaml-cpp_FOUND AND INSTALL_EXTERNAL_DEPENDENCIES) OR FORCE_YAMLCPP_INSTALL)
	include(AddYamlCpp)
	set(YamlCpp_FOUND_HOW "Added")
elseif(NOT yaml-cpp_FOUND)
	message(FATAL_ERROR
		"yaml-cpp >= ${YamlCpp_MIN_VERSION} not found. "
		"Install it or configure with -DINSTALL_EXTERNAL_DEPENDENCIES=ON.")
else()
	set(YamlCpp_FOUND_HOW "New")
endif()

# --- HepMC3 (required only for examples) ---
set(HAVE_HEPMC3 FALSE)
if(WITH_EXAMPLES)
	find_package(HepMC3 3.2.0 CONFIG)
	if((NOT HepMC3_FOUND AND INSTALL_EXTERNAL_DEPENDENCIES) OR FORCE_HEPMC_INSTALL)
		include(AddHepMC)
		set(HepMC_Found_HOW "Added")
		set(HAVE_HEPMC3 TRUE)
	elseif(NOT HepMC3_FOUND)
		message(WARNING
			"HepMC3 not found. Examples using HepMC3 will not be built. "
			"Install HepMC3 or configure with -DINSTALL_EXTERNAL_DEPENDENCIES=ON.")
	else()
		set(HepMC_Found_HOW "New")
		set(HAVE_HEPMC3 TRUE)
	endif()
endif()

# --- Python ---
if(WITH_PYTHON)
	set(Python3_FIND_VIRTUALENV FIRST)

	if(WITH_ROOT AND HAVE_ROOT)
		# Python version must exactly match the one used to build PyROOT
		string(REGEX MATCH "([0-9]+)\\.([0-9]+)" _Python_Exact_VERSION "${ROOT_PYTHON_VERSION}")
		find_package(Python3 ${_Python_Exact_VERSION} EXACT COMPONENTS Interpreter Development)
		unset(_Python_Exact_VERSION)
	else()
		find_package(Python3 "3.10" COMPONENTS Interpreter Development)
	endif()

	if(NOT (Python3_Interpreter_FOUND AND Python3_Development_FOUND))
		set(WITH_PYTHON OFF CACHE BOOL "" FORCE)
		message(WARNING
			"Python3 >= 3.10 not found. Hammer will only be accessible as a C++ library.")
	endif()
endif()

if(WITH_PYTHON)
	include(AddPythonModule)

	if(INSTALL_EXTERNAL_DEPENDENCIES)
		find_pip()
	endif()

	# build and setuptools are required regardless of the binding mode
	foreach(_pymod IN ITEMS setuptools build)
		find_python_module(${_pymod} QUIET)
		if(NOT Python3_${_pymod}_FOUND)
			if(Python3_pip_FOUND)
				install_python_module(${_pymod})
			endif()
			if(NOT Python3_${_pymod}_INSTALLED)
				message(WARNING "${_pymod} not found/installed. Disabling Python bindings.")
				set(WITH_PYTHON OFF CACHE BOOL "" FORCE)
				break()
			endif()
		endif()
	endforeach()
	unset(_pymod)
endif()

if(WITH_PYTHON)
	if(PYTHON_USE_CPPYY)
		if(WITH_ROOT AND HAVE_ROOT)
			message(STATUS "Using ROOT cppyy interface for python bindings.")
			if(NOT DEFINED ROOT_ROOTTPython_LIBRARY)
				message(WARNING "ROOT found but without PyROOT configuration. Disabling Python.")
				set(WITH_PYTHON OFF CACHE BOOL "" FORCE)
			endif()
		else()
			# Standalone cppyy (no ROOT)
			set(CPPYY_MODULE_PATH "" CACHE PATH "Path to cppyy CMake modules (tried via cling-config, fall back on manually)")
			if(NOT CPPYY_MODULE_PATH)
				execute_process(
					COMMAND cling-config --cmake
					OUTPUT_VARIABLE CPPYY_MODULE_PATH
					OUTPUT_STRIP_TRAILING_WHITESPACE)
			endif()
			if(NOT CPPYY_MODULE_PATH AND DEFINED ROOT_DIR)
				set(CPPYY_MODULE_PATH "${ROOT_DIR}")
				message(STATUS
					"cling-config not found; derived CPPYY_MODULE_PATH from ROOT_DIR: ${CPPYY_MODULE_PATH}")
			else()
				message(STATUS "CPPYY_MODULE_PATH: ${CPPYY_MODULE_PATH}")
			endif()
			list(INSERT CMAKE_MODULE_PATH 0 "${CPPYY_MODULE_PATH}")
			find_package(Cppyy)

			if(NOT Cppyy_FOUND)
				if(Python3_pip_FOUND)
					install_python_module(cppyy)
				endif()
				if(Python3_cppyy_INSTALLED)
					# Re-detect CMake package after pip install
					execute_process(
						COMMAND cling-config --cmake
						OUTPUT_VARIABLE CPPYY_MODULE_PATH
						OUTPUT_STRIP_TRAILING_WHITESPACE)
					if(CPPYY_MODULE_PATH)
						list(INSERT CMAKE_MODULE_PATH 0 "${CPPYY_MODULE_PATH}")
					endif()
					find_package(Cppyy)
				endif()
				if(NOT Cppyy_FOUND)
					message(WARNING "cppyy not found/installed. Disabling Python.")
					set(WITH_PYTHON OFF CACHE BOOL "" FORCE)
				endif()
			endif()
		endif()

	else()
		# Cython bindings
		include(FindCython)

		if(NOT Cython_FOUND)
			if(Python3_pip_FOUND)
				install_python_module(Cython)
			endif()
			if(NOT Python3_Cython_INSTALLED)
				message(WARNING "Cython compiler not found. Disabling Python.")
				set(WITH_PYTHON OFF CACHE BOOL "" FORCE)
			else()
				include(FindCython)
			endif()
		endif()

		if(WITH_PYTHON)
			find_cython_module(cymove)
			if(NOT Python3_cymove_FOUND)
				if(Python3_pip_FOUND)
					install_python_module(cymove "1.0.0")
				endif()
				if(NOT Python3_cymove_INSTALLED)
					message(WARNING "cymove not found/installed. Disabling Python.")
					set(WITH_PYTHON OFF CACHE BOOL "" FORCE)
				endif()
			endif()
		endif()

		if(WITH_PYTHON)
			find_package(Python3 COMPONENTS NumPy)
			if(NOT Python3_NumPy_FOUND)
				if(Python3_pip_FOUND)
					install_python_module(numpy)
				endif()
				if(NOT Python3_numpy_INSTALLED)
					message(WARNING "numpy not found/installed. Disabling Python.")
					set(WITH_PYTHON OFF CACHE BOOL "" FORCE)
				endif()
			endif()
		endif()
	endif()
endif()

# --- Doxygen ---
set(PYTHON_DOCS FALSE)
if(BUILD_DOCUMENTATION)
	find_package(Doxygen OPTIONAL_COMPONENTS dot)

	if(NOT DOXYGEN_FOUND)
		message(WARNING "BUILD_DOCUMENTATION=ON but doxygen not found. Disabling documentation.")
		set(BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
	elseif(NOT DOXYGEN_DOT_FOUND)
		message(WARNING "doxygen found but 'dot' (graphviz) not found. Disabling documentation.")
		set(BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
	endif()
endif()

if(BUILD_DOCUMENTATION)
	if(WITH_PYTHON AND Python3_Interpreter_FOUND)
		find_python_module(doxypypy QUIET)
		if(NOT Python3_doxypypy_FOUND)
			if(Python3_pip_FOUND)
				install_python_module(doxypypy)
			endif()
			if(Python3_doxypypy_INSTALLED)
				set(PYTHON_DOCS TRUE)
			else()
				message(WARNING
					"doxypypy not found/installed. Python API documentation will not be built.")
			endif()
		else()
			set(PYTHON_DOCS TRUE)
		endif()
	endif()

	configure_file(
		${PROJECT_SOURCE_DIR}/Utils/hammer.doxy.in
		${PROJECT_BINARY_DIR}/hammer.doxy
		ESCAPE_QUOTES @ONLY)
	add_custom_target(doc ALL
		COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/hammer.doxy
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
		DEPENDS ${PROJECT_BINARY_DIR}/hammer.doxy)
	install(DIRECTORY "${PROJECT_BINARY_DIR}/html"
		COMPONENT documentation
		DESTINATION ${CMAKE_INSTALL_DOCDIR}
		OPTIONAL)
endif()

# -------------------------------------------------------------------------------
# Build targets
# -------------------------------------------------------------------------------
add_subdirectory(src)

set(PYHEPMC3_SITE_PACKAGES "")
if(WITH_EXAMPLES)
	if(HAVE_HEPMC3 AND WITH_PYTHON AND Python3_Interpreter_FOUND)
		execute_process(
			COMMAND ${Python3_EXECUTABLE} -c
				"import pyHepMC3, os; print(os.path.dirname(os.path.dirname(pyHepMC3.__file__)))"
			OUTPUT_VARIABLE PYHEPMC3_SITE_PACKAGES
			OUTPUT_STRIP_TRAILING_WHITESPACE
			ERROR_QUIET
			RESULT_VARIABLE _pyHepMC3_result)
		if(NOT _pyHepMC3_result EQUAL 0 OR NOT PYHEPMC3_SITE_PACKAGES)
			set(PYHEPMC3_SITE_PACKAGES "")
			execute_process(
				COMMAND ${Python3_EXECUTABLE} -c "import pyhepmc"
				ERROR_QUIET
				RESULT_VARIABLE _pyhepmc_result)
			if(_pyhepmc_result EQUAL 0)
				message(STATUS
					"pyHepMC3 not found for Python ${Python3_VERSION}; "
					"pyhepmc is available as fallback.")
			else()
				message(STATUS
					"Neither pyHepMC3 nor pyhepmc found for Python ${Python3_VERSION}. "
					"Python HepMC3 examples will use the internal hammer.hepmc parser.")
			endif()
		else()
			message(STATUS "Found pyHepMC3 at: ${PYHEPMC3_SITE_PACKAGES}")
		endif()
	endif()

	add_subdirectory(Examples)
endif()

if(WITH_PYTHON)
	set(PYTHON_INSTALL_DIR
		"${CMAKE_INSTALL_FULL_LIBDIR}/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages")
	add_subdirectory(pyext)
endif()

add_subdirectory(include)

if(ENABLE_TESTS)
	enable_testing()
	add_subdirectory(Tests)
endif()

add_subdirectory(scripts)
add_subdirectory(CMakeModules)

# -------------------------------------------------------------------------------
# Packaging (CPack)
# -------------------------------------------------------------------------------
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${Hammer_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${Hammer_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${Hammer_VERSION_PATCH}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
	"The Hammer Package: Helicity Amplitude Module for Matrix Element Reweighting")
set(CPACK_PACKAGE_VENDOR "The Hammer Collaboration")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README)
set(CPACK_GENERATOR TGZ)
set(CPACK_SOURCE_GENERATOR TGZ)
set(CPACK_SOURCE_IGNORE_FILES
	"/[.]AppleDouble/"
	"/[.]DS_Store"
	"^${CMAKE_CURRENT_SOURCE_DIR}/[.][amvgc]"
	"/build/"
	"/devel/"
	"/python/"
	"/master/"
	"/sandbox/"
	"/other/"
	"/docs/"
	"/pyext/wrapper_cppyy/old/"
	"/graphify-out/"
	"^${CMAKE_CURRENT_SOURCE_DIR}/internal/"
	"/Talks/"
	"CLAUDE.md"
	"^${CMAKE_CURRENT_SOURCE_DIR}/Utils/[cfb]"
)
include(CPack)
