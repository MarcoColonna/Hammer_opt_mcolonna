###############################################################################
# Find a python module
#
# This sets the following variables:
# PY<MODULE_NAME>_FOUND - True if the package named <MODULE_NAME> was found.
# PY<MODULE_NAME> containts the path
#

# Find if a Python module is installed
# based on function found at http://www.cmake.org/pipermail/cmake/2011-January/041666.html
# To use do: find_python_module(PyQt4 REQUIRED)

function(find_python_module module_name)
	cmake_parse_arguments(PARSE_ARGV 1 arg "REQUIRED;QUIET" "" "")
	# A module's location is usually a directory, but for binary modules
	# it's a .so file.
	string(TOLOWER "${module_name}" module_lower)
	execute_process(
		COMMAND "${Python3_EXECUTABLE}" "-c"
			"import importlib.util; spec = importlib.util.find_spec('${module_lower}'); print(spec.origin if spec and spec.origin else '')"
		RESULT_VARIABLE _${module_name}_status
		OUTPUT_VARIABLE _${module_name}_location
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	if(NOT _${module_name}_status AND _${module_name}_location)
		set(Python3_${module_name}_PATH "${_${module_name}_location}" CACHE STRING
			"Location of Python module ${module_name}")
		set(Python3_${module_name}_FOUND TRUE)
		if(NOT arg_QUIET)
			message(STATUS "Found Python module ${module_name}: ${_${module_name}_location}")
		endif()
	else()
		set(Python3_${module_name}_FOUND FALSE)
		if(arg_REQUIRED)
			message(FATAL_ERROR "Required Python module ${module_name} not found")
		elseif(NOT arg_QUIET)
			message(STATUS "Python module ${module_name} not found")
		endif()
	endif()
	set(Python3_${module_name}_FOUND "${Python3_${module_name}_FOUND}" PARENT_SCOPE)
	mark_as_advanced(_${module_name}_status _${module_name}_location)
endfunction()

function(find_cython_module module_name)
	file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/CMakeFiles/cython_probe")
	file(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/cython_probe/tmp_${module_name}.pyx" "cimport ${module_name}\n")
	# A module's location is usually a directory, but for binary modules
	# it's a .so file.
	execute_process(COMMAND "${Cython_EXECUTABLE}" "-3" "--cplus"
		"${CMAKE_BINARY_DIR}/CMakeFiles/cython_probe/tmp_${module_name}.pyx"
		RESULT_VARIABLE _${module_name}_status
		OUTPUT_VARIABLE _${module_name}_output
		ERROR_QUIET
		OUTPUT_STRIP_TRAILING_WHITESPACE
		ERROR_STRIP_TRAILING_WHITESPACE)
	file(REMOVE
		"${CMAKE_BINARY_DIR}/CMakeFiles/cython_probe/tmp_${module_name}.pyx"
		"${CMAKE_BINARY_DIR}/CMakeFiles/cython_probe/tmp_${module_name}.cpp"
		"${CMAKE_BINARY_DIR}/CMakeFiles/cython_probe/tmp_${module_name}.h"
	)
	if(NOT _${module_name}_status)
		set(Python3_${module_name}_FOUND TRUE)
	else()
		set(Python3_${module_name}_FOUND FALSE)
	endif()
	set(Python3_${module_name}_FOUND "${Python3_${module_name}_FOUND}" PARENT_SCOPE)
endfunction()
