# Find the Cython compiler.
#
# This code sets the following variables:
#
#  Cython_EXECUTABLE
#  Cython_VERSION
#  Cython_FOUND
#
# See also UseCython.cmake

#=============================================================================
# Copyright 2011 Kitware, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#=============================================================================

# Try the CMake config package first (Cython 3.x ships one).
# If it succeeds, Cython_EXECUTABLE and Cython_VERSION are already set.
find_package(Cython CONFIG QUIET)

if(NOT Cython_FOUND)
  # Fallback: locate the cython executable next to the Python interpreter.
  get_filename_component(_cython_hints "${Python3_EXECUTABLE}" DIRECTORY)
  find_program(Cython_EXECUTABLE
    NAMES cython cython3
    HINTS "${_cython_hints}"
  )
  unset(_cython_hints)
endif()

# Extract the version when not already provided by the CONFIG package.
if(Cython_EXECUTABLE AND NOT Cython_VERSION)
  execute_process(
    COMMAND "${Cython_EXECUTABLE}" --version
    OUTPUT_VARIABLE _cython_version_output
    ERROR_VARIABLE _cython_version_output
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
  )
  string(REGEX MATCH "[0-9]+\\.[0-9]+(\\.[0-9]+)?" Cython_VERSION "${_cython_version_output}")
  unset(_cython_version_output)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cython
  FOUND_VAR Cython_FOUND
  REQUIRED_VARS Cython_EXECUTABLE
  VERSION_VAR Cython_VERSION
)

mark_as_advanced(Cython_EXECUTABLE)
