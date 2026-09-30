include(FindPythonModule)

function(find_pip)
    message(STATUS "Looking for pip availability.")
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip -V
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE
    )
    if(result EQUAL 0)
        set(Python3_pip_FOUND TRUE PARENT_SCOPE)
        message(STATUS "Found pip: ${output}.")
    else()
        message(WARNING "Pip not found: unable to install missing packages.")
    endif()
endfunction()

function(install_python_module Package)
    # ARGV0 = pip package name
    # ARGV1 = optional version constraint (empty string or omit = no pin)
    # ARGV2 = optional Python import name (defaults to Package when pip name != import name)
    if(ARGC GREATER 2 AND NOT "${ARGV2}" STREQUAL "")
        set(_import_name "${ARGV2}")
    else()
        set(_import_name "${Package}")
    endif()

    if(ARGC GREATER 1 AND NOT "${ARGV1}" STREQUAL "")
        set(_install_string "${Package}==${ARGV1}")
    else()
        set(_install_string "${Package}")
    endif()

    message(STATUS "Trying to install python package '${Package}' with pip.")
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip install ${_install_string}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error
    )
    if(result EQUAL 0)
        find_python_module(${_import_name} QUIET)
        if(Python3_${_import_name}_FOUND)
            set(Python3_${_import_name}_INSTALLED TRUE PARENT_SCOPE)
            message(STATUS "Done.")
        else()
            # cimport probe — only valid for Cython-API packages, not the compiler itself
            if(NOT "${_import_name}" STREQUAL "Cython")
                find_cython_module(${_import_name})
                if(Python3_${_import_name}_FOUND)
                    set(Python3_${_import_name}_INSTALLED TRUE PARENT_SCOPE)
                    message(STATUS "Done.")
                else()
                    message(WARNING "Could not find ${_import_name} after installation with pip.")
                endif()
            endif()
        endif()
    else()
        message(WARNING "Could not install ${Package}: ${result}")
    endif()
endfunction()
