message(STATUS "Building python hammerlib package")

if(Cython_FOUND)
    include(UseCython)
    set(CYTHON_FLAGS "-X language_level=3")

    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
    set_source_files_properties(pyHammer.pyx PROPERTIES
        CYTHON_IS_CXX TRUE)

    set(CYTHON_WARNINGS "-Wno-unused-macros -Wno-old-style-cast -Wno-sign-conversion")
    set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-float-equal -Wno-unused-function -Wno-sign-compare")

    if(CMAKE_COMPILER_IS_GNUCXX)
        set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-float-equal -Wno-shadow -Wno-cast-qual")
        set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-conversion")
    else(CMAKE_COMPILER_IS_GNUCXX)
        set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-unneeded-internal-declaration")
        set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-reserved-id-macro -Wno-missing-variable-declarations")
        set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-unreachable-code-return -Wno-implicit-fallthrough")
        set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-used-but-marked-unused -Wno-shorten-64-to-32")
        set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-zero-as-null-pointer-constant -Wno-c++17-extensions")
        set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-double-promotion -Wno-cast-qual -Wno-cast-align")
        set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-conditional-uninitialized -Wno-shadow -Wno-missing-noreturn ")
        set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-comma -Wno-format-nonliteral -Wno-implicit-int-conversion")
        set(CYTHON_WARNINGS "${CYTHON_WARNINGS} -Wno-disabled-macro-expansion")
    endif(CMAKE_COMPILER_IS_GNUCXX)

    separate_arguments(CYTHON_WARNINGS)

    cython_add_module(pyHammer pyHammer.pyx)
    target_link_libraries(pyHammer PRIVATE HammerBase HammerTools HammerCore Amplitudes Rates FormFactors)

    set_target_properties(pyHammer PROPERTIES POSITION_INDEPENDENT_CODE ${Hammer_PIC})

    if(VERBOSE_DEBUG)
        target_compile_options(pyHammer PRIVATE ${Hammer_VerboseOptions})
    endif()

    if(SANITIZE)
        target_link_libraries(pyHammer PRIVATE HammerSanitizers)
    endif()

    target_compile_options(pyHammer PRIVATE ${CYTHON_WARNINGS})

    set(MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/bindings")
    set(HAMMER_PYTHON_DEPS "[]")
    set(HAMMER_PYTHON_OPTIONAL_DEPS "[\"numpy\"]")
    set(HAMMER_BINDING_CLASSIFIERS "'Programming Language :: Cython',")

    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../pyproject.toml.in"
        "${CMAKE_CURRENT_BINARY_DIR}/pyproject.toml"
        @ONLY)

    add_custom_target(SetUpModules DEPENDS pyHammer
        COMMAND ${CMAKE_COMMAND} -E make_directory ${MODULE_DIR}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${MODULE_DIR}/hammer)

    set(ModuleFiles "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/pdg.py"
        "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/$<TARGET_FILE_NAME:pyHammer>"
    )
    separate_arguments(ModuleFiles)

    if(WITH_EXAMPLES)
        set(ModuleFiles "${ModuleFiles};${CMAKE_CURRENT_SOURCE_DIR}/../hepmc.py")
    endif()

    set(PackageFiles "${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in"
        "${CMAKE_SOURCE_DIR}/GUIDELINES"
        "${CMAKE_SOURCE_DIR}/COPYING"
        "${CMAKE_SOURCE_DIR}/README.md"
        "${CMAKE_CURRENT_SOURCE_DIR}/../setup.cfg"
        "${CMAKE_CURRENT_BINARY_DIR}/pyproject.toml"
    )
    separate_arguments(PackageFiles)

    foreach(ModuleFile ${ModuleFiles})
        add_custom_command(TARGET SetUpModules POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E
            copy ${ModuleFile} "${MODULE_DIR}/hammer")
    endforeach()

    foreach(PackageFile ${PackageFiles})
        add_custom_command(TARGET SetUpModules PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E
            copy ${PackageFile} "${MODULE_DIR}")
    endforeach()

    add_custom_target(PythonPackages ALL DEPENDS SetUpModules pyHammer)

else()
    add_custom_target(PythonPackages)
endif()
