if(WITH_ROOT)
    message(STATUS "Preparing bindings with ROOT.")
else()
    message(STATUS "Preparing bindings with cppyy standalone.")
endif()

add_library(pyHammer SHARED EXCLUDE_FROM_ALL)

# need to include again because ROOT doesn't get include dirs from linked libraries :(
target_include_directories(pyHammer PUBLIC ${PROJECT_SOURCE_DIR}/include
                                           ${PROJECT_BINARY_DIR}/include)
set_target_properties(pyHammer PROPERTIES POSITION_INDEPENDENT_CODE ${Hammer_PIC})
SET_TARGET_PROPERTIES(pyHammer PROPERTIES VERSION ${Hammer_VERSION})
SET_TARGET_PROPERTIES(pyHammer PROPERTIES SOVERSION ${Hammer_VERSION})


if(NOT WITH_ROOT)
    add_library(Cppyy::Core SHARED IMPORTED)

    set_target_properties(Cppyy::Core PROPERTIES
        INTERFACE_COMPILE_FEATURES "cxx_std_17"
        INTERFACE_INCLUDE_DIRECTORIES "${Cppyy_INCLUDE_DIRS}"
        IMPORTED_LOCATION "${Cppyy_DIR}/lib/libCoreLegacy.so"
        IMPORTED_SONAME "@rpath/libCoreLegacy.so"
    )

    target_link_libraries(pyHammer PUBLIC Cppyy::Core)
endif()

# install(TARGETS pyHammer LIBRARY DESTINATION lib ARCHIVE DESTINATION lib RUNTIME DESTINATION bin) 

if(VERBOSE_DEBUG)
    target_compile_options(pyHammer PRIVATE "${Hammer_VerboseOptions};-Wno-sign-conversion;-Wno-shorten-64-to-32;-Wno-implicit-int-conversion")
endif()

target_link_libraries(pyHammer PRIVATE HammerCore HammerConfig HammerMath HammerBase HammerTools Amplitudes FormFactors Rates)

if(SANITIZE)
    target_link_libraries(pyHammer PRIVATE HammerSanitizers)
endif()

if(WITH_ROOT)
    target_link_libraries(pyHammer PUBLIC ROOT::Physics ROOT::Hist)
endif()

# Hammer requires C++17; force it on pyHammer and (via COMPILE_OPTIONS inheritance)
# on the rootcling-generated dict OBJECT library.
target_compile_options(pyHammer PRIVATE -std=c++17)

set(LINKDEF_INCLUDES "Hammer/Math/FourMomentum.hh"
    "Hammer/Tools/Pdg.hh"
    "Hammer/Particle.hh"
    "Hammer/Process.hh"
    "Hammer/Hammer.hh"
    "Hammer/Tools/IOTypes.hh"
    "Hammer/IndexTypes.hh"
)
separate_arguments(LINKDEF_INCLUDES)


root_generate_dictionary(pyHammer_dict "${LINKDEF_INCLUDES}" MODULE pyHammer NOINSTALL LINKDEF HammerLinkDef.h)

set(MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/bindings")

if(WITH_ROOT)
    set(HAMMER_USE_ROOT_PYTHON "True")
    set(HAMMER_PYTHON_DEPS "[]")
    set(HAMMER_PYTHON_OPTIONAL_DEPS "[]")
    set(HAMMER_BINDING_CLASSIFIERS "")
else()
    set(HAMMER_USE_ROOT_PYTHON "False")
    set(HAMMER_PYTHON_DEPS "[\"cppyy\"]")
    set(HAMMER_PYTHON_OPTIONAL_DEPS "[]")
    set(HAMMER_BINDING_CLASSIFIERS "")
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../pyproject.toml.in"
    "${CMAKE_CURRENT_BINARY_DIR}/pyproject.toml"
    @ONLY)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/_hammer_config.py.in"
    "${CMAKE_CURRENT_BINARY_DIR}/_hammer_config.py"
    @ONLY)

add_custom_target(SetUpModules DEPENDS pyHammer
    COMMAND ${CMAKE_COMMAND} -E make_directory ${MODULE_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${MODULE_DIR}/hammer)

set(ModuleFiles "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
    "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}pyHammer_rdict.pcm"
    "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${CMAKE_SHARED_LIBRARY_PREFIX}pyHammer.rootmap"
    "$<TARGET_LINKER_FILE:pyHammer>"
)
separate_arguments(ModuleFiles)

if(WITH_EXAMPLES)
    set(ModuleFiles "${ModuleFiles};${CMAKE_CURRENT_SOURCE_DIR}/../hepmc.py")
endif()

if(WITH_ROOT)
    set(ModuleFiles "${ModuleFiles};${CMAKE_CURRENT_SOURCE_DIR}/root_pythonizer.py")
else()
    set(ModuleFiles "${ModuleFiles};${CMAKE_CURRENT_SOURCE_DIR}/cppyy_pythonizer.py")
endif()

set(ModuleFiles "${ModuleFiles};${CMAKE_CURRENT_BINARY_DIR}/_hammer_config.py")

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

# install(TARGETS pyHammer DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}/Hammer)
# install(FILES "${CMAKE_BINARY_DIR}/lib/libpyHammer_rdict.pcm" DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}/Hammer)
# install(FILES "${CMAKE_BINARY_DIR}/lib/libpyHammer.rootmap" DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}/Hammer)

# target_sources(pyHammer PUBLIC "G__pyHammer.cxx")