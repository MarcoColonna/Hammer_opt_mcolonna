
message(STATUS "Preparing unit test support. Build the tests with 'make tests' and run them with 'ctest'.")

find_package(Threads REQUIRED)

include(ExternalProject)

set(GTEST_CXXFLAGS "")

if(SANITIZE)
    string(REPLACE ";" " " GTEST_CXXFLAGS "${Hammer_SanitizeOptions}")
endif()

# Download and install GoogleTest
ExternalProject_Add(
    gtest
    URL file://${PROJECT_SOURCE_DIR}/Tests/GTest
    PREFIX ${PROJECT_BINARY_DIR}/gtest

    # Disable install step
    INSTALL_COMMAND ""
    CMAKE_ARGS
    -DCMAKE_CXX_FLAGS=${GTEST_CXXFLAGS} -DGOOGLETEST_VERSION=1.7.0
)
ExternalProject_Get_Property(gtest source_dir binary_dir)

add_library(ThirdParty::GTest IMPORTED STATIC GLOBAL)
add_dependencies(ThirdParty::GTest gtest)

if(CMAKE_GENERATOR STREQUAL Xcode)
    set_target_properties(ThirdParty::GTest PROPERTIES
        "IMPORTED_LOCATION_DEBUG" "${binary_dir}/lib/Debug/libgtest.a"
        "IMPORTED_LOCATION_RELEASE" "${binary_dir}/lib/Release/libgtest.a"
        "IMPORTED_LOCATION_RELWITHDEBINFO" "${binary_dir}/lib/Debug/libgtest.a"
        "IMPORTED_LOCATION_MINSIZEREL" "${binary_dir}/lib/Release/libgtest.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
    )
else()
    set_target_properties(ThirdParty::GTest PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/lib/libgtest.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
    )
endif()

add_library(ThirdParty::GTestMain IMPORTED STATIC GLOBAL)
add_dependencies(ThirdParty::GTestMain gtest)

if(CMAKE_GENERATOR STREQUAL Xcode)
    set_target_properties(ThirdParty::GTestMain PROPERTIES
        "IMPORTED_LOCATION_DEBUG" "${binary_dir}/lib/Debug/libgtest_main.a"
        "IMPORTED_LOCATION_RELEASE" "${binary_dir}/lib/Release/libgtest_main.a"
        "IMPORTED_LOCATION_RELWITHDEBINFO" "${binary_dir}/lib/Debug/libgtest_main.a"
        "IMPORTED_LOCATION_MINSIZEREL" "${binary_dir}/lib/Release/libgtest_main.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
    )
else()
    set_target_properties(ThirdParty::GTestMain PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/lib/libgtest_main.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
    )
endif()

function(new_test name libs)
    add_executable(Test${name}
        Test${name}.cc)

    target_include_directories(Test${name} PRIVATE ${PROJECT_SOURCE_DIR}/include)

    target_include_directories(Test${name} SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/Tests/GTest/include)

    if(VERBOSE_DEBUG)
        target_compile_options(Test${name} PRIVATE ${Hammer_VerboseOptions})
    endif()

    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
        target_compile_options(Test${name} PRIVATE "-Wno-global-constructors;-Wno-covered-switch-default;-Wno-used-but-marked-unused")
    endif()

    set_target_properties(Test${name} PROPERTIES
        POSITION_INDEPENDENT_CODE ${Hammer_PIC})

    target_link_libraries(Test${name} PRIVATE ${libs} ThirdParty::GTest ThirdParty::GTestMain)
    if(BUILD_SHARED_LIBS)
        set_target_properties(Test${name} PROPERTIES
            BUILD_RPATH "${PROJECT_BINARY_DIR}/lib")
    endif()

    if(${ARGC} GREATER 2)
        SET(ARGN_COPY ${ARGN})

        foreach(MY_FILE ${ARGN_COPY})
            get_filename_component(MY_EXT "${MY_FILE}" EXT)
            get_filename_component(MY_NAME "${MY_FILE}" NAME_WE)

            if("${MY_EXT}" STREQUAL ".yaml.in")
                configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${MY_FILE} ${CMAKE_CURRENT_BINARY_DIR}/${MY_NAME}.yaml ESCAPE_QUOTES @ONLY)
            else()
                file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${MY_FILE}
                    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
            endif()
        endforeach()
    endif()

    if(SANITIZE)
        target_link_libraries(Test${name} PRIVATE HammerSanitizers)
    endif()

    add_test(NAME Test${name}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND Test${name} --gtest_color=yes)

    get_filename_component(namedir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    set_tests_properties(Test${name} PROPERTIES
        LABELS "${namedir}")

endfunction()

# add tests in subdirectories
add_subdirectory(Core)
add_subdirectory(Math)
add_subdirectory(Amplitudes)
add_subdirectory(FormFactors)
add_subdirectory(Rates)
add_subdirectory(Tools)
