message(STATUS "Preparing Hammer Examples.")

add_custom_target(ExampleFiles)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
set(remote_Files "")
set(extra_Files "")

function(download_example_file file_name)
    if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/data/${file_name})
        message(STATUS "Downloading ${file_name}.bz2")
        file(DOWNLOAD https://hammer.physics.lbl.gov/mc_data/${file_name}.bz2 ${CMAKE_CURRENT_BINARY_DIR}/data/${file_name}.bz2 STATUS DOWNLOAD_STATUS)
        list(GET DOWNLOAD_STATUS 0 STATUS_CODE)
        list(GET DOWNLOAD_STATUS 1 ERROR_MESSAGE)

        # Check if download was successful.
        if(${STATUS_CODE} EQUAL 0)
            message(STATUS "Extracting ${file_name}.bz2")
            execute_process(COMMAND bunzip2 ${CMAKE_CURRENT_BINARY_DIR}/data/${file_name}.bz2
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
            message(STATUS "Done.")
            # list(APPEND remote_Files "${CMAKE_CURRENT_BINARY_DIR}/data/${file_name}")
            # set(remote_Files ${remote_Files} PARENT_SCOPE)
        else()
            message(WARNING "------------------------------------------------")
            message(WARNING "Error occurred during download: ${ERROR_MESSAGE}")
            message(WARNING "You should download ${file_name}.bz2 from https://hammer.physics.lbl.gov/mc_data/")
            message(WARNING "extract it and save it in ${CMAKE_INSTALL_DATADIR}/Hammer/examples/data")
            message(WARNING "after installing Hammer to be able to run the examples")
            message(WARNING "------------------------------------------------")
        endif()
    endif()
    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/data/${file_name})
        list(APPEND remote_Files "${CMAKE_CURRENT_BINARY_DIR}/data/${file_name}")
        set(remote_Files ${remote_Files} PARENT_SCOPE)
    endif()
endfunction()

function(new_example)
    set(options ROOT HEPMC PYTHON)
    set(oneValueArgs NAME)
    set(multiValueArgs REMOTE_FILES LOCAL_FILES)
    cmake_parse_arguments(PARSE_ARGV 0 arg "${options}" "${oneValueArgs}" "${multiValueArgs}")

    if(NOT arg_PYTHON)
        set(target_name ${arg_NAME})
        add_executable(${target_name})
        target_sources(${target_name} PRIVATE ${arg_NAME}.cc)

        if(arg_HEPMC)
           target_include_directories(${target_name} PRIVATE "${PROJECT_SOURCE_DIR}/Examples")
           #target_sources(${target_name} PRIVATE "HepMCParser.hh")   
        endif()

        set_target_properties(${target_name} PROPERTIES
            POSITION_INDEPENDENT_CODE ${Hammer_PIC})

        target_compile_features(${target_name} PRIVATE cxx_std_17)

        if(VERBOSE_DEBUG)
            target_compile_options(${target_name} PRIVATE ${Hammer_VerboseOptions})
        endif()

        if(SANITIZE)
            target_link_libraries(${target_name} PRIVATE HammerSanitizers)
        endif()

        # target_link_libraries(${target_name} INTERFACE "-Wl,-rpath,${PROJECT_BINARY_DIR}/lib")
        target_link_libraries(${target_name} PRIVATE HammerCore Amplitudes FormFactors HammerTools HammerBase Rates)

        if(arg_HEPMC)
            target_link_libraries(${target_name} PRIVATE HepMC3::All)
        endif()

        if(arg_ROOT)
            target_link_libraries(${target_name} PRIVATE ROOT::Physics ROOT::Hist ROOT::Tree ROOT::Gpad)
        endif()

    else()
        set(target_name "${arg_NAME}_python")
        add_custom_target(${target_name})
        target_sources(${target_name} PRIVATE ${arg_NAME}.py)

        # configure_file(${name}.py ${CMAKE_CURRENT_BINARY_DIR}/${name}.py @ONLY)
    endif()

    if(DEFINED arg_LOCAL_FILES AND arg_LOCAL_FILES)
        foreach(MY_FILE ${arg_LOCAL_FILES})
            list(APPEND extra_Files "${CMAKE_SOURCE_DIR}/Examples/data/${MY_FILE}")
        endforeach()

        add_dependencies(${target_name} ExampleFiles)
    endif()

    if(DEFINED arg_REMOTE_FILES AND arg_REMOTE_FILES)
        foreach(MY_FILE ${arg_REMOTE_FILES})
            download_example_file(${MY_FILE})
        endforeach()

        add_dependencies(${target_name} ExampleFiles)
    endif()

    set(extra_Files ${extra_Files} PARENT_SCOPE)
    set(remote_Files ${remote_Files} PARENT_SCOPE)

    if(NOT arg_PYTHON)
        install(TARGETS ${target_name}
            RUNTIME DESTINATION ${CMAKE_INSTALL_DATADIR}/Hammer/examples
            COMPONENT Hammer_Runtime
        )

        install(FILES ${arg_NAME}.cc
            DESTINATION ${CMAKE_INSTALL_DATADIR}/Hammer/examples
            COMPONENT Hammer_Development)
    else()
        install(PROGRAMS ${arg_NAME}.py
            DESTINATION ${CMAKE_INSTALL_DATADIR}/Hammer/examples
            COMPONENT Hammer_Runtime
        )
    endif()
endfunction()

if(HAVE_HEPMC3)
    install(FILES HepMCParser.hh
        DESTINATION ${CMAKE_INSTALL_DATADIR}/Hammer/examples
        COMPONENT Hammer_Development
    )

    new_example(NAME demo01 HEPMC REMOTE_FILES "BCLepNu.hepmc")
    new_example(NAME demo01card HEPMC REMOTE_FILES "BCLepNu.hepmc" LOCAL_FILES "Options.yml" "Processes.yml")
    new_example(NAME demo03 HEPMC REMOTE_FILES "BCLepNu.hepmc")
    new_example(NAME demo05 HEPMC REMOTE_FILES "BCTauNu_PS.hepmc")
    new_example(NAME demo07 HEPMC REMOTE_FILES "BDEllNu.hepmc")
    new_example(NAME demo09 HEPMC REMOTE_FILES "BCLepNu.hepmc")
    new_example(NAME demo11 HEPMC REMOTE_FILES "BDEllNu.hepmc")
    new_example(NAME demo13 HEPMC REMOTE_FILES "BCLepNu.hepmc")
    new_example(NAME demo13rl HEPMC REMOTE_FILES "BCLepNu.hepmc")

    if(HAVE_ROOT)
        new_example(NAME demo01root ROOT HEPMC REMOTE_FILES "BCLepNu.hepmc")
        new_example(NAME demo03root ROOT HEPMC REMOTE_FILES "BCLepNu.hepmc")
    endif()

    if(WITH_PYTHON)
        if(PYTHON_USE_CPPYY)
            install(FILES HepMCParser_pyhepmc.py HepMCParser_HepMC3.py HepMCParser_hepmc.py
                DESTINATION ${CMAKE_INSTALL_DATADIR}/Hammer/examples
                COMPONENT Hammer_Runtime
            )
            new_example(NAME demo01 PYTHON HEPMC REMOTE_FILES "BCLepNu.hepmc")
            new_example(NAME demo01_io PYTHON HEPMC REMOTE_FILES "BCLepNu.hepmc")
            new_example(NAME demo03 PYTHON HEPMC REMOTE_FILES "BCLepNu.hepmc")
            if(HAVE_ROOT)
                new_example(NAME demo01root PYTHON ROOT HEPMC REMOTE_FILES "BCLepNu.hepmc")
            endif()
        else()
            new_example(NAME demo01_cython PYTHON HEPMC REMOTE_FILES "BCLepNu.hepmc")
            new_example(NAME demo03_cython PYTHON HEPMC REMOTE_FILES "BCLepNu.hepmc")
        endif()
    endif()
endif()

new_example(NAME demo02)
new_example(NAME demo04)
new_example(NAME demo04parallel)
new_example(NAME demo04mma)
new_example(NAME demo08)
new_example(NAME demo08mma)
new_example(NAME demo10)
new_example(NAME demo12)
new_example(NAME demo14)

if(HAVE_ROOT)
    new_example(NAME demo02root ROOT)
    new_example(NAME demo04root ROOT)
    new_example(NAME demo04root2 ROOT)
endif()

if(WITH_PYTHON)
    if(PYTHON_USE_CPPYY)
        new_example(NAME demo02 PYTHON)
        new_example(NAME demo04 PYTHON)
    else()
        new_example(NAME demo02_cython PYTHON)
        new_example(NAME demo04_cython PYTHON)
    endif()
endif()

target_sources(ExampleFiles PRIVATE ${extra_Files} ${remote_Files})
install(FILES ${extra_Files}
    DESTINATION ${CMAKE_INSTALL_DATADIR}/Hammer/examples/data
    COMPONENT Hammer_Runtime
)
install(FILES ${remote_Files}
    DESTINATION ${CMAKE_INSTALL_DATADIR}/Hammer/examples/data
    COMPONENT Hammer_Runtime
)
