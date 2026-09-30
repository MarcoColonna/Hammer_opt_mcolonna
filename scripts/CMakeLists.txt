message(STATUS "Preparing scripts")

if(BUILT_EXTERNAL_DEPENDENCIES)
    set(HAMMER_EXTRA_LIBPATH "${CMAKE_INSTALL_PREFIX}/lib:")
    set(HAMMER_EXTERNAL_LIBDIR "${CMAKE_INSTALL_PREFIX}/lib")
else()
    set(HAMMER_EXTRA_LIBPATH "")
    set(HAMMER_EXTERNAL_LIBDIR "")
endif()

configure_file(thishammer.sh.in ${CMAKE_CURRENT_BINARY_DIR}/thishammer.sh @ONLY)
configure_file(thishammer.csh.in ${CMAKE_CURRENT_BINARY_DIR}/thishammer.csh @ONLY)
configure_file(thishammer.fish.in ${CMAKE_CURRENT_BINARY_DIR}/thishammer.fish @ONLY)

add_custom_target(HammerScripts)

target_sources(HammerScripts PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/thishammer.sh
    ${CMAKE_CURRENT_BINARY_DIR}/thishammer.csh
    ${CMAKE_CURRENT_BINARY_DIR}/thishammer.fish)

install(PROGRAMS
    ${CMAKE_CURRENT_BINARY_DIR}/thishammer.sh
    ${CMAKE_CURRENT_BINARY_DIR}/thishammer.csh
    ${CMAKE_CURRENT_BINARY_DIR}/thishammer.fish
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT Hammer_Runtime
)

if(HAVE_ROOT)
    configure_file(thisrootalias.sh.in ${CMAKE_CURRENT_BINARY_DIR}/thisrootalias.sh @ONLY)
    
    add_custom_target(HammerRootScripts)
    
    target_sources(HammerRootScripts PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/thisrootalias.sh)
 
    install(PROGRAMS
        ${CMAKE_CURRENT_BINARY_DIR}/thisrootalias.sh
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT Hammer_Runtime
    )
    
endif()

