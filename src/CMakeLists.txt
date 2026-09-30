message(STATUS "Preparing Hammer libraries.")

add_subdirectory(Tools)

add_subdirectory(Math)

add_subdirectory(Core)
add_subdirectory(Amplitudes)
add_subdirectory(FormFactors)
add_subdirectory(Rates)

add_library(Hammer INTERFACE)
target_link_libraries(Hammer INTERFACE HammerBase HammerTools HammerConfig HammerMath HammerTensor HammerCore Amplitudes FormFactors Rates)

install(TARGETS Hammer
    EXPORT HammerTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT Hammer_Runtime
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
    COMPONENT Hammer_Runtime
    NAMELINK_COMPONENT Hammer_Development
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/Hammer
    COMPONENT Hammer_Development
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Hammer
    COMPONENT Hammer_Development
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
