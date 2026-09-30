message(STATUS "Preparing Python extensions")

if(PYTHON_USE_CPPYY)
  add_subdirectory(wrapper_cppyy)
  set(WDIR "${CMAKE_CURRENT_BINARY_DIR}/wrapper_cppyy/bindings")

# cppyy_add_bindings(
# "pyhammer" "${PROJECT_VERSION}" "HAMMER Collaboration" ""
# LICENSE "GPL3; MCnet"
# LANGUAGE_STANDARD "17"
# SELECTION_XML ${CMAKE_CURRENT_SOURCE_DIR}/selection.xml
# INTERFACE_FILE ${CMAKE_SOURCE_DIR}/interface.hh
# HEADERS ${LIB_HEADERS} ${CMAKE_SOURCE_DIR}/interface.hh
# INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include
# LINK_LIBRARIES "HammerBase HammerCore HammerTools Amplitudes FormFactors Rates"
# NAMESPACES Hammer
# )

# install(TARGETS pyhammer
# LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
# PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Hammer
# )

# install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install ${PY_WHEEL_FILE})")
else()
  add_subdirectory(wrapper_cython)
  set(WDIR "${CMAKE_CURRENT_BINARY_DIR}/wrapper_cython/bindings")

  # set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
  # set(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  # set(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

  # FILE(GLOB DEPS "${CMAKE_CURRENT_SOURCE_DIR}/hammer")
  # SEPARATE_ARGUMENTS(DEPS)

  # file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/hammer DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  # configure_file(${SETUP_PY_IN} ${SETUP_PY} ESCAPE_QUOTES @ONLY)

  # add_custom_command(OUTPUT ${OUTPUT}
  # COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} build
  # COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
  # DEPENDS ${DEPS} hammerlib)

  # add_custom_target(PythonPackages ALL DEPENDS ${OUTPUT})

  # install(CODE "execute_process(
  # COMMAND ${Python3_EXECUTABLE} -m pip install .
  # --config-settings prefix=${CMAKE_INSTALL_PREFIX}
  # --config-settings install-lib=${PYTHON_INSTALL_DIR}
  # --config-settings single-version-externally-managed=1
  # --config-settings record=${CMAKE_CURRENT_BINARY_DIR}/build/filelist.txt
  # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
endif()

add_custom_target(SetupTools ALL DEPENDS PythonPackages)
add_custom_command(TARGET SetupTools POST_BUILD
  COMMAND ${Python3_EXECUTABLE} "-m" "build" "--wheel" WORKING_DIRECTORY "${WDIR}")

# need to find the actual path without relying on pip to avoid the local non-local snafu on certain linux platforms.
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c
        "import sysconfig, os; print(os.path.relpath(sysconfig.get_path('purelib', 'posix_prefix', {'base': '/p', 'platbase': '/p'}), '/p'))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES_SUBDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

install(CODE "execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip
        install --no-deps --no-index -f ${WDIR}/dist --target=\${CMAKE_INSTALL_PREFIX}/${PYTHON_SITE_PACKAGES_SUBDIR} hammer
        WORKING_DIRECTORY ${WDIR})")
