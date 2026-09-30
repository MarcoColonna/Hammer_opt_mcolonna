set(YamlCppVersion 0.8.0)
set(YamlCppSHA256 fbe74bbdcee21d656715688706da3c8becfd946d92cd44705cc6098bb23b3a16)

include(FetchContent)

set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
set(YAML_CPP_INSTALL ON CACHE BOOL "" FORCE)

FetchContent_Declare(yaml-cpp
    URL https://github.com/jbeder/yaml-cpp/archive/refs/tags/${YamlCppVersion}.tar.gz
    URL_HASH SHA256=${YamlCppSHA256}
    DOWNLOAD_DIR ${PROJECT_BINARY_DIR}/ThirdParty/tmp
    SOURCE_DIR ${PROJECT_BINARY_DIR}/ThirdParty/YamlCpp
    BINARY_DIR ${PROJECT_BINARY_DIR}/ThirdParty/YamlCpp-build
)
FetchContent_MakeAvailable(yaml-cpp)

set(BUILT_EXTERNAL_DEPENDENCIES ON CACHE BOOL "" FORCE)
