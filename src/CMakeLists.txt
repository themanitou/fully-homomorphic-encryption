project(FULLHOMCRYPT)

cmake_minimum_required(VERSION 3.3)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(FULLHOMCRYPT_VERSION 1.0.0)
string(REGEX MATCH "^[0-9]+" FULLHOMCRYPT_SOVERSION ${FULLHOMCRYPT_VERSION})

set(QTPATH /home/ngu145/Workspace/git/workspace/Tools/ThirdParty)

find_package(MPI)
find_package(Qt5Core PATHS ${QTPATH})
find_package(Qt5Network PATHS ${QTPATH})
set(QT_LIBRARIES Qt5::Core;Qt5::Network)

include_directories(${MPI_CXX_INCLUDE_DIRS})

set(CMAKE_AUTOMOC ON)

set(FULLHOMCRYPT_HEADERS
  ${FULLHOMCRYPT_SOURCE_DIR}/calculation/fullhom.h
  ${FULLHOMCRYPT_SOURCE_DIR}/calculation/multiplication.h
  ${FULLHOMCRYPT_SOURCE_DIR}/calculation/tools.h
  ${FULLHOMCRYPT_SOURCE_DIR}/communication/fhe_mpi.h
  ${FULLHOMCRYPT_SOURCE_DIR}/communication/commserver.h
  ${FULLHOMCRYPT_SOURCE_DIR}/communication/commclient.h
)

set(FULLHOMCRYPT_SOURCES
  ${FULLHOMCRYPT_SOURCE_DIR}/calculation/fullhom.cpp
  ${FULLHOMCRYPT_SOURCE_DIR}/calculation/multiplication.cpp
  ${FULLHOMCRYPT_SOURCE_DIR}/calculation/tools.cpp
  ${FULLHOMCRYPT_SOURCE_DIR}/communication/fhe_mpi.cpp
  ${FULLHOMCRYPT_SOURCE_DIR}/communication/commserver.cpp
  ${FULLHOMCRYPT_SOURCE_DIR}/communication/commclient.cpp
  ${FULLHOMCRYPT_SOURCE_DIR}/FullHomCrypt.cpp
)

add_definitions(-DFULLHOMCRYPT_VERSION=${FULLHOMCRYPT_VERSION} -DFHE_TEST -DFHE_MPI_LOG)

add_executable(fullhomcrypt ${FULLHOMCRYPT_SOURCES} ${FULLHOMCRYPT_HEADERS} ${FULLHOMCRYPT_MOC_SOURCES})
target_link_libraries(fullhomcrypt ntl ${MPI_CXX_LIBRARIES} ${QT_LIBRARIES})

set_target_properties(fullhomcrypt PROPERTIES
    VERSION ${FULLHOMCRYPT_VERSION}
    SOVERSION ${FULLHOMCRYPT_SOVERSION}
)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_LIST_DIR}/../bin CACHE STRING " " FORCE )
install(TARGETS fullhomcrypt DESTINATION ${CMAKE_INSTALL_PREFIX})

source_group(Sources FILES ${FULLHOMCRYPT_SOURCES})
