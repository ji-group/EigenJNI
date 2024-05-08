cmake_minimum_required(VERSION 3.23)
project(EigenJNI)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)
set(EIGEN_VERSION 3.4)
find_package(Eigen3 ${EIGEN_VERSION} QUIET)
if(NOT EIGEN3_FOUND)
    message("   Eigen3: library not found - downloading.")
    set(BUILD_TESTING OFF CACHE INTERNAL "")
    FetchContent_Declare(eigen
            GIT_REPOSITORY  https://gitlab.com/libeigen/eigen.git
            GIT_TAG         ${EIGEN_VERSION}
            GIT_SHALLOW     ON)
    FetchContent_MakeAvailable(eigen)
    unset(BUILD_TESTING CACHE)
    message("   Eigen3: downloaded.")
endif()

set(CMAKE_CXX_STANDARD 17)

find_package(JNI)
if(JNI_FOUND)
    message (STATUS "JAVA_HOME=$ENV{JAVA_HOME}")
    message (STATUS "JNI_INCLUDE_DIRS=${JNI_INCLUDE_DIRS}")
    message (STATUS "JNI_LIBRARIES=${JNI_LIBRARIES}")

    include_directories(
            $ENV{JAVA_HOME}/include
            ${JNI_INCLUDE_DIRS}
    )
else(JNI_FOUND)
    if(BUILD_JNI)
        message (FATAL_ERROR "No JNI includes and libraries found")
    endif(BUILD_JNI)
endif(JNI_FOUND)

add_library(eigen-jni SHARED
        dr_evomodel_substmodel_eigen_EigenJNIWrapper.cpp
        dr_evomodel_substmodel_eigen_EigenJNIWrapper.h)
target_link_libraries(eigen-jni Eigen3::Eigen)

add_executable(EigenJNI main.cpp)
target_link_libraries(EigenJNI Eigen3::Eigen)

if(APPLE)
    set_target_properties(eigen-jni PROPERTIES PREFIX "lib" SUFFIX ".jnilib")
elseif(WIN32)
    set_target_properties(eigen-jni PROPERTIES PREFIX "" SUFFIX "64.dll")
elseif(UNIX)
    set_target_properties(eigen-jni PROPERTIES PREFIX "lib" SUFFIX ".so")
endif(APPLE)

set(EIGEN_JNI_INSTALL_DIR lib)
install(TARGETS eigen-jni
        DESTINATION ${EIGEN_JNI_INSTALL_DIR}
        )


