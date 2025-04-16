cmake_minimum_required(VERSION 3.14)
project(DDG_Tutorial_Fast)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")

# Per https://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html eigen should be used
# with MKL LAPACK only with lp64 set. The default may be ILP64, which we don't want.
set(MKL_INTERFACE_FULL intel_lp64)

# Define a flag (string) to choose which case to build
set(TARGET_BUILD "sim" CACHE STRING "Choose which target to build: 2d_curve, 2d_surface, 3d_curve, 3d_surface, hollow_net")


# List of allowed values
set(ALLOWED_BUILDS 2d_curve 2d_surface 3d_curve 3d_surface hollow_net)

# Validate user input
list(FIND ALLOWED_BUILDS ${TARGET_BUILD} build_index)
if(build_index EQUAL -1)
    message(FATAL_ERROR "Invalid TARGET_BUILD: ${TARGET_BUILD}. Choose one of: 2d_curve, 2d_surface, 3d_curve, 3d_surface, hollow_net.")
endif()

add_definitions(-DEIGEN_USE_MKL_ALL)

# Dependencies
find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(MKL CONFIG REQUIRED)


include_directories(
        ${EIGEN_INCLUDE_DIRS}
        ${OPENGL_INCLUDE_DIRS}
        ${GLUT_INCLUDE_DIRS}
        )



add_executable(simDER
        src/${TARGET_BUILD}/main.cpp
        src/${TARGET_BUILD}/world.cpp
        src/${TARGET_BUILD}/setInput.cpp
        src/${TARGET_BUILD}/timeStepper.cpp

        src/${TARGET_BUILD}/dampingForce.cpp
        src/${TARGET_BUILD}/elasticPlate.cpp
        src/${TARGET_BUILD}/elasticStretchingForce.cpp
        src/${TARGET_BUILD}/externalGravityForce.cpp
        src/${TARGET_BUILD}/inertialForce.cpp
        )

# Conditionally add twistingForce.cpp if TARGET_BUILD equals "xxx"
if(TARGET_BUILD STREQUAL "3d_curve")
    target_sources(simDER PRIVATE src/${TARGET_BUILD}/elasticTwistingForce.cpp)
endif()

if (NOT TARGET_BUILD STREQUAL "hollow_net")
    target_sources(simDER PRIVATE src/${TARGET_BUILD}/elasticBendingForce.cpp)
endif()


target_include_directories(simDER PUBLIC
        ${CMAKE_SOURCE_DIR}/src/${TARGET_BUILD}
)


target_link_libraries(simDER
        ${OPENGL_LIBRARIES}
        ${GLUT_LIBRARIES}
        MKL::MKL
        Eigen3::Eigen
        lapack
        gfortran
        pthread
        rt
        m)