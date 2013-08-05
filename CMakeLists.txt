cmake_minimum_required ( VERSION 2.6 )
project (raytracer)

include(FindROOT.cmake)
include(FindISPC.cmake)
include(FindGSL.cmake)

set(alignment 32)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath,./src/ROOT_Ttree_io")

set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -std=c++0x")
include_directories(include ${ROOT_INCLUDE_DIR} ${GSL_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR} ${GSL_LINK_DIRECTORIES} /home/vitillo/Dropbox/ispc_test/lib)

if(NOT DEFINED USE_DOUBLE)
  set(use_double -DUSE_DOUBLE=0)
else()
  set(use_double -DUSE_DOUBLE=${USE_DOUBLE})
endif()
add_definitions(-DALIGNMENT=${alignment} -DHAS_MATERIAL_EFFECTS=0 ${use_double})

#TODO: align memory to 32 bytes independently of the numbers of tracks
#ispc_compile(${CMAKE_SOURCE_DIR}/ispc/KalmanFilter.ispc "--target=avx-x2 ${use_double} --opt=force-aligned-memory --force-alignment=${alignment}" KalmanFilterISPC)

ispc_compile(${CMAKE_SOURCE_DIR}/ispc/KalmanFilter.ispc "--target=avx-x2 ${use_double}" KalmanFilterISPC)

add_executable(unit_tests src/unit_tests.cpp src/serial/KalmanFilterSerial.cpp include/ispc/KalmanFilter.ispc.h)
target_link_libraries(unit_tests ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${KalmanFilterISPC})

add_executable(kalman_benchmark src/kalman_benchmark.cpp src/serial/KalmanFilterSerial.cpp include/ispc/KalmanFilter.ispc.h)
target_link_libraries(kalman_benchmark ${GSL_LIBRARIES} ${ROOT_LIBRARIES} RootIoWrapper ${KalmanFilterISPC})

add_executable(synthetic_benchmark src/synthetic_benchmark.cpp src/serial/gemm.cpp src/serial/KalmanFilterSerial.cpp include/ispc/KalmanFilter.ispc.h)
target_link_libraries(synthetic_benchmark ${GSL_LIBRARIES} ${ROOT_LIBRARIES} RootIoWrapper ${KalmanFilterISPC})

add_subdirectory(src/ROOT_Ttree_io)

enable_testing()

add_test(unit_tests ./unit_tests)
add_test(synthetic_benchmark ./synthetic_benchmark)
