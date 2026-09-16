cmake_minimum_required(VERSION 3.10)
project(BosonTest)

# std=c++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# default Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED)

# 查找 SuperLU（可能需要手动指定）
# find_package(SuperLU)  # 如果 CMake 能找到的话

# 添加可执行文件
add_executable(FMKH
    FMKH.cpp
    src/Boson.cpp
    src/lapack.cpp
    src/SpinWave.cpp
    src/spin.cpp
    src/MeanField.cpp
    src/SelfEnergy.cpp
)

# 设置头文件目录
target_include_directories(FMKH PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# 链接库
target_link_libraries(FMKH
    LAPACK::LAPACK
    BLAS::BLAS
    OpenMP::OpenMP_CXX
    gfortran
    m
)

# 关键：添加优化选项
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(FMKH PRIVATE 
        -Wall 
        -Wextra 
        -Ofast
        -march=native    
        -fopenmp         
        -ffast-math      
    )
    
    target_link_options(FMKH PRIVATE -fopenmp)
    
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(FMKH PRIVATE /W4 /O3 /openmp)
endif()