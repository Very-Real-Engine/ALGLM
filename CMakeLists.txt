cmake_minimum_required(VERSION 3.20)
project(alglm LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(alglm STATIC
    src/vec2.cpp
    src/vec3.cpp
    src/vec4.cpp
    src/mat3.cpp
	src/mat4.cpp
	src/quat.cpp
	src/math_utils.cpp
)

# include 디렉토리 설정
target_include_directories(alglm PUBLIC include)

# 필요하면 C++ 표준 설정
target_compile_features(alglm PUBLIC cxx_std17)
