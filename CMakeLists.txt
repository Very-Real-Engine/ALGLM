cmake_minimum_required(VERSION 3.20)
project(alglm LANGUAGES CXX)

add_library(alglm STATIC
    src/simd_avx2/avx2_vec2.cpp
    src/simd_avx2/avx2_vec3.cpp
    src/simd_avx2/avx2_vec4.cpp
    src/simd_avx2/avx2_mat3.cpp
	src/simd_avx2/avx2_mat4.cpp
	src/simd_avx2/avx2_quat.cpp
	src/simd_avx2/avx2_math_utils.cpp
)

# Debug 빌드 시 라이브러리 이름에 'd' 접미사를 추가합니다.
set_target_properties(alglm PROPERTIES
    DEBUG_POSTFIX "d"
)

# alglm 설치 구성
install(TARGETS alglm
    EXPORT alglmTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION include/alglm
    FILES_MATCHING PATTERN "*.h"
)

# 설치용 CMake 파일 추가
install(EXPORT alglmTargets
    FILE alglmConfig.cmake
    DESTINATION lib/cmake/alglm
)