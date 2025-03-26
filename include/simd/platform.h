#ifndef PLATFORM_H
#define PLATFORM_H

// CPU 아키텍처에 따라 적절한 SIMD 헤더 파일을 포함

// #ifdef __AVX2__
// #include "simd_avx2.h"
// #else
//     #include "simd_scalar.h"
// #endif
#include "simd_scalar.h"

#endif // PLATFORM_H