#ifndef PLATFORM_H
#define PLATFORM_H

// CPU 아키텍처에 따라 적절한 SIMD 헤더 파일을 포함
// #if defined(__AVX2__)
//     #define CPU_AVX2
//     #include "simd_avx2.h"
// #else
//     #define CPU_SCALAR
//     #include "simd_scalar.h"
// #endif
#include "simd_scalar.h"

#endif // PLATFORM_H