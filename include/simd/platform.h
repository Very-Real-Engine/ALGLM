#ifndef PLATFORM_H
#define PLATFORM_H

// CPU 아키텍처에 따라 적절한 SIMD 헤더 파일을 포함
#if defined(__AVX512F__)
    #define CPU_AVX512
    #include "simd_avx512.h"
#elif defined(__AVX2__)
    #define CPU_AVX2
    #include "simd_avx2.h"
#elif defined(__SSE4_2__)
    #define CPU_SSE42
    #include "simd_sse42.h"
#elif defined(__ARM_NEON)
    #define CPU_NEON
    #include "simd_neon.h"
#else
    #define CPU_SCALAR
    #include "simd_scalar.h"
#endif

#endif // PLATFORM_H