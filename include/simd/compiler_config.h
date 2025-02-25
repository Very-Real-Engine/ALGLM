#ifndef COMPILER_CONFIG_H
#define COMPILER_CONFIG_H

// 컴파일러별 인라인 처리
#ifdef _MSC_VER
    // MSVC에서는 강제 인라인: __forceinline
    #define FORCE_INLINE __forceinline
    // MSVC에서는 alignas 대신 __declspec(align(N)) 사용
    #define ALIGN4 __declspec(align(4))
    #define ALIGN8 __declspec(align(8))
    #define ALIGN16 __declspec(align(16))
    #define ALIGN32 __declspec(align(32))
    #define ALIGN64 __declspec(align(64))
#elif defined(__GNUC__) || defined(__clang__)
    // GCC, Clang에서는 __attribute__((always_inline))와 inline
    #define FORCE_INLINE __attribute__((always_inline)) inline
    // GCC, Clang에서는 __attribute__((aligned(N))) 사용
    #define ALIGN16 __attribute__((aligned(16)))
    #define ALIGN32 __attribute__((aligned(32)))
    #define ALIGN64 __attribute__((aligned(64)))
#else
    #define FORCE_INLINE inline
    #define ALIGN16 alignas(16)
    #define ALIGN32 alignas(32)
    #define ALIGN64 alignas(64)
#endif

#endif // COMPILER_CONFIG_H