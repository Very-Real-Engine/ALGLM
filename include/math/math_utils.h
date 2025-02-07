#pragma once

#include "simd/compiler_config.h"

namespace alglm
{
constexpr float pi = 3.14159f;

FORCE_INLINE float radians(float degree)
{
	return degree * (3.14159f / 180.0f);
}

} // namespace alglm
