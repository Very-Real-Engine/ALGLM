#pragma once

#include <cmath>
#include <cstdint>
#include <algorithm>

namespace alglm
{

template <typename T> constexpr T pi()
{
	if constexpr (std::is_same_v<T, float>)
	{
		return 3.14159265358979323846f; // float 타입
	}
	else
	{
		return 3.14159265358979323846; // 기본 double 타입
	}
}

float radians(float degree)
{
	return degree * (3.14159f / 180.0f);
}

float degrees(float radian)
{
	return radian * (180.0f / 3.14159f);
}

float abs(float value)
{
	return (value < 0 ? -value : value);
}

float clamp(float x, float minVal, float maxVal)
{
	return std::max(minVal, std::min(x, maxVal));
}

vec3 radians(vec3 v);
vec3 degrees(vec3 v);

} // namespace alglm
