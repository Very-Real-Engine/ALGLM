#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace alglm
{

struct vec3;
struct vec4;

template <typename T> constexpr T pi()
{
	if constexpr (std::is_same<T, float>::value)
	{
		return 3.1415927f;
	}
	else
	{
		return 3.141592653589793;
	}
}

float abs(float value);
float clamp(float x, float minVal, float maxVal);

float radians(float degree);
vec3 radians(const vec3 &v);

float degrees(float radian);
vec3 degrees(const vec3 &v);

float mix(float x, float y, float a);
vec4 mix(const vec4& x, const vec4& y, float a);

} // namespace alglm
