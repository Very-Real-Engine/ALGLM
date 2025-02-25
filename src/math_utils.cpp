#include "../include/alglm.h"

namespace alglm
{

float radians(float degree)
{
	return degree * (3.1415927f / 180.0f);
}

float degrees(float radian)
{
	return radian * (180.0f / 3.1415927f);
}

float abs(float value)
{
	return (value < 0 ? -value : value);
}

float clamp(float x, float minVal, float maxVal)
{
	return std::max(minVal, std::min(x, maxVal));
}

vec3 radians(const vec3 &v)
{
	return vec3(radians(v.x), radians(v.y), radians(v.z));
}

vec3 degrees(const vec3 &v)
{
	return vec3(degrees(v.x), degrees(v.y), degrees(v.z));
}

float mix(float x, float y, float a)
{
	return x * (1.0f - a) + y * a;
}

vec3 mix(const vec3 &x, const vec3 &y, float a)
{
	return x * (1.0f - a) + y * a;
}

} // namespace alglm
