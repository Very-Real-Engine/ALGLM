#pragma once

#include <cmath>
#include <cstdint>

namespace alglm
{

struct vec2
{
	float x, y;

	vec2();
	vec2(float x);
	vec2(float x, float y);
	vec2(const vec2 &copy);
	vec2 &operator=(const vec2 &copy);
	vec2 operator+(const vec2 &rhs) const;
	vec2 operator-(const vec2 &rhs) const;
	vec2 operator*(const vec2 &rhs) const;
	float &operator[](int idx);
	float operator[](int idx) const;
};

vec2 operator*(float scalar, const vec2 &vector);
vec2 operator*(const vec2 &vector, float scalar);

float dot(const vec2 &vector1, const vec2 &vector2);
float length(const vec2 &vector);
float length2(const vec2 &vector);
float *value_ptr(vec2 &vector);

vec2 normalize(const vec2 &vector);

} // namespace alglm