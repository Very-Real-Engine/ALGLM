#pragma once

#include <cmath>
#include <cstdint>
#include "../simd/compiler_config.h"

namespace alglm
{
struct mat4;
struct vec3;

struct ALIGN16 vec4
{
	union
	{
		struct {float x, y, z, w;};
		struct {float r, g, b, a;};
	};

	vec4();
	vec4(float x);
	vec4(float x, float y, float z, float w);
	vec4(const vec3 &copy, float w);
	vec4(const vec4 &copy);
	vec4 &operator=(const vec4 &copy);
	vec4 operator+(const vec4 &rhs) const;
	vec4 operator-(const vec4 &rhs) const;
	vec4 operator*(const vec4 &rhs) const;
	float &operator[](int idx);
	float operator[](int idx) const;
};

vec4 normalize(const vec4 &vector);
vec4 operator*(float scalar, const vec4 &vector);
vec4 operator*(const vec4 &vector, float scalar);
vec4 operator*(const mat4 &matrix, const vec4 &vector);

float dot(const vec4 &vector1, const vec4 &vector2);
float length(const vec4 &vector);
float length2(const vec4 &vector);
float *value_ptr(vec4 &vector);



struct ALIGN16 ivec4
{
	int32_t x, y, z, w;

	ivec4();
	ivec4(int32_t x);
	ivec4(int32_t x, int32_t y, int32_t z, int32_t w);
	ivec4(const ivec4 &copy);
	ivec4 &operator=(const ivec4 &copy);
	ivec4 operator+(const ivec4 &rhs) const;
	ivec4 operator-(const ivec4 &rhs) const;
	ivec4 operator*(const ivec4 &rhs) const;
};

} // namespace alglm