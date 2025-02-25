#pragma once

#include <cmath>
#include <cstdint>

namespace alglm
{

struct vec4;

struct vec3
{
	float x, y, z;

	vec3();
	vec3(float x);
	vec3(float x, float y, float z);
	vec3(const vec3 &copy);
	vec3(const vec4 &copy);
	vec3 &operator=(const vec3 &copy);
	vec3 &operator=(float copy);
	vec3 operator+(const vec3 &rhs) const;
	vec3 operator-(const vec3 &rhs) const;
	vec3 operator*(const vec3 &rhs) const;
	float &operator[](int idx);
	float operator[](int idx) const;
};

vec3 operator*(float scalar, const vec3 &vector);
vec3 operator*(const vec3 &vector, float scalar);
vec3 operator/(float scalar, const vec3 &vector);
vec3 operator/(const vec3 &vector, float scalar);

float dot(const vec3 &vector1, const vec3 &vector2);
float length(const vec3 &vector);
float length2(const vec3 &vector);
float *value_ptr(vec3 &vector);

vec3 cross(const vec3 &vector1, const vec3 &vector2);
vec3 normalize(const vec3 &vector);
vec3 min(const vec3 &v1, const vec3 &v2);
vec3 max(const vec3 &v1, const vec3 &v2);
vec3 mix(const vec3 &x, const vec3 &y, float a);

} // namespace alglm