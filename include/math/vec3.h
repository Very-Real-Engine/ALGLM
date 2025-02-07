#pragma once

#include <cmath>
#include <stdexcept>

namespace alglm
{

class vec4;

class vec3
{
  public:
	float x;
	float y;
	float z;

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

float dot(const vec3 &vector1, const vec3 &vector2);
float length(const vec3 &vector);
float *value_ptr(vec3 &vector);

vec3 cross(const vec3 &vector1, const vec3 &vector2);
vec3 normalize(const vec3 &vector);

} // namespace alglm