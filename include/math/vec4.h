#pragma once

#include <cmath>
#include <stdexcept>

namespace alglm
{
class mat4;

class vec4
{
  public:
	float x;
	float y;
	float z;
	float w;

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

vec4 operator*(float scalar, const vec4 &vector);
vec4 operator*(const vec4 &vector, float scalar);
vec4 operator*(const mat4 &matrix, const vec4 &vector);

float *value_ptr(vec4 &vector);

float dot(const vec4 &vector1, const vec4 &vector2);
vec4 normalize(const vec4 &vector);


} // namespace alglm