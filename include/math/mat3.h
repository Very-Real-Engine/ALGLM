#pragma once

#include <cmath>
#include <cstdint>

namespace alglm
{

struct mat4;
struct vec3;

struct mat3
{
	float data[9];

	mat3();
	mat3(float x);
	mat3(const mat3 &copy);
	mat3(const mat4 &copy);
	mat3(float f0, float f1, float f2, float f3, float f4, float f5, float f6, float f7, float f8);
	mat3(const vec3 &v1, const vec3 &v2, const vec3 &v3);
	mat3 &operator=(const mat3 &copy);
	mat3 operator+(const mat3 &rhs) const;
	mat3 operator-(const mat3 &rhs) const;
	mat3 operator*(const mat3 &rhs) const;
	float *operator[](int idx);
	const float *operator[](int idx) const;
};

mat3 operator*(float scalar, const mat3 &matrix);
mat3 operator*(const mat3 &matrix, float scalar);
mat3 inverse(const mat3 &matrix);
mat3 transpose(const mat3 &matrix);
float *value_ptr(mat3 &matrix);
float determinant(const mat3 &m);

} // namespace alglm