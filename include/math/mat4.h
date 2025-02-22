#pragma once

#include <cmath>
#include <cstdint>

namespace alglm
{

struct quat;
struct vec3;
struct vec4;

struct mat4
{
	float data[16];

	mat4();
	mat4(float x);
	mat4(const mat4 &copy);
	mat4(const quat &quat);
	mat4(const vec4 &v1, const vec4 &v2, const vec4 &v3, const vec4 &v4);
	mat4(float f0, float f1, float f2, float f3, float f4, float f5, float f6, float f7, float f8, float f9, float f10,
		 float f11, float f12, float f13, float f14, float f15);
	mat4 &operator=(const mat4 &copy);
	mat4 operator+(const mat4 &rhs) const;
	mat4 operator-(const mat4 &rhs) const;
	mat4 operator*(const mat4 &rhs) const;
	float *operator[](int idx);
	const float *operator[](int idx) const;
};

mat4 operator*(float scalar, const mat4 &matrix);
mat4 operator*(const mat4 &matrix, float scalar);

float *value_ptr(mat4 &matrix);

mat4 inverse(const mat4 &matrix);
mat4 scale(const mat4 &matrix, const vec3 &vector);
mat4 translate(const mat4 &matrix, const vec3 &vector);
mat4 rotate(const mat4 &matrix, float theta, const vec3 &vector);
mat4 perspective(float fovy, float aspect, float zNear, float zFar);
mat4 lookAt(const vec3 &cameraPos, const vec3 &cameraTarget, const vec3 &cameraUp);
mat4 toMat4(const quat &quat);
mat4 ortho(float left, float right, float bottom, float top, float zNear, float zFar);
bool decompose(const mat4 &modelMatrix, vec3 &translation, quat &rotation, vec3 &scale, vec3 &skew, vec4 &perspective);
quat quat_cast(const mat4 &matrix);

} // namespace alglm