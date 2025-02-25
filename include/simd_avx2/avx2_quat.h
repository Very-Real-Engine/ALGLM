#pragma once

#include <cmath>
#include <cstdint>
#include "../simd/compiler_config.h"
#include <immintrin.h> 

namespace alglm
{

struct vec3;

struct ALIGN16 quat
{
	float x, y, z, w;

	quat();
	quat(float w, float x, float y, float z);
	quat(float w, const vec3 &v);
	quat(const vec3 &eulerAngle);
	quat(const quat &copy);
	quat &operator=(const quat &copy);
	quat operator*(const quat &rhs) const;
	quat operator+(const quat &rhs) const;
	quat operator-(const quat &rhs) const;
	float &operator[](int idx);
	float operator[](int idx) const;
};

float dot(const quat &q1, const quat &q2);
float length(const quat &q);
float length2(const quat &q);
float getPitch(const quat &q);
float getYaw(const quat &q);
float getRoll(const quat &q);
vec3 eulerAngles(const quat &q);
quat slerp(const quat &x, const quat &y, float a);
quat operator*(float a, const quat &q);
quat operator*(const quat &q, float a);
quat operator/(float a, const quat &q);
quat operator/(const quat &q, float a);
quat normalize(const quat &q);
quat angleAxis(float angle, const vec3 &axis);

} // namespace alglm