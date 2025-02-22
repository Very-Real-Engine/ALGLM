#pragma once

#include <cmath>
#include <cstdint>

namespace alglm
{

struct vec3;

struct quat
{
	float x, y, z, w;

	quat();
	quat(float x);
	quat(float x, float y, float z, float w);
	quat(const vec3 &v, float w);
	quat(float angle, const vec3 &axis);
	quat(const vec3 &eulerAngle);
	quat(const quat &copy);
	quat &operator=(const quat &copy);
	quat operator*(const quat &rhs) const;
	quat operator+(const quat &rhs) const;
	quat operator-(const quat &rhs) const;
};

float dot(const quat &q1, const quat &q2);
float length(const quat &q);
float length2(const quat &q);
vec3 eulerAngles(const quat &q);
quat slerp(const quat &x, const quat &y, float a);
quat operator*(float a, const quat &q);
quat operator*(const quat &q, float a);
quat normalize(const quat &q);
quat angleAxis(float angle, const vec3 &axis);

} // namespace alglm