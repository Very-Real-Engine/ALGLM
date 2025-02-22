#include "alglm.h"

namespace alglm
{
// member function

quat::quat() : x(0), y(0), z(0), w(1) {};

quat::quat(float x) : x(x), y(x), z(x), w(x) {};

quat::quat(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {};

quat::quat(const vec3 &v, float w) : x(v.x), y(v.y), z(v.z), w(w) {};

quat::quat(const quat &copy) : x(copy.x), y(copy.y), z(copy.z), w(copy.w) {};

quat &quat::operator=(const quat &copy)
{
	this->x = copy.x;
	this->y = copy.y;
	this->z = copy.z;
	this->w = copy.w;
	return *this;
}

quat::quat(float angle, const vec3 &axis)
{
	vec3 normalAxis = normalize(axis);
	float halfAngle = angle / 2.0f;
	float s = sin(halfAngle);
	x = normalAxis.x * s;
	y = normalAxis.y * s;
	z = normalAxis.z * s;
	w = cos(halfAngle);
}

quat::quat(const vec3 &eulerAngle)
{
	quat xQuat(vec3(1.0f, 0.0f, 0.0f), radians(eulerAngle.x));
	quat yQuat(vec3(0.0f, 1.0f, 0.0f), radians(eulerAngle.y));
	quat zQuat(vec3(0.0f, 0.0f, 1.0f), radians(eulerAngle.z));
	quat mulQuat = xQuat * yQuat * zQuat;

	x = mulQuat.x;
	y = mulQuat.y;
	z = mulQuat.z;
	w = mulQuat.w;
}

quat quat::operator*(const quat &rhs) const
{
	return quat(rhs.w * x + rhs.x * w + rhs.y * z - rhs.z * y, rhs.w * y - rhs.x * z + rhs.y * w + rhs.z * x,
				rhs.w * z + rhs.x * y - rhs.y * x + rhs.z * w, rhs.w * w - rhs.x * x - rhs.y * y - rhs.z * z);
}

quat quat::operator+(const quat &rhs) const
{
	return quat(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w);
}

quat quat::operator-(const quat &rhs) const
{
	return quat(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
}

quat normalize(const quat &q)
{
	quat ret(q);

	float length = sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
	ret.x = ret.x / length;
	ret.y = ret.y / length;
	ret.z = ret.z / length;
	ret.w = ret.w / length;

	return ret;
}

float dot(const quat &q1, const quat &q2)
{
	return q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;
}

float length(const quat &q)
{
	return sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
}

float length2(const quat &q)
{
	return q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w;
}

quat angleAxis(float angle, const vec3 &axis)
{
	vec3 normAxis = normalize(axis);
	float halfAngle = angle * 0.5f;
	float s = sin(halfAngle);

	return quat(normAxis.x * s, normAxis.y * s, normAxis.z * s, cos(halfAngle));
}

quat operator*(float a, const quat &q)
{
	return quat(a * q.x, a * q.y, a * q.z, a * q.w);
}

quat operator*(const quat &q, float a)
{
	return quat(a * q.x, a * q.y, a * q.z, a * q.w);
}

quat slerp(const quat& x, const quat& y, float a) {

	float dotRes = dot(x, y);

    quat yAdjusted = y;
    if (dotRes < 0.0f) {
        yAdjusted = -1 * y;
        dotRes = -dotRes;
    }

    const float THRESHOLD = 0.9995f;
    if (dotRes > THRESHOLD) {
        return normalize(x + a * (yAdjusted - x)); // LERP 후 정규화
    }

    float theta = std::acos(dotRes);
    float sinTheta = std::sin(theta);

    float w1 = std::sin((1.0f - a) * theta) / sinTheta;
    float w2 = std::sin(a * theta) / sinTheta;

    return normalize(w1 * x + w2 * yAdjusted);
}


} // namespace alglm
