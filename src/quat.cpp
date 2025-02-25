#include "../include/alglm.h"

namespace alglm
{
// member function

quat::quat() : x(0), y(0), z(0), w(1) {};

quat::quat(float w, float x, float y, float z) : x(x), y(y), z(z), w(w) {};

quat::quat(float w, const vec3 &v) : x(v.x), y(v.y), z(v.z), w(w) {};

quat::quat(const quat &copy) : x(copy.x), y(copy.y), z(copy.z), w(copy.w) {};

quat &quat::operator=(const quat &copy)
{
	this->x = copy.x;
	this->y = copy.y;
	this->z = copy.z;
	this->w = copy.w;
	return *this;
}

quat::quat(const vec3 &eulerAngle)
{
    vec3 half = eulerAngle * 0.5f;

    vec3 c(std::cos(half.x), std::cos(half.y), std::cos(half.z));
    vec3 s(std::sin(half.x), std::sin(half.y), std::sin(half.z));

    this->w = c.x * c.y * c.z + s.x * s.y * s.z;
    this->x = s.x * c.y * c.z - c.x * s.y * s.z;
    this->y = c.x * s.y * c.z + s.x * c.y * s.z;
    this->z = c.x * c.y * s.z - s.x * s.y * c.z;
}

quat quat::operator*(const quat &rhs) const
{
	return quat(rhs.w * w - rhs.x * x - rhs.y * y - rhs.z * z, rhs.w * x + rhs.x * w - rhs.y * z + rhs.z * y,
				rhs.w * y + rhs.x * z + rhs.y * w - rhs.z * x, rhs.w * z - rhs.x * y + rhs.y * x + rhs.z * w);
}

quat quat::operator+(const quat &rhs) const
{
	return quat(w + rhs.w, x + rhs.x, y + rhs.y, z + rhs.z);
}

quat quat::operator-(const quat &rhs) const
{
	return quat(w - rhs.w, x - rhs.x, y - rhs.y, z - rhs.z);
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
	float halfAngle = angle * 0.5f;
	float s = sin(halfAngle);

	return quat(cos(halfAngle), axis.x * s, axis.y * s, axis.z * s);
}

quat operator*(float a, const quat &q)
{
	return quat(a * q.w, a * q.x, a * q.y, a * q.z);
}

quat operator*(const quat &q, float a)
{
	return quat(a * q.w, a * q.x, a * q.y, a * q.z);
}

quat operator/(float a, const quat &q)
{
	return quat(q.w / a, q.x / a, q.y / a, q.z / a);
}

quat operator/(const quat &q, float a)
{
	return quat(q.w / a, q.x / a, q.y / a, q.z / a);
}

float getPitch(const quat &q)
{
	float y = 2.0f * (q.w * q.x + q.y * q.z);
    float x = (q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z);

	if (abs(y - x) <= 1e-6f)
	{
		return 2 * atan2(q.x, q.w);
	}
	return atan2(y, x);
}

float getYaw(const quat &q)
{
	return asin(clamp(-2.0f * (q.x * q.z - q.w * q.y), -1.0f, 1.0f));
}

float getRoll(const quat &q)
{
	float y = 2.0f * (q.x * q.y + q.w * q.z);
	float x = q.w * q.w + q.x * q.x - q.y * q.y - q.z * q.z;

	if (abs(y - x) < 1e-4f)
		return 0.0f;
	return atan2(y, x);
}

float &quat::operator[](int idx)
{
	switch (idx)
	{
		case 0:
			return this->x;
		case 1:
			return this->y;
		case 2:
			return this->z;
		case 3:
			return this->w;
		default:
			return this->x;
	}
}

float quat::operator[](int idx) const
{
	switch (idx)
	{
		case 0:
			return this->x;
		case 1:
			return this->y;
		case 2:
			return this->z;
		case 3:
			return this->w;
		default:
			return this->x;
	}
}

vec3 eulerAngles(const quat &q)
{
	return vec3(getPitch(q), getYaw(q), getRoll(q));
}

quat slerp(const quat &x, const quat &y, float a)
{
	float cosTheta = dot(x, y);

	quat z = y;
	if (cosTheta < 0.0f)
	{
		z = -1 * y;
		cosTheta = -cosTheta;
	}

	const float THRESHOLD = 0.9995f;
	if (cosTheta > THRESHOLD)
	{
		return quat(mix(x.w, z.w, a), mix(x.x, z.x, a), mix(x.y, z.y, a), mix(x.z, z.z, a));
	}
	else
	{
		float angle = acos(cosTheta);
		return (sin((1.0f - a) * angle) * x + sin(a * angle) * z) / sin(angle);
	}
}

} // namespace alglm
