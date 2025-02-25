#include "../../include/alglm.h"

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
	__m128 q1 = _mm_load_ps(&this->x); // (x1, y1, z1, w1)
	__m128 q2 = _mm_load_ps(&rhs.x);   // (x2, y2, z2, w2)

	// Compute term1: q1.w * q2
	__m128 q1_w = _mm_shuffle_ps(q1, q1, _MM_SHUFFLE(3, 3, 3, 3)); // (w1, w1, w1, w1)
	__m128 term1 = _mm_mul_ps(q1_w, q2);						   // (w1*x2, w1*y2, w1*z2, w1*w2)

	// Compute term2: q2.w * q1
	__m128 q2_w = _mm_shuffle_ps(q2, q2, _MM_SHUFFLE(3, 3, 3, 3)); // (w2, w2, w2, w2)
	__m128 term2 = _mm_mul_ps(q2_w, q1);						   // (w2*x1, w2*y1, w2*z1, w2*w1)

	// Compute cross product for the vector part:
	// cross(q1.xyz, q2.xyz) = (q1.y*q2.z - q1.z*q2.y,
	//                           q1.z*q2.x - q1.x*q2.z,
	//                           q1.x*q2.y - q1.y*q2.x)
	__m128 q1_yzx = _mm_shuffle_ps(q1, q1, _MM_SHUFFLE(3, 0, 2, 1)); // (y1, z1, x1, w1)
	__m128 q2_zxy = _mm_shuffle_ps(q2, q2, _MM_SHUFFLE(3, 1, 0, 2)); // (z2, x2, y2, w2)
	__m128 q1_zxy = _mm_shuffle_ps(q1, q1, _MM_SHUFFLE(3, 1, 0, 2)); // (z1, x1, y1, w1)
	__m128 q2_yzx = _mm_shuffle_ps(q2, q2, _MM_SHUFFLE(3, 0, 2, 1)); // (y2, z2, x2, w2)
	__m128 cross = _mm_sub_ps(_mm_mul_ps(q1_yzx, q2_zxy), _mm_mul_ps(q1_zxy, q2_yzx));

	// Compute xyz part: term1_xyz + term2_xyz + cross.
	__m128 xyz = _mm_add_ps(_mm_add_ps(term1, term2), cross);

	// For the scalar part, compute: w1*w2 - dot(q1.xyz, q2.xyz)
	// Compute dot(q1.xyz, q2.xyz) using _mm_dp_ps with mask 0x71.
	__m128 dp = _mm_dp_ps(q1, q2, 0x71);
	float dot_xyz = _mm_cvtss_f32(dp);
	float w1 = this->w, w2 = rhs.w;
	float res_w = w1 * w2 - dot_xyz;

	// Store xyz into an array.
	float result_array[4];
	_mm_store_ps(result_array, xyz);
	// result_array[0] = x, [1] = y, [2] = z.
	// Construct quaternion using constructor: quat(w, x, y, z).
	return quat(res_w, result_array[0], result_array[1], result_array[2]);
}

quat quat::operator+(const quat &rhs) const
{
	__m128 a = _mm_load_ps(&this->x);
	__m128 b = _mm_load_ps(&rhs.x);
	__m128 r = _mm_add_ps(a, b);
	quat ret;
	_mm_store_ps(&ret.x, r);
	return ret;
}

quat quat::operator-(const quat &rhs) const
{
	__m128 a = _mm_load_ps(&this->x);
	__m128 b = _mm_load_ps(&rhs.x);
	__m128 r = _mm_sub_ps(a, b);
	quat ret;
	_mm_store_ps(&ret.x, r);
	return ret;
}

quat normalize(const quat &q)
{
	__m128 a = _mm_load_ps(&q.x);
	float len = length(q);
	__m128 m = _mm_set1_ps(len);
	__m128 r = _mm_div_ps(a, m);
	quat ret;
	_mm_store_ps(&ret.x, r);
	return ret;
}

float dot(const quat &q1, const quat &q2)
{
	__m128 a = _mm_load_ps(&q1.x);
	__m128 b = _mm_load_ps(&q2.x);
	__m128 dp = _mm_dp_ps(a, b, 0xFF); // sum all four components
	return _mm_cvtss_f32(dp);
}

float length(const quat &q)
{
	float dp = dot(q, q);
	__m128 d = _mm_set_ss(dp);
	__m128 sqrtVal = _mm_sqrt_ss(d);
	return _mm_cvtss_f32(sqrtVal);
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
	__m128 v = _mm_load_ps(&q.x);
	__m128 s = _mm_set1_ps(a);
	__m128 r = _mm_mul_ps(v, s);
	quat ret;
	_mm_store_ps(&ret.x, r);
	return ret;
}

quat operator*(const quat &q, float a)
{
	__m128 v = _mm_load_ps(&q.x);
	__m128 s = _mm_set1_ps(a);
	__m128 r = _mm_mul_ps(v, s);
	quat ret;
	_mm_store_ps(&ret.x, r);
	return ret;
}

quat operator/(float a, const quat &q)
{
	__m128 v = _mm_load_ps(&q.x);
	__m128 s = _mm_set1_ps(a);
	__m128 r = _mm_div_ps(v, s);
	quat ret;
	_mm_store_ps(&ret.x, r);
	return ret;
}

quat operator/(const quat &q, float a)
{
	__m128 v = _mm_load_ps(&q.x);
	__m128 s = _mm_set1_ps(a);
	__m128 r = _mm_div_ps(v, s);
	quat ret;
	_mm_store_ps(&ret.x, r);
	return ret;
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

quat &quat::operator+=(const quat &rhs)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&this->x));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs.x));
	__m128 result = _mm_add_ps(a, b);
	_mm_store_ps(reinterpret_cast<float *>(&this->x), result);
	return *this;
}

} // namespace alglm
