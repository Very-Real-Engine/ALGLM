#include "../../include/alglm.h"

namespace alglm
{
// member function

vec3::vec3() : x(0.0f), y(0.0f), z(0.0f), w(0.0f) {};

vec3::vec3(float x) : x(x), y(x), z(x), w(0.0f) {};

vec3::vec3(float x, float y, float z) : x(x), y(y), z(z), w(0.0f) {};

vec3::vec3(const vec3 &copy) : x(copy.x), y(copy.y), z(copy.z), w(0.0f) {};

vec3::vec3(const vec4 &copy) : x(copy.x), y(copy.y), z(copy.z), w(0.0f) {};

vec3 &vec3::operator=(const vec3 &copy)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&copy));
	_mm_store_ps(reinterpret_cast<float *>(this), a);
	return *this;
}

vec3 &vec3::operator=(float copy)
{
	__m128 val = _mm_set1_ps(copy);
	_mm_store_ps(reinterpret_cast<float *>(this), val);
	return *this;
}

vec3 vec3::operator+(const vec3 &rhs) const
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 r = _mm_add_ps(a, b);
	vec3 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec3 vec3::operator-(const vec3 &rhs) const
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 r = _mm_sub_ps(a, b);
	vec3 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec3 vec3::operator*(const vec3 &rhs) const
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 r = _mm_mul_ps(a, b);
	vec3 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

float &vec3::operator[](int idx)
{
	switch (idx)
	{
	case 0:
		return this->x;
	case 1:
		return this->y;
	case 2:
		return this->z;
	default:
		return this->x;
	}
}

float vec3::operator[](int idx) const
{
	switch (idx)
	{
	case 0:
		return this->x;
	case 1:
		return this->y;
	case 2:
		return this->z;
	default:
		return this->x;
	}
}

// non-member function

vec3 operator*(float scalar, const vec3 &vector)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector));
	__m128 s = _mm_set1_ps(scalar);
	__m128 r = _mm_mul_ps(a, s);
	vec3 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec3 operator*(const vec3 &vector, float scalar)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector));
	__m128 s = _mm_set1_ps(scalar);
	__m128 r = _mm_mul_ps(a, s);
	vec3 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec3 operator/(float scalar, const vec3 &vector)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector));
	__m128 s = _mm_set1_ps(scalar);
	__m128 r = _mm_div_ps(a, s);
	vec3 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec3 operator/(const vec3 &vector, float scalar)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector));
	__m128 s = _mm_set1_ps(scalar);
	__m128 r = _mm_div_ps(a, s);
	vec3 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

float dot(const vec3 &vector1, const vec3 &vector2)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector1));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&vector2));
	__m128 dp = _mm_dp_ps(a, b, 0x71);
	return _mm_cvtss_f32(dp);
}

float length(const vec3 &vector)
{
	float d = dot(vector, vector);
	__m128 dVal = _mm_set_ss(d);
	__m128 sqrtVal = _mm_sqrt_ss(dVal);
	return _mm_cvtss_f32(sqrtVal);
}

float length2(const vec3 &vector)
{
	return dot(vector, vector);
}

float *value_ptr(vec3 &vector)
{
	return &vector[0];
}

vec3 cross(const vec3 &vector1, const vec3 &vector2)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector1));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&vector2));
	// Shuffle vectors: a_yzx, b_yzx
	__m128 a_yzx = _mm_shuffle_ps(a, a, _MM_SHUFFLE(3, 0, 2, 1));
	__m128 b_yzx = _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 0, 2, 1));
	__m128 c = _mm_sub_ps(_mm_mul_ps(a, b_yzx), _mm_mul_ps(a_yzx, b));
	__m128 result = _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 0, 2, 1));
	vec3 ret;
	_mm_store_ps(reinterpret_cast<float *>(&ret), result);
	return ret;
}

vec3 normalize(const vec3 &vector)
{
	float mag = length(vector);
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector));
	__m128 m = _mm_set1_ps(mag);
	__m128 r = _mm_div_ps(a, m);
	vec3 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec3 min(const vec3 &v1, const vec3 &v2)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&v1));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&v2));
	__m128 m = _mm_min_ps(a, b);
	vec3 ret;
	_mm_store_ps(reinterpret_cast<float *>(&ret), m);
	return ret;
}

vec3 max(const vec3 &v1, const vec3 &v2)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&v1));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&v2));
	__m128 m = _mm_max_ps(a, b);
	vec3 ret;
	_mm_store_ps(reinterpret_cast<float *>(&ret), m);
	return ret;
}

vec3 &vec3::operator+=(const vec3 &rhs)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 sum = _mm_add_ps(a, b);
	_mm_store_ps(reinterpret_cast<float *>(this), sum);
	return *this;
}

vec3 &vec3::operator-=(const vec3 &rhs)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 diff = _mm_sub_ps(a, b);
	_mm_store_ps(reinterpret_cast<float *>(this), diff);
	return *this;
}

vec3 &vec3::operator*=(const vec3 &rhs)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 result = _mm_mul_ps(a, b);
	_mm_store_ps(reinterpret_cast<float *>(this), result);
	return *this;
}

vec3 vec3::operator-() const
{
	__m128 v = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 neg = _mm_set1_ps(-1.0f);
	__m128 result = _mm_mul_ps(v, neg);
	vec3 ret;
	_mm_store_ps(reinterpret_cast<float *>(&ret), result);
	return ret;
}

} // namespace alglm