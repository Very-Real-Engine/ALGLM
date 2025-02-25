#include "../../include/alglm.h"

namespace alglm
{
// member function

vec2::vec2() : x(0.0f), y(0.0f), z(0.0f), w(0.0f) {};

vec2::vec2(float x) : x(x), y(x), z(0.0f), w(0.0f){};

vec2::vec2(float x, float y) : x(x), y(y), z(0.0f), w(0.0f) {};

vec2::vec2(const vec2 &copy) : x(copy.x), y(copy.y), z(0.0f), w(0.0f) {};

vec2 &vec2::operator=(const vec2 &copy)
{
    __m128 a = _mm_load_ps(reinterpret_cast<const float*>(&copy));
    _mm_store_ps(reinterpret_cast<float*>(this), a);
    return *this;
}

vec2 vec2::operator+(const vec2 &rhs) const
{
	// this와 rhs를 128비트 레지스터로 로드
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 r = _mm_add_ps(a, b);
	vec2 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec2 vec2::operator-(const vec2 &rhs) const
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 r = _mm_sub_ps(a, b);
	vec2 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec2 vec2::operator*(const vec2 &rhs) const
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 r = _mm_mul_ps(a, b);
	vec2 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

float &vec2::operator[](int idx)
{
	switch (idx)
	{
	case 0:
		return this->x;
	case 1:
		return this->y;
	default:
		return this->x;
	}
}

float vec2::operator[](int idx) const
{
	switch (idx)
	{
	case 0:
		return this->x;
	case 1:
		return this->y;
	default:
		return this->x;
	}
}

// non-member function

vec2 operator*(float scalar, const vec2 &vector)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector));
	__m128 s = _mm_set1_ps(scalar);
	__m128 r = _mm_mul_ps(a, s);
	vec2 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec2 operator*(const vec2 &vector, float scalar)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector));
	__m128 s = _mm_set1_ps(scalar);
	__m128 r = _mm_mul_ps(a, s);
	vec2 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

float dot(const vec2 &vector1, const vec2 &vector2)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector1));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&vector2));
	__m128 mul = _mm_mul_ps(a, b);
	__m128 sum = _mm_hadd_ps(mul, mul);
	return _mm_cvtss_f32(sum);
}

float length(const vec2 &vector)
{
	float dotVal = dot(vector, vector);
	__m128 d = _mm_set_ss(dotVal);
	__m128 sqrtVal = _mm_sqrt_ss(d);
	return _mm_cvtss_f32(sqrtVal);
}

float length2(const vec2 &vector)
{
	return dot(vector, vector);
}

float *value_ptr(vec2 &vector)
{
	return &vector[0];
}

vec2 normalize(const vec2 &vector)
{
	float mag = length(vector);
	// 브로드캐스트 후 나눗셈
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector));
	__m128 m = _mm_set1_ps(mag);
	__m128 r = _mm_div_ps(a, m);
	vec2 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}
} // namespace alglm