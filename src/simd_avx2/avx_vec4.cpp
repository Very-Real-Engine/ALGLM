#include "../../include/alglm.h"

namespace alglm
{
// member function

vec4::vec4() : x(0.0f), y(0.0f), z(0.0f), w(0.0f) {};

vec4::vec4(float x) : x(x), y(x), z(x), w(x) {};

vec4::vec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {};

vec4::vec4(const vec3 &copy, float w) : x(copy.x), y(copy.y), z(copy.z), w(w) {};

vec4::vec4(const vec4 &copy) : x(copy.x), y(copy.y), z(copy.z), w(copy.w) {};

vec4 &vec4::operator=(const vec4 &copy)
{
	// SIMD 방식으로 4개의 float를 한 번에 복사
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&copy));
	_mm_store_ps(reinterpret_cast<float *>(this), a);
	return *this;
}

float length(const vec4 &vector)
{
	float d = dot(vector, vector);
	__m128 dVal = _mm_set_ss(d);
	__m128 sqrtVal = _mm_sqrt_ss(dVal);
	return _mm_cvtss_f32(sqrtVal);
}

float length2(const vec4 &vector)
{
	return dot(vector, vector);
}

vec4 vec4::operator+(const vec4 &rhs) const
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 r = _mm_add_ps(a, b);
	vec4 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec4 vec4::operator-(const vec4 &rhs) const
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 r = _mm_sub_ps(a, b);
	vec4 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec4 vec4::operator*(const vec4 &rhs) const
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(this));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&rhs));
	__m128 r = _mm_mul_ps(a, b);
	vec4 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

float &vec4::operator[](int idx)
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

float vec4::operator[](int idx) const
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

// non-member function

vec4 operator*(float scalar, const vec4 &vector)
{
	__m128 v = _mm_load_ps(reinterpret_cast<const float *>(&vector));
	__m128 s = _mm_set1_ps(scalar);
	__m128 r = _mm_mul_ps(v, s);
	vec4 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec4 operator*(const vec4 &vector, float scalar)
{
	__m128 v = _mm_load_ps(reinterpret_cast<const float *>(&vector));
	__m128 s = _mm_set1_ps(scalar);
	__m128 r = _mm_mul_ps(v, s);
	vec4 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec4 operator*(const mat4 &matrix, const vec4 &vector)
{
	// 가정: mat4는 16바이트 정렬된 16 float가 연속된 data[]로 저장됨.
	// mat4가 column-major로 저장되었다고 가정하면,
	// ret = col0 * vector.x + col1 * vector.y + col2 * vector.z + col3 * vector.w;
	__m128 v = _mm_load_ps(&vector.x);
	__m128 col0 = _mm_load_ps(&matrix.data[0]);
	__m128 col1 = _mm_load_ps(&matrix.data[4]);
	__m128 col2 = _mm_load_ps(&matrix.data[8]);
	__m128 col3 = _mm_load_ps(&matrix.data[12]);

	__m128 vx = _mm_set1_ps(vector.x);
	__m128 vy = _mm_set1_ps(vector.y);
	__m128 vz = _mm_set1_ps(vector.z);
	__m128 vw = _mm_set1_ps(vector.w);

	__m128 res = _mm_add_ps(_mm_add_ps(_mm_mul_ps(col0, vx), _mm_mul_ps(col1, vy)),
							_mm_add_ps(_mm_mul_ps(col2, vz), _mm_mul_ps(col3, vw)));
	vec4 ret;
	_mm_store_ps(&ret.x, res);
	return ret;
}

float *value_ptr(vec4 &vector)
{
	return &vector[0];
}

float dot(const vec4 &vector1, const vec4 &vector2)
{
	__m128 a = _mm_load_ps(reinterpret_cast<const float *>(&vector1));
	__m128 b = _mm_load_ps(reinterpret_cast<const float *>(&vector2));
	// _mm_dp_ps (SSE4.1) computes dot product; mask 0xFF sums all 4 elements.
	__m128 dp = _mm_dp_ps(a, b, 0xFF);
	return _mm_cvtss_f32(dp);
}

vec4 normalize(const vec4 &vector)
{
	float mag = length(vector);
	__m128 v = _mm_load_ps(reinterpret_cast<const float *>(&vector));
	__m128 m = _mm_set1_ps(mag);
	__m128 r = _mm_div_ps(v, m);
	vec4 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

vec4 mix(const vec4 &x, const vec4 &y, float a)
{
	// mix(x,y,a) = x*(1-a) + y*a;
	__m128 vx = _mm_load_ps(reinterpret_cast<const float *>(&x));
	__m128 vy = _mm_load_ps(reinterpret_cast<const float *>(&y));
	__m128 factor = _mm_set1_ps(a);
	__m128 invFactor = _mm_set1_ps(1.0f - a);
	__m128 r = _mm_add_ps(_mm_mul_ps(vx, invFactor), _mm_mul_ps(vy, factor));
	vec4 result;
	_mm_store_ps(reinterpret_cast<float *>(&result), r);
	return result;
}

// ivec4

ivec4::ivec4() : x(0.0f), y(0.0f), z(0.0f), w(0.0f) {};

ivec4::ivec4(int32_t x) : x(x), y(x), z(x), w(x) {};

ivec4::ivec4(int32_t x, int32_t y, int32_t z, int32_t w) : x(x), y(y), z(z), w(w) {};

ivec4::ivec4(const ivec4 &copy) : x(copy.x), y(copy.y), z(copy.z), w(copy.w) {};

ivec4 &ivec4::operator=(const ivec4 &copy)
{
	// 정수형의 경우, _mm_load_si128 / _mm_store_si128를 사용합니다.
	__m128i a = _mm_load_si128(reinterpret_cast<const __m128i *>(&copy));
	_mm_store_si128(reinterpret_cast<__m128i *>(this), a);
	return *this;
}

ivec4 ivec4::operator+(const ivec4 &rhs) const
{
	__m128i a = _mm_load_si128(reinterpret_cast<const __m128i *>(this));
	__m128i b = _mm_load_si128(reinterpret_cast<const __m128i *>(&rhs));
	__m128i r = _mm_add_epi32(a, b);
	ivec4 result;
	_mm_store_si128(reinterpret_cast<__m128i *>(&result), r);
	return result;
}

ivec4 ivec4::operator-(const ivec4 &rhs) const
{
	__m128i a = _mm_load_si128(reinterpret_cast<const __m128i *>(this));
	__m128i b = _mm_load_si128(reinterpret_cast<const __m128i *>(&rhs));
	__m128i r = _mm_sub_epi32(a, b);
	ivec4 result;
	_mm_store_si128(reinterpret_cast<__m128i *>(&result), r);
	return result;
}

ivec4 ivec4::operator*(const ivec4 &rhs) const
{
	// _mm_mullo_epi32 multiplies 32-bit integers
	__m128i a = _mm_load_si128(reinterpret_cast<const __m128i *>(this));
	__m128i b = _mm_load_si128(reinterpret_cast<const __m128i *>(&rhs));
	__m128i r = _mm_mullo_epi32(a, b);
	ivec4 result;
	_mm_store_si128(reinterpret_cast<__m128i *>(&result), r);
	return result;
}

} // namespace alglm
