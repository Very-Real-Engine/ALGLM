#include "../../include/alglm.h"

namespace alglm
{

// member function
mat3::mat3()
{
	data[0] = 0.0f;
	data[1] = 0.0f;
	data[2] = 0.0f;
	data[3] = 0.0f;
	data[4] = 0.0f;
	data[5] = 0.0f;
	data[6] = 0.0f;
	data[7] = 0.0f;
	data[8] = 0.0f;
	data[9] = 0.0f;
	data[10] = 0.0f;
	data[11] = 0.0f;
}

mat3::mat3(float x)
{
	data[0] = x;
	data[1] = 0.0f;
	data[2] = 0.0f;
	data[3] = 0.0f;
	data[4] = 0.0f;
	data[5] = x;
	data[6] = 0.0f;
	data[7] = 0.0f;
	data[8] = 0.0f;
	data[9] = 0.0f;
	data[10] = x;
	data[11] = 0.0f;
}

mat3::mat3(float f0, float f1, float f2, float f3, float f4, float f5, float f6, float f7, float f8)
{
	data[0] = f0;
	data[1] = f1;
	data[2] = f2;
	data[3] = 0.0f;
	data[4] = f3;
	data[5] = f4;
	data[6] = f5;
	data[7] = 0.0f;
	data[8] = f6;
	data[9] = f7;
	data[10] = f8;
	data[11] = 0.0f;
}

mat3::mat3(const mat3 &copy)
{
	for (int32_t i = 0; i < 12; i++)
	{
		data[i] = copy.data[i];
	}
}

mat3::mat3(const mat4 &copy)
{
	data[0] = copy.data[0];
	data[1] = copy.data[1];
	data[2] = copy.data[2];
	data[3] = 0.0f;
	data[4] = copy.data[4];
	data[5] = copy.data[5];
	data[6] = copy.data[6];
	data[7] = 0.0f;
	data[8] = copy.data[8];
	data[9] = copy.data[9];
	data[10] = copy.data[10];
	data[11] = 0.0f;
}

mat3::mat3(const vec3 &v1, const vec3 &v2, const vec3 &v3)
{
	for (int32_t i = 0; i < 3; ++i)
	{
		data[i] = v1[i];
		data[4 + i] = v2[i];
		data[8 + i] = v3[i];
	}
	data[3] = 0.0f;
	data[7] = 0.0f;
	data[11] = 0.0f;
}

mat3 &mat3::operator=(const mat3 &copy)
{
	for (int32_t i = 0; i < 12; i++)
	{
		data[i] = copy.data[i];
	}
	return *this;
}

mat3 mat3::operator+(const mat3 &rhs) const
{
	mat3 ret;
	for (int i = 0; i < 3; i++)
	{
		__m128 a = _mm_load_ps(&data[i * 4]);
		__m128 b = _mm_load_ps(&rhs.data[i * 4]);
		__m128 r = _mm_add_ps(a, b);
		_mm_store_ps(&ret.data[i * 4], r);
	}
	return ret;
}

mat3 mat3::operator-(const mat3 &rhs) const
{
	mat3 ret;
	for (int i = 0; i < 3; i++)
	{
		__m128 a = _mm_load_ps(&data[i * 4]);
		__m128 b = _mm_load_ps(&rhs.data[i * 4]);
		__m128 r = _mm_sub_ps(a, b);
		_mm_store_ps(&ret.data[i * 4], r);
	}
	return ret;
}

mat3 mat3::operator*(const mat3 &rhs) const
{
    mat3 ret;
    for (int r = 0; r < 3; r++) {
        __m128 lhs_row = _mm_set_ps(0.0f, data[r + 8], data[r + 4], data[r]);
        for (int c = 0; c < 3; c++) {
            __m128 rhs_col = _mm_load_ps(&rhs.data[c * 4]);
            __m128 dp = _mm_dp_ps(lhs_row, rhs_col, 0x71);
            ret.data[r + c * 4] = _mm_cvtss_f32(dp);
        }
    }
    for (int c = 0; c < 3; c++) {
        ret.data[3 + c * 4] = 0.0f;
    }
    return ret;
}

float *mat3::operator[](int idx)
{
	return &data[idx * 4];
}

const float *mat3::operator[](int idx) const
{
	return &data[idx * 4];
}

// non-member function

mat3 operator*(float scalar, const mat3 &matrix)
{
    mat3 ret;
    __m128 s = _mm_set1_ps(scalar);
    for (int i = 0; i < 3; i++) {
        __m128 a = _mm_load_ps(&matrix.data[i * 4]);
        __m128 r = _mm_mul_ps(a, s);
        _mm_store_ps(&ret.data[i * 4], r);
    }
    return ret;
}

mat3 operator*(const mat3 &matrix, float scalar)
{
    mat3 ret;
    __m128 s = _mm_set1_ps(scalar);
    for (int i = 0; i < 3; i++) {
        __m128 a = _mm_load_ps(&matrix.data[i * 4]);
        __m128 r = _mm_mul_ps(a, s);
        _mm_store_ps(&ret.data[i * 4], r);
    }
    return ret;
}

float *value_ptr(mat3 &matrix)
{
	return matrix[0];
}

mat3 inverse(const mat3 &matrix)
{
	float det = determinant(matrix);

	float invDet = 1.0f / det;

	mat3 inv;
	inv[0][0] = (matrix.data[5] * matrix.data[10] - matrix.data[6] * matrix.data[9]) * invDet;
	inv[0][1] = (matrix.data[2] * matrix.data[9] - matrix.data[1] * matrix.data[10]) * invDet;
	inv[0][2] = (matrix.data[1] * matrix.data[6] - matrix.data[2] * matrix.data[5]) * invDet;
	inv[1][0] = (matrix.data[6] * matrix.data[8] - matrix.data[4] * matrix.data[10]) * invDet;
	inv[1][1] = (matrix.data[0] * matrix.data[10] - matrix.data[2] * matrix.data[8]) * invDet;
	inv[1][2] = (matrix.data[2] * matrix.data[4] - matrix.data[0] * matrix.data[6]) * invDet;
	inv[2][0] = (matrix.data[4] * matrix.data[9] - matrix.data[5] * matrix.data[8]) * invDet;
	inv[2][1] = (matrix.data[1] * matrix.data[8] - matrix.data[0] * matrix.data[9]) * invDet;
	inv[2][2] = (matrix.data[0] * matrix.data[5] - matrix.data[1] * matrix.data[4]) * invDet;

	return inv;
}

mat3 transpose(const mat3 &matrix)
{
	return mat3(matrix[0][0], matrix[1][0], matrix[2][0], matrix[0][1], matrix[1][1], matrix[2][1], matrix[0][2],
				matrix[1][2], matrix[2][2]);
}

float determinant(const mat3 &m)
{
	return m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) - m[1][0] * (m[0][1] * m[2][2] - m[2][1] * m[0][2]) +
		   m[2][0] * (m[0][1] * m[1][2] - m[1][1] * m[0][2]);
}

vec3 operator*(const mat3 &matrix, const vec3 &v)
{
    __m128 x = _mm_set1_ps(v.x); 
    __m128 y = _mm_set1_ps(v.y);
    __m128 z = _mm_set1_ps(v.z);

    __m128 col0 = _mm_load_ps(matrix[0]);
    __m128 col1 = _mm_load_ps(matrix[1]);
    __m128 col2 = _mm_load_ps(matrix[2]);

    __m128 sum = _mm_add_ps(_mm_add_ps(_mm_mul_ps(x, col0), _mm_mul_ps(y, col1)), _mm_mul_ps(z, col2));

    float result[4];
    _mm_storeu_ps(result, sum);
    return vec3(result[0], result[1], result[2]);
}

} // namespace alglm
