#include "../include/alglm.h"

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
}

mat3::mat3(float x)
{
	data[0] = x;
	data[1] = 0.0f;
	data[2] = 0.0f;
	data[3] = 0.0f;
	data[4] = x;
	data[5] = 0.0f;
	data[6] = 0.0f;
	data[7] = 0.0f;
	data[8] = x;
}

mat3::mat3(float f0, float f1, float f2, float f3, float f4, float f5, float f6, float f7, float f8)
{
	data[0] = f0;
	data[1] = f1;
	data[2] = f2;
	data[3] = f3;
	data[4] = f4;
	data[5] = f5;
	data[6] = f6;
	data[7] = f7;
	data[8] = f8;
}

mat3::mat3(const mat3 &copy)
{
	for (int32_t i = 0; i < 9; i++)
	{
		data[i] = copy.data[i];
	}
}

mat3::mat3(const mat4 &copy)
{
	for (int32_t i = 0; i < 9; i++)
	{
		for (int32_t j = 0; j < 9; j++)
		{
			data[i * 3 + j] = copy.data[i * 4 + j];
		}
	}
}

mat3::mat3(const vec3 &v1, const vec3 &v2, const vec3 &v3)
{
	for (int32_t i = 0; i < 3; ++i)
	{
		data[i] = v1[i];
		data[3 + i] = v2[i];
		data[6 + i] = v3[i];
	}
}

mat3 &mat3::operator=(const mat3 &copy)
{
	for (int32_t i = 0; i < 9; i++)
	{
		data[i] = copy.data[i];
	}
	return *this;
}

mat3 mat3::operator+(const mat3 &rhs) const
{
	mat3 ret;

	for (int32_t i = 0; i < 9; ++i)
	{
		ret.data[i] = data[i] + rhs.data[i];
	}

	return ret;
}

mat3 mat3::operator-(const mat3 &rhs) const
{
	mat3 ret;

	for (int32_t i = 0; i < 9; ++i)
	{
		ret.data[i] = data[i] - rhs.data[i];
	}

	return ret;
}

mat3 mat3::operator*(const mat3 &rhs) const
{
	mat3 ret;

	ret[0][0] = data[0] * rhs.data[0] + data[3] * rhs.data[1] + data[6] * rhs.data[2];
	ret[0][1] = data[1] * rhs.data[0] + data[4] * rhs.data[1] + data[7] * rhs.data[2];
	ret[0][2] = data[2] * rhs.data[0] + data[5] * rhs.data[1] + data[8] * rhs.data[2];

	ret[1][0] = data[0] * rhs.data[3] + data[3] * rhs.data[4] + data[6] * rhs.data[5];
	ret[1][1] = data[1] * rhs.data[3] + data[4] * rhs.data[4] + data[7] * rhs.data[5];
	ret[1][2] = data[2] * rhs.data[3] + data[5] * rhs.data[4] + data[8] * rhs.data[5];

	ret[2][0] = data[0] * rhs.data[6] + data[3] * rhs.data[7] + data[6] * rhs.data[8];
	ret[2][1] = data[1] * rhs.data[6] + data[4] * rhs.data[7] + data[7] * rhs.data[8];
	ret[2][2] = data[2] * rhs.data[6] + data[5] * rhs.data[7] + data[8] * rhs.data[8];

	return ret;
}

float *mat3::operator[](int idx)
{
	return &data[idx * 3];
}

const float *mat3::operator[](int idx) const
{
	return &data[idx * 3];
}

// non-member function

mat3 operator*(float scalar, const mat3 &matrix)
{
	mat3 ret;

	for (int32_t i = 0; i < 9; ++i)
	{
		ret.data[i] = scalar * matrix.data[i];
	}

	return ret;
}

mat3 operator*(const mat3 &matrix, float scalar)
{
	mat3 ret;

	for (int32_t i = 0; i < 9; ++i)
	{
		ret.data[i] = scalar * matrix.data[i];
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
	inv[0][0] = (matrix.data[4] * matrix.data[8] - matrix.data[5] * matrix.data[7]) * invDet;
	inv[0][1] = (matrix.data[2] * matrix.data[7] - matrix.data[1] * matrix.data[8]) * invDet;
	inv[0][2] = (matrix.data[1] * matrix.data[5] - matrix.data[2] * matrix.data[4]) * invDet;
	inv[1][0] = (matrix.data[5] * matrix.data[6] - matrix.data[3] * matrix.data[8]) * invDet;
	inv[1][1] = (matrix.data[0] * matrix.data[8] - matrix.data[2] * matrix.data[6]) * invDet;
	inv[1][2] = (matrix.data[2] * matrix.data[3] - matrix.data[0] * matrix.data[5]) * invDet;
	inv[2][0] = (matrix.data[3] * matrix.data[7] - matrix.data[4] * matrix.data[6]) * invDet;
	inv[2][1] = (matrix.data[1] * matrix.data[6] - matrix.data[0] * matrix.data[7]) * invDet;
	inv[2][2] = (matrix.data[0] * matrix.data[4] - matrix.data[1] * matrix.data[3]) * invDet;

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

} // namespace alglm
