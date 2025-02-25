#include "../../include/alglm.h"

namespace alglm
{

// member function

mat4::mat4()
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
	data[12] = 0.0f;
	data[13] = 0.0f;
	data[14] = 0.0f;
	data[15] = 0.0f;
}

mat4::mat4(float x)
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
	data[12] = 0.0f;
	data[13] = 0.0f;
	data[14] = 0.0f;
	data[15] = x;
}

mat4::mat4(const vec4 &v1, const vec4 &v2, const vec4 &v3, const vec4 &v4)
{
	for (int32_t i = 0; i < 4; ++i)
	{
		data[i] = v1[i];
		data[4 + i] = v2[i];
		data[8 + i] = v3[i];
		data[12 + i] = v4[i];
	}
}

mat4::mat4(float f0, float f1, float f2, float f3, float f4, float f5, float f6, float f7, float f8, float f9,
		   float f10, float f11, float f12, float f13, float f14, float f15)
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
	data[9] = f9;
	data[10] = f10;
	data[11] = f11;
	data[12] = f12;
	data[13] = f13;
	data[14] = f14;
	data[15] = f15;
}

mat4::mat4(const mat4 &copy)
{
	for (int32_t i = 0; i < 16; i++)
	{
		data[i] = copy.data[i];
	}
}

mat4::mat4(const quat &quat)
{
	data[0] = 1.0f - 2.0f * quat.y * quat.y - 2.0f * quat.z * quat.z;
	data[4] = 2.0f * quat.x * quat.y - 2.0f * quat.w * quat.z;
	data[8] = 2.0f * quat.x * quat.z + 2.0f * quat.w * quat.y;
	data[12] = 0.0f;

	data[1] = 2.0f * quat.x * quat.y + 2.0f * quat.w * quat.z;
	data[5] = 1.0f - 2.0f * quat.x * quat.x - 2.0f * quat.z * quat.z;
	data[9] = 2.0f * quat.y * quat.z - 2.0f * quat.w * quat.x;
	data[13] = 0.0f;

	data[2] = 2.0f * quat.x * quat.z - 2.0f * quat.w * quat.y;
	data[6] = 2.0f * quat.y * quat.z + 2.0f * quat.w * quat.x;
	data[10] = 1.0f - 2.0f * quat.x * quat.x - 2.0f * quat.y * quat.y;
	data[14] = 0.0f;

	data[3] = 0.0f;
	data[7] = 0.0f;
	data[11] = 0.0f;
	data[15] = 1.0f;
}

mat4 &mat4::operator=(const mat4 &copy)
{
	for (int32_t i = 0; i < 16; i++)
	{
		data[i] = copy.data[i];
	}
	return *this;
}

mat4 mat4::operator+(const mat4 &rhs) const
{
	mat4 ret;

	for (int32_t i = 0; i < 16; ++i)
	{
		ret.data[i] = data[i] + rhs.data[i];
	}

	return ret;
}

mat4 mat4::operator-(const mat4 &rhs) const
{
	mat4 ret;

	for (int32_t i = 0; i < 16; ++i)
	{
		ret.data[i] = data[i] - rhs.data[i];
	}

	return ret;
}

mat4 mat4::operator*(const mat4 &rhs) const
{
	mat4 ret;

	ret[0][0] = data[0] * rhs.data[0] + data[4] * rhs.data[1] + data[8] * rhs.data[2] + data[12] * rhs.data[3];
	ret[0][1] = data[1] * rhs.data[0] + data[5] * rhs.data[1] + data[9] * rhs.data[2] + data[13] * rhs.data[3];
	ret[0][2] = data[2] * rhs.data[0] + data[6] * rhs.data[1] + data[10] * rhs.data[2] + data[14] * rhs.data[3];
	ret[0][3] = data[3] * rhs.data[0] + data[7] * rhs.data[1] + data[11] * rhs.data[2] + data[15] * rhs.data[3];

	ret[1][0] = data[0] * rhs.data[4] + data[4] * rhs.data[5] + data[8] * rhs.data[6] + data[12] * rhs.data[7];
	ret[1][1] = data[1] * rhs.data[4] + data[5] * rhs.data[5] + data[9] * rhs.data[6] + data[13] * rhs.data[7];
	ret[1][2] = data[2] * rhs.data[4] + data[6] * rhs.data[5] + data[10] * rhs.data[6] + data[14] * rhs.data[7];
	ret[1][3] = data[3] * rhs.data[4] + data[7] * rhs.data[5] + data[11] * rhs.data[6] + data[15] * rhs.data[7];

	ret[2][0] = data[0] * rhs.data[8] + data[4] * rhs.data[9] + data[8] * rhs.data[10] + data[12] * rhs.data[11];
	ret[2][1] = data[1] * rhs.data[8] + data[5] * rhs.data[9] + data[9] * rhs.data[10] + data[13] * rhs.data[11];
	ret[2][2] = data[2] * rhs.data[8] + data[6] * rhs.data[9] + data[10] * rhs.data[10] + data[14] * rhs.data[11];
	ret[2][3] = data[3] * rhs.data[8] + data[7] * rhs.data[9] + data[11] * rhs.data[10] + data[15] * rhs.data[11];

	ret[3][0] = data[0] * rhs.data[12] + data[4] * rhs.data[13] + data[8] * rhs.data[14] + data[12] * rhs.data[15];
	ret[3][1] = data[1] * rhs.data[12] + data[5] * rhs.data[13] + data[9] * rhs.data[14] + data[13] * rhs.data[15];
	ret[3][2] = data[2] * rhs.data[12] + data[6] * rhs.data[13] + data[10] * rhs.data[14] + data[14] * rhs.data[15];
	ret[3][3] = data[3] * rhs.data[12] + data[7] * rhs.data[13] + data[11] * rhs.data[14] + data[15] * rhs.data[15];

	return ret;
}

float *mat4::operator[](int idx)
{
	return &data[idx * 4];
}

const float *mat4::operator[](int idx) const
{
	return &data[idx * 4];
}

// non-member function

mat4 operator*(float scalar, const mat4 &matrix)
{
	mat4 ret;

	for (int32_t i = 0; i < 16; ++i)
	{
		ret.data[i] = scalar * matrix.data[i];
	}

	return ret;
}

mat4 operator*(const mat4 &matrix, float scalar)
{
	mat4 ret;

	for (int32_t i = 0; i < 16; ++i)
	{
		ret.data[i] = scalar * matrix.data[i];
	}

	return ret;
}

float *value_ptr(mat4 &matrix)
{
	return &matrix[0][0];
}

mat4 inverse(const mat4 &matrix)
{
	mat4 inv;
	float det;
	int i;

	inv[0][0] = matrix[1][1] * matrix[2][2] * matrix[3][3] - matrix[1][1] * matrix[2][3] * matrix[3][2] -
				matrix[2][1] * matrix[1][2] * matrix[3][3] + matrix[2][1] * matrix[1][3] * matrix[3][2] +
				matrix[3][1] * matrix[1][2] * matrix[2][3] - matrix[3][1] * matrix[1][3] * matrix[2][2];

	inv[0][1] = -matrix[0][1] * matrix[2][2] * matrix[3][3] + matrix[0][1] * matrix[2][3] * matrix[3][2] +
				matrix[2][1] * matrix[0][2] * matrix[3][3] - matrix[2][1] * matrix[0][3] * matrix[3][2] -
				matrix[3][1] * matrix[0][2] * matrix[2][3] + matrix[3][1] * matrix[0][3] * matrix[2][2];

	inv[0][2] = matrix[0][1] * matrix[1][2] * matrix[3][3] - matrix[0][1] * matrix[1][3] * matrix[3][2] -
				matrix[1][1] * matrix[0][2] * matrix[3][3] + matrix[1][1] * matrix[0][3] * matrix[3][2] +
				matrix[3][1] * matrix[0][2] * matrix[1][3] - matrix[3][1] * matrix[0][3] * matrix[1][2];

	inv[0][3] = -matrix[0][1] * matrix[1][2] * matrix[2][3] + matrix[0][1] * matrix[1][3] * matrix[2][2] +
				matrix[1][1] * matrix[0][2] * matrix[2][3] - matrix[1][1] * matrix[0][3] * matrix[2][2] -
				matrix[2][1] * matrix[0][2] * matrix[1][3] + matrix[2][1] * matrix[0][3] * matrix[1][2];

	inv[1][0] = -matrix[1][0] * matrix[2][2] * matrix[3][3] + matrix[1][0] * matrix[2][3] * matrix[3][2] +
				matrix[2][0] * matrix[1][2] * matrix[3][3] - matrix[2][0] * matrix[1][3] * matrix[3][2] -
				matrix[3][0] * matrix[1][2] * matrix[2][3] + matrix[3][0] * matrix[1][3] * matrix[2][2];

	inv[1][1] = matrix[0][0] * matrix[2][2] * matrix[3][3] - matrix[0][0] * matrix[2][3] * matrix[3][2] -
				matrix[2][0] * matrix[0][2] * matrix[3][3] + matrix[2][0] * matrix[0][3] * matrix[3][2] +
				matrix[3][0] * matrix[0][2] * matrix[2][3] - matrix[3][0] * matrix[0][3] * matrix[2][2];

	inv[1][2] = -matrix[0][0] * matrix[1][2] * matrix[3][3] + matrix[0][0] * matrix[1][3] * matrix[3][2] +
				matrix[1][0] * matrix[0][2] * matrix[3][3] - matrix[1][0] * matrix[0][3] * matrix[3][2] -
				matrix[3][0] * matrix[0][2] * matrix[1][3] + matrix[3][0] * matrix[0][3] * matrix[1][2];

	inv[1][3] = matrix[0][0] * matrix[1][2] * matrix[2][3] - matrix[0][0] * matrix[1][3] * matrix[2][2] -
				matrix[1][0] * matrix[0][2] * matrix[2][3] + matrix[1][0] * matrix[0][3] * matrix[2][2] +
				matrix[2][0] * matrix[0][2] * matrix[1][3] - matrix[2][0] * matrix[0][3] * matrix[1][2];

	inv[2][0] = matrix[1][0] * matrix[2][1] * matrix[3][3] - matrix[1][0] * matrix[2][3] * matrix[3][1] -
				matrix[2][0] * matrix[1][1] * matrix[3][3] + matrix[2][0] * matrix[1][3] * matrix[3][1] +
				matrix[3][0] * matrix[1][1] * matrix[2][3] - matrix[3][0] * matrix[1][3] * matrix[2][1];

	inv[2][1] = -matrix[0][0] * matrix[2][1] * matrix[3][3] + matrix[0][0] * matrix[2][3] * matrix[3][1] +
				matrix[2][0] * matrix[0][1] * matrix[3][3] - matrix[2][0] * matrix[0][3] * matrix[3][1] -
				matrix[3][0] * matrix[0][1] * matrix[2][3] + matrix[3][0] * matrix[0][3] * matrix[2][1];

	inv[2][2] = matrix[0][0] * matrix[1][1] * matrix[3][3] - matrix[0][0] * matrix[1][3] * matrix[3][1] -
				matrix[1][0] * matrix[0][1] * matrix[3][3] + matrix[1][0] * matrix[0][3] * matrix[3][1] +
				matrix[3][0] * matrix[0][1] * matrix[1][3] - matrix[3][0] * matrix[0][3] * matrix[1][1];

	inv[2][3] = -matrix[0][0] * matrix[1][1] * matrix[2][3] + matrix[0][0] * matrix[1][3] * matrix[2][1] +
				matrix[1][0] * matrix[0][1] * matrix[2][3] - matrix[1][0] * matrix[0][3] * matrix[2][1] -
				matrix[2][0] * matrix[0][1] * matrix[1][3] + matrix[2][0] * matrix[0][3] * matrix[1][1];

	inv[3][0] = -matrix[1][0] * matrix[2][1] * matrix[3][2] + matrix[1][0] * matrix[2][2] * matrix[3][1] +
				matrix[2][0] * matrix[1][1] * matrix[3][2] - matrix[2][0] * matrix[1][2] * matrix[3][1] -
				matrix[3][0] * matrix[1][1] * matrix[2][2] + matrix[3][0] * matrix[1][2] * matrix[2][1];

	inv[3][1] = matrix[0][0] * matrix[2][1] * matrix[3][2] - matrix[0][0] * matrix[2][2] * matrix[3][1] -
				matrix[2][0] * matrix[0][1] * matrix[3][2] + matrix[2][0] * matrix[0][2] * matrix[3][1] +
				matrix[3][0] * matrix[0][1] * matrix[2][2] - matrix[3][0] * matrix[0][2] * matrix[2][1];

	inv[3][2] = -matrix[0][0] * matrix[1][1] * matrix[3][2] + matrix[0][0] * matrix[1][2] * matrix[3][1] +
				matrix[1][0] * matrix[0][1] * matrix[3][2] - matrix[1][0] * matrix[0][2] * matrix[3][1] -
				matrix[3][0] * matrix[0][1] * matrix[1][2] + matrix[3][0] * matrix[0][2] * matrix[1][1];

	inv[3][3] = matrix[0][0] * matrix[1][1] * matrix[2][2] - matrix[0][0] * matrix[1][2] * matrix[2][1] -
				matrix[1][0] * matrix[0][1] * matrix[2][2] + matrix[1][0] * matrix[0][2] * matrix[2][1] +
				matrix[2][0] * matrix[0][1] * matrix[1][2] - matrix[2][0] * matrix[0][2] * matrix[1][1];

	det = matrix[0][0] * inv[0][0] + matrix[0][1] * inv[1][0] + matrix[0][2] * inv[2][0] + matrix[0][3] * inv[3][0];

	det = 1.0f / det;
	for (i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			inv[i][j] *= det;
		}
	}

	return inv;
}

mat4 scale(const mat4 &matrix, const vec3 &vector)
{
	mat4 scaleMatrix(vec4(vector.x, 0.0f, 0.0f, 0.0f), vec4(0.0f, vector.y, 0.0f, 0.0f),
					 vec4(0.0f, 0.0f, vector.z, 0.0f), vec4(0.0f, 0.0f, 0.0f, 1.0f));
	return matrix * scaleMatrix;
}

mat4 translate(const mat4 &matrix, const vec3 &vector)
{
	mat4 translateMatrix(vec4(1.0f, 0.0f, 0.0f, 0.0f), vec4(0.0f, 1.0f, 0.0f, 0.0f), vec4(0.0f, 0.0f, 1.0f, 0.0f),
						 vec4(vector.x, vector.y, vector.z, 1.0f));
	return matrix * translateMatrix;
}

mat4 rotate(const mat4 &matrix, float theta, const vec3 &vector)
{
	vec3 a = normalize(vector);
	mat4 rotateMatrix(vec4(cos(theta) + (1 - cos(theta)) * a.x * a.x, (1 - cos(theta)) * a.y * a.x + sin(theta) * a.z,
						   (1 - cos(theta)) * a.z * a.x - sin(theta) * a.y, 0.0f),
					  vec4((1 - cos(theta)) * a.x * a.y - sin(theta) * a.z, cos(theta) + (1 - cos(theta)) * a.y * a.y,
						   (1 - cos(theta)) * a.z * a.y + sin(theta) * a.x, 0.0f),
					  vec4((1 - cos(theta)) * a.x * a.z + sin(theta) * a.y,
						   (1 - cos(theta)) * a.y * a.z - sin(theta) * a.x, cos(theta) + (1 - cos(theta)) * a.z * a.z,
						   0.0f),
					  vec4(0.0f, 0.0f, 0.0f, 1.0f));
	return matrix * rotateMatrix;
}

mat4 perspective(float fovy, float aspect, float zNear, float zFar)
{
	return mat4(vec4(1 / (aspect * tan(fovy / 2)), 0.0f, 0.0f, 0.0f), vec4(0.0f, 1 / tan(fovy / 2), 0.0f, 0.0f),
				vec4(0.0f, 0.0f, (zFar + zNear) / (zNear - zFar), -1.0f),
				vec4(0.0f, 0.0f, (2 * zFar * zNear) / (zNear - zFar), 0.0f));
}

mat4 lookAt(const vec3 &cameraPos, const vec3 &cameraTarget, const vec3 &cameraUp)
{
	vec3 cameraZ = normalize(cameraPos - cameraTarget);
	vec3 cameraX = normalize(cross(cameraUp, cameraZ));
	vec3 cameraY = normalize(cross(cameraZ, cameraX));
	mat4 matrix(vec4(cameraX.x, cameraX.y, cameraX.z, 0.0f), vec4(cameraY.x, cameraY.y, cameraY.z, 0.0f),
				vec4(cameraZ.x, cameraZ.y, cameraZ.z, 0.0f), vec4(cameraPos.x, cameraPos.y, cameraPos.z, 1.0f));
	return inverse(matrix);
}

mat4 toMat4(const quat &quat)
{
	mat4 rotate(1.0f);

	rotate[0][0] = 1.0f - 2.0f * quat.y * quat.y - 2.0f * quat.z * quat.z;
	rotate[1][0] = 2.0f * quat.x * quat.y - 2.0f * quat.w * quat.z;
	rotate[2][0] = 2.0f * quat.x * quat.z + 2.0f * quat.w * quat.y;

	rotate[0][1] = 2.0f * quat.x * quat.y + 2.0f * quat.w * quat.z;
	rotate[1][1] = 1.0f - 2.0f * quat.x * quat.x - 2.0f * quat.z * quat.z;
	rotate[2][1] = 2.0f * quat.y * quat.z - 2.0f * quat.w * quat.x;

	rotate[0][2] = 2.0f * quat.x * quat.z - 2.0f * quat.w * quat.y;
	rotate[1][2] = 2.0f * quat.y * quat.z + 2.0f * quat.w * quat.x;
	rotate[2][2] = 1.0f - 2.0f * quat.x * quat.x - 2.0f * quat.y * quat.y;

	return rotate;
}

mat4 ortho(float left, float right, float bottom, float top, float zNear, float zFar)
{
	mat4 result(1.0f);

	result[0][0] = 2.0f / (right - left);
	result[1][1] = 2.0f / (top - bottom);
	result[2][2] = -2.0f / (zFar - zNear);

	result[3][0] = -(right + left) / (right - left);
	result[3][1] = -(top + bottom) / (top - bottom);
	result[3][2] = -(zFar + zNear) / (zFar - zNear);
	result[3][3] = 1.0f;

	return result;
}

float determinant(const mat4 &m)
{
	float SubFactor00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
	float SubFactor01 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
	float SubFactor02 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
	float SubFactor03 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
	float SubFactor04 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
	float SubFactor05 = m[2][0] * m[3][1] - m[3][0] * m[2][1];

	vec4 detCof(+(m[1][1] * SubFactor00 - m[1][2] * SubFactor01 + m[1][3] * SubFactor02),
				-(m[1][0] * SubFactor00 - m[1][2] * SubFactor03 + m[1][3] * SubFactor04),
				+(m[1][0] * SubFactor01 - m[1][1] * SubFactor03 + m[1][3] * SubFactor05),
				-(m[1][0] * SubFactor02 - m[1][1] * SubFactor04 + m[1][2] * SubFactor05));

	return m[0][0] * detCof[0] + m[0][1] * detCof[1] + m[0][2] * detCof[2] + m[0][3] * detCof[3];
}

bool decompose(const mat4 &modelMatrix, vec3 &scale, quat &rotation, vec3 &translation, vec4 &perspective)
{
	mat4 localMatrix(modelMatrix);
	if (abs(localMatrix[3][3]) <= 1e-6f)
	{
		return false;
	}

	for (int32_t i = 0; i < 4; ++i)
	{
		for (int32_t j = 0; j < 4; ++j)
		{
			localMatrix[i][j] /= localMatrix[3][3];
		}
	}

	mat4 perspectiveMatrix(localMatrix);
	for (int32_t i = 0; i < 3; ++i)
	{
		perspectiveMatrix[i][3] = 0.0f;
	}
	perspectiveMatrix[3][3] = 1.0f;

	if (abs(determinant(perspectiveMatrix)) <= 1e-6f)
	{
		return false;
	}

	if (abs(localMatrix[0][3]) > 1e-6f || abs(localMatrix[1][3]) > 1e-6f || abs(localMatrix[2][3]) > 1e-6f)
	{
		vec4 rhs;
		rhs[0] = localMatrix[0][3];
		rhs[1] = localMatrix[1][3];
		rhs[2] = localMatrix[2][3];
		rhs[3] = localMatrix[3][3];

		mat4 inversePerspectiveMatrix = inverse(perspectiveMatrix);
		mat4 transposedInversePerspectiveMatrix = transpose(inversePerspectiveMatrix);
		
		perspective = transposedInversePerspectiveMatrix * rhs;

		localMatrix[0][3] = localMatrix[1][3] = localMatrix[2][3] = 0.0f;
		localMatrix[3][3] = 1.0f;
	}
	else
	{
		perspective = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	}

	translation = vec3(localMatrix[3][0], localMatrix[3][1], localMatrix[3][2]);
	localMatrix[3][0] = 0.0f;
	localMatrix[3][1] = 0.0f;
	localMatrix[3][2] = 0.0f;

	vec3 row[3], pDum3;

	for (int32_t i = 0; i < 3; ++i)
	{
		for (int32_t j = 0; j < 3; ++j)
		{
			row[i][j] = localMatrix[i][j];
		}
	}
	scale.x = length(row[0]);
	row[0] = row[0] / scale.x;

	scale.y = length(row[1]);
	row[1] = row[1] / scale.y;

	scale.z = length(row[2]);
	row[2] = row[2] / scale.z;

	pDum3 = cross(row[1], row[2]);
	if (dot(row[0], pDum3) < 0)
	{
		for (int32_t i = 0; i < 3; ++i)
		{
			scale[i] = scale[i] * -1.0f;
			row[i] = row[i] * -1.0f;
		}
	}

	int32_t i, j, k = 0;

	float root, trace = row[0].x + row[1].y + row[2].z;
	if (trace > 0.0f)
	{
		root = sqrt(trace + 1.0f);
		rotation.w = 0.5f * root;
		root = 0.5f / root;
		rotation.x = root * (row[1].z - row[2].y);
		rotation.y = root * (row[2].x - row[0].z);
		rotation.z = root * (row[0].y - row[1].x);
	}
	else
	{
		static int next[3] = {1, 2, 0};
		i = 0;
		if (row[1].y > row[0].x)
		{
			i = 1;
		}
		if (row[2].z > row[i][i])
		{
			i = 2;
		}
		j = next[i];
		k = next[j];

		root = sqrt(row[i][i] - row[j][j] - row[k][k] + 1.0f);

		rotation[i] = 0.5f * root;
		root = 0.5f / root;
		rotation[j] = root * (row[i][j] + row[j][i]);
		rotation[k] = root * (row[i][k] + row[k][i]);
		rotation.w = root * (row[j][k] - row[k][j]);
	}
	return true;
}

quat quat_cast(const mat3 &matrix)
{
	float x = matrix[0][0] - matrix[1][1] - matrix[2][2];
	float y = matrix[1][1] - matrix[0][0] - matrix[2][2];
	float z = matrix[2][2] - matrix[0][0] - matrix[1][1];
	float w = matrix[0][0] + matrix[1][1] + matrix[2][2];

	int bigIdx = 0;
	float v = w;
	if (x > v)
	{
		v = x;
		bigIdx = 1;
	}
	if (y > v)
	{
		v = y;
		bigIdx = 2;
	}
	if (z > v)
	{
		v = z;
		bigIdx = 3;
	}

	float bigValue = sqrt(v + 1.0f) * 0.5f;
	float mult = 0.25f / bigValue;

	switch(bigIdx)
	{
		case 0:
			return quat(bigValue, (matrix[1][2] - matrix[2][1]) * mult, (matrix[2][0] - matrix[0][2]) * mult, (matrix[0][1] - matrix[1][0]) * mult);
		case 1:
			return quat((matrix[1][2] - matrix[2][1]) * mult, bigValue, (matrix[0][1] + matrix[1][0]) * mult, (matrix[2][0] + matrix[0][2]) * mult);
		case 2:
			return quat((matrix[2][0] - matrix[0][2]) * mult, (matrix[0][1] + matrix[1][0]) * mult, bigValue, (matrix[1][2] + matrix[2][1]) * mult);
		case 3:
			return quat((matrix[0][1] - matrix[1][0]) * mult, (matrix[2][0] + matrix[0][2]) * mult, (matrix[1][2] + matrix[2][1]) * mult, bigValue);
		default:
			return quat(1, 0, 0, 0);
	}
}

mat4 transpose(const mat4 &matrix)
{
	return mat4(matrix[0][0], matrix[1][0], matrix[2][0], matrix[3][0], matrix[0][1], matrix[1][1], matrix[2][1],
				matrix[3][1], matrix[0][2], matrix[1][2], matrix[2][2], matrix[3][2], matrix[0][3], matrix[1][3],
				matrix[2][3], matrix[3][3]);
}
} // namespace alglm
