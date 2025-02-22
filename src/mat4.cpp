#include "alglm.h"

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
	(*this)[0][0] = 1.0f - 2.0f * quat.y * quat.y - 2.0f * quat.z * quat.z;
	(*this)[1][0] = 2.0f * quat.x * quat.y - 2.0f * quat.w * quat.z;
	(*this)[2][0] = 2.0f * quat.x * quat.z + 2.0f * quat.w * quat.y;
	(*this)[3][0] = 0.0f;

	(*this)[0][1] = 2.0f * quat.x * quat.y + 2.0f * quat.w * quat.z;
	(*this)[1][1] = 1.0f - 2.0f * quat.x * quat.x - 2.0f * quat.z * quat.z;
	(*this)[2][1] = 2.0f * quat.y * quat.z - 2.0f * quat.w * quat.x;
	(*this)[3][1] = 0.0f;

	(*this)[0][2] = 2.0f * quat.x * quat.z - 2.0f * quat.w * quat.y;
	(*this)[1][2] = 2.0f * quat.y * quat.z + 2.0f * quat.w * quat.x;
	(*this)[2][2] = 1.0f - 2.0f * quat.x * quat.x - 2.0f * quat.y * quat.y;
	(*this)[3][2] = 0.0f;

	(*this)[0][3] = 0.0f;
	(*this)[1][3] = 0.0f;
	(*this)[2][3] = 0.0f;
	(*this)[3][3] = 1.0f;
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

	for (int32_t i = 0; i < 16; ++i)
	{
		ret.data[i] = data[i] * rhs.data[i];
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
	float coef00 = matrix[2][2] * matrix[3][3] - matrix[3][2] * matrix[2][3];
	float coef02 = matrix[1][2] * matrix[3][3] - matrix[3][2] * matrix[1][3];
	float coef03 = matrix[1][2] * matrix[2][3] - matrix[2][2] * matrix[1][3];

	float coef04 = matrix[2][1] * matrix[3][3] - matrix[3][1] * matrix[2][3];
	float coef06 = matrix[1][1] * matrix[3][3] - matrix[3][1] * matrix[1][3];
	float coef07 = matrix[1][1] * matrix[2][3] - matrix[2][1] * matrix[1][3];

	float coef08 = matrix[2][1] * matrix[3][2] - matrix[3][1] * matrix[2][2];
	float coef10 = matrix[1][1] * matrix[3][2] - matrix[3][1] * matrix[1][2];
	float coef11 = matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2];

	float coef12 = matrix[2][0] * matrix[3][3] - matrix[3][0] * matrix[2][3];
	float coef14 = matrix[1][0] * matrix[3][3] - matrix[3][0] * matrix[1][3];
	float coef15 = matrix[1][0] * matrix[2][3] - matrix[2][0] * matrix[1][3];

	float coef16 = matrix[2][0] * matrix[3][2] - matrix[3][0] * matrix[2][2];
	float coef18 = matrix[1][0] * matrix[3][2] - matrix[3][0] * matrix[1][2];
	float coef19 = matrix[1][0] * matrix[2][2] - matrix[2][0] * matrix[1][2];

	float coef20 = matrix[2][0] * matrix[3][1] - matrix[3][0] * matrix[2][1];
	float coef22 = matrix[1][0] * matrix[3][1] - matrix[3][0] * matrix[1][1];
	float coef23 = matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1];

	vec4 fac0(coef00, coef00, coef02, coef03);
	vec4 fac1(coef04, coef04, coef06, coef07);
	vec4 fac2(coef08, coef08, coef10, coef11);
	vec4 fac3(coef12, coef12, coef14, coef15);
	vec4 fac4(coef16, coef16, coef18, coef19);
	vec4 fac5(coef20, coef20, coef22, coef23);

	vec4 vec0(matrix[1][0], matrix[0][0], matrix[0][0], matrix[0][0]);
	vec4 vec1(matrix[1][1], matrix[0][1], matrix[0][1], matrix[0][1]);
	vec4 vec2(matrix[1][2], matrix[0][2], matrix[0][2], matrix[0][2]);
	vec4 vec3(matrix[1][3], matrix[0][3], matrix[0][3], matrix[0][3]);

	vec4 inv0(vec1 * fac0 - vec2 * fac1 + vec3 * fac2);
	vec4 inv1(vec0 * fac0 - vec2 * fac3 + vec3 * fac4);
	vec4 inv2(vec0 * fac1 - vec1 * fac3 + vec3 * fac5);
	vec4 inv3(vec0 * fac2 - vec1 * fac4 + vec2 * fac5);

	vec4 signA(+1, -1, +1, -1);
	vec4 signB(-1, +1, -1, +1);
	mat4 inverse(inv0 * signA, inv1 * signB, inv2 * signA, inv3 * signB);

	vec4 row0(inverse[0][0], inverse[1][0], inverse[2][0], inverse[3][0]);

	vec4 dot0((*matrix[0]) * row0);
	float det = (dot0.x + dot0.y) + (dot0.z + dot0.w);

	float oneOverDeterminant = 1.0f / det;

	return inverse * oneOverDeterminant;
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

bool decompose(const mat4 &modelMatrix, vec3 &translation, quat &rotation, vec3 &scale, vec3 &skew, vec4 &perspective)
{
	mat4 localMatrix = modelMatrix;

	if (abs(localMatrix[3][3]) < std::numeric_limits<float>::epsilon())
	{
		return false; // 분해 실패
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			localMatrix[i][j] /= localMatrix[3][3];
		}
	}

	translation.x = localMatrix[3][0];
	translation.y = localMatrix[3][1];
	translation.z = localMatrix[3][2];

	perspective.x = localMatrix[3][0];
	perspective.y = localMatrix[3][1];
	perspective.z = localMatrix[3][2];
	perspective.w = localMatrix[3][3];

	localMatrix[3][0] = 0.0f;
	localMatrix[3][1] = 0.0f;
	localMatrix[3][2] = 0.0f;
	localMatrix[3][3] = 1.0f;

	for (int i = 0; i < 3; i++)
	{
		vec3 tmp(localMatrix[i][0], localMatrix[i][1], localMatrix[i][2]);
		scale[i] = length(tmp);
		localMatrix[i][0] = tmp.x / scale[i];
		localMatrix[i][1] = tmp.y / scale[i];
		localMatrix[i][2] = tmp.z / scale[i];
		localMatrix[i][3] = 0.0f;
	}

	vec3 col1(localMatrix[0][0], localMatrix[0][1], localMatrix[0][2]);
	vec3 col2(localMatrix[1][0], localMatrix[1][1], localMatrix[1][2]);
	vec3 col3(localMatrix[2][0], localMatrix[2][1], localMatrix[2][2]);
	skew.x = dot(col1, col2);
	skew.y = dot(col1, col3);
	skew.z = dot(col2, col3);

	rotation = quat_cast(localMatrix);

	return true;
}

quat quat_cast(const mat3 &matrix)
{
	float trace = matrix[0][0] + matrix[1][1] + matrix[2][2];
	quat q;

	if (trace > 0.0f)
	{
		float s = 0.5f / std::sqrt(trace + 1.0f);
		q.w = 0.25f / s;
		q.x = (matrix[2][1] - matrix[1][2]) * s;
		q.y = (matrix[0][2] - matrix[2][0]) * s;
		q.z = (matrix[1][0] - matrix[0][1]) * s;
	}
	else
	{
		if (matrix[0][0] > matrix[1][1] && matrix[0][0] > matrix[2][2])
		{
			float s = 2.0f * std::sqrt(1.0f + matrix[0][0] - matrix[1][1] - matrix[2][2]);
			q.w = (matrix[2][1] - matrix[1][2]) / s;
			q.x = 0.25f * s;
			q.y = (matrix[0][1] + matrix[1][0]) / s;
			q.z = (matrix[0][2] + matrix[2][0]) / s;
		}
		else if (matrix[1][1] > matrix[2][2])
		{
			float s = 2.0f * std::sqrt(1.0f + matrix[1][1] - matrix[0][0] - matrix[2][2]);
			q.w = (matrix[0][2] - matrix[2][0]) / s;
			q.x = (matrix[0][1] + matrix[1][0]) / s;
			q.y = 0.25f * s;
			q.z = (matrix[1][2] + matrix[2][1]) / s;
		}
		else
		{
			float s = 2.0f * std::sqrt(1.0f + matrix[2][2] - matrix[0][0] - matrix[1][1]);
			q.w = (matrix[1][0] - matrix[0][1]) / s;
			q.x = (matrix[0][2] + matrix[2][0]) / s;
			q.y = (matrix[1][2] + matrix[2][1]) / s;
			q.z = 0.25f * s;
		}
	}

	return normalize(q);
}

vec3 eulerAngles(const quat &q)
{
	vec3 angles;

	float sinrCosp = 2.0f * (q.w * q.x + q.y * q.z);
	float cosrCosp = 1.0f - 2.0f * (q.x * q.x + q.y * q.y);
	angles.x = std::atan2(sinrCosp, cosrCosp);

	float sinp = 2.0f * (q.w * q.y - q.z * q.x);
	if (std::abs(sinp) >= 1.0f)
	{
		angles.y = std::copysign(pi<float>() / 2.0f, sinp);
	}
	else
	{
		angles.y = std::asin(sinp);
	}

	// Yaw (z-axis rotation)
	float sinyCosp = 2.0f * (q.w * q.z + q.x * q.y);
	float cosyCosp = 1.0f - 2.0f * (q.y * q.y + q.z * q.z);
	angles.z = std::atan2(sinyCosp, cosyCosp);

	return angles;
}

} // namespace alglm
