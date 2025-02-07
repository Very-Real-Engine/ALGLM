#include "include/math/mat4.h"
#include "include/math/quat.h"
#include "include/math/vec3.h"
#include "include/math/vec4.h"

namespace alglm
{

// member function

mat4::mat4() : data{} {};

mat4::mat4(float x)
	: data{vec4(x, 0.0f, 0.0f, 0.0f), vec4(0.0f, x, 0.0f, 0.0f), vec4(0.0f, 0.0f, x, 0.0f), vec4(0.0f, 0.0f, 0.0f, x)} {
	  };

mat4::mat4(const vec4 &c1, const vec4 &c2, const vec4 &c3, const vec4 &c4) : data{c1, c2, c3, c4} {};

mat4::mat4(const mat4 &copy)
{
	for (int i = 0; i < 4; i++)
	{
		(*this)[i] = copy[i];
	}
}

mat4 &mat4::operator=(const mat4 &copy)
{
	for (int i = 0; i < 4; i++)
	{
		(*this)[i] = copy[i];
	}
	return *this;
}

mat4 mat4::operator+(const mat4 &rhs) const
{
	return mat4(vec4((*this)[0] + rhs[0]), vec4((*this)[1] + rhs[1]), vec4((*this)[2] + rhs[2]),
				vec4((*this)[3] + rhs[3]));
}

mat4 mat4::operator-(const mat4 &rhs) const
{
	return mat4(vec4((*this)[0] - rhs[0]), vec4((*this)[1] - rhs[1]), vec4((*this)[2] - rhs[2]),
				vec4((*this)[3] - rhs[3]));
}

mat4 mat4::operator*(const mat4 &rhs) const
{
	mat4 ret(0.0f);
	for (int i = 0; i < 4; i++)
	{
		ret[i] = (*this) * rhs[i];
	}
	return ret;
}

vec4 &mat4::operator[](int idx)
{
	if (0 <= idx && idx < 4)
		return this->data[idx];
	throw std::out_of_range("Index out of range");
}

vec4 mat4::operator[](int idx) const
{
	if (0 <= idx && idx < 4)
		return this->data[idx];
	throw std::out_of_range("Index out of range");
}

// non-member function

mat4 operator*(float scalar, const mat4 &matrix)
{
	return mat4({scalar * matrix[0], scalar * matrix[1], scalar * matrix[2], scalar * matrix[3]});
}

mat4 operator*(const mat4 &matrix, float scalar)
{
	return mat4({scalar * matrix[0], scalar * matrix[1], scalar * matrix[2], scalar * matrix[3]});
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

	vec4 dot0(matrix[0] * row0);
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

mat4 mat4_cast(const quat &quat)
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

} // namespace alglm
