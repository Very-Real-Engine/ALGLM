#pragma once

namespace alglm
{

class quat;
class vec3;
class vec4;

class mat4
{
  public:
	vec4 data[4];

	mat4();
	mat4(float x);
	mat4(const vec4 &c1, const vec4 &c2, const vec4 &c3, const vec4 &c4);
	mat4(const mat4 &copy);
	mat4 &operator=(const mat4 &copy);
	mat4 operator+(const mat4 &rhs) const;
	mat4 operator-(const mat4 &rhs) const;
	mat4 operator*(const mat4 &rhs) const;
	vec4 &operator[](int idx);
	vec4 operator[](int idx) const;
};

mat4 operator*(float scalar, const mat4 &matrix);
mat4 operator*(const mat4 &matrix, float scalar);

float *value_ptr(mat4 &matrix);

mat4 inverse(const mat4 &matrix);
mat4 scale(const mat4 &matrix, const vec3 &vector);
mat4 translate(const mat4 &matrix, const vec3 &vector);
mat4 rotate(const mat4 &matrix, float theta, const vec3 &vector);
mat4 perspective(float fovy, float aspect, float zNear, float zFar);
mat4 lookAt(const vec3 &cameraPos, const vec3 &cameraTarget, const vec3 &cameraUp);
mat4 mat4_cast(const quat &quat);



} // namespace alglm