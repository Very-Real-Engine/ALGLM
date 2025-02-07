#include "include/math/vec4.h"
#include "include/math/vec3.h"
#include "include/math/mat4.h"

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
	this->x = copy.x;
	this->y = copy.y;
	this->z = copy.z;
	this->w = copy.w;
	return *this;
}

vec4 vec4::operator+(const vec4 &rhs) const
{
	return vec4(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z, this->w + rhs.w);
}

vec4 vec4::operator-(const vec4 &rhs) const
{
	return vec4(this->x - rhs.x, this->y - rhs.y, this->z - rhs.z, this->w - rhs.w);
}

vec4 vec4::operator*(const vec4 &rhs) const
{
	return vec4(this->x * rhs.x, this->y * rhs.y, this->z * rhs.z, this->w * rhs.w);
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
		throw std::out_of_range("Index out of range");
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
		throw std::out_of_range("Index out of range");
	}
}

// non-member function

vec4 operator*(float scalar, const vec4 &vector)
{
	return vec4(scalar * vector.x, scalar * vector.y, scalar * vector.z, scalar * vector.w);
}

vec4 operator*(const vec4 &vector, float scalar)
{
	return vec4(scalar * vector.x, scalar * vector.y, scalar * vector.z, scalar * vector.w);
}

vec4 operator*(const mat4 &matrix, const vec4 &vector)
{
	vec4 ret(0.0f);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			ret[i] += matrix[j][i] * vector[j];
		}
	}
	return ret;
}

float *value_ptr(vec4 &vector)
{
	return &vector[0];
}

float dot(const vec4 &vector1, const vec4 &vector2)
{
	float ret = 0.0f;
	for (int i = 0; i < 4; i++)
	{
		ret += vector1[i] * vector2[i];
	}
	return ret;
}

vec4 normalize(const vec4 &vector)
{
	float magnitude = 0;
	for (int i = 0; i < 4; i++)
	{
		magnitude += pow(vector[i], 2);
	}
	magnitude = sqrt(magnitude);
	return vec4(vector.x / magnitude, vector.y / magnitude, vector.z / magnitude, vector.w / magnitude);
}

} // namespace alglm
