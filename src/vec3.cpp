#include "include/math/vec3.h"
#include "include/math/vec4.h"

namespace alglm
{
// member function

vec3::vec3() : x(0.0f), y(0.0f), z(0.0f) {};

vec3::vec3(float x) : x(x), y(x), z(x) {};

vec3::vec3(float x, float y, float z) : x(x), y(y), z(z) {};

vec3::vec3(const vec3 &copy) : x(copy.x), y(copy.y), z(copy.z) {};

vec3::vec3(const vec4 &copy) : x(copy.x), y(copy.y), z(copy.z) {};

vec3 &vec3::operator=(const vec3 &copy)
{
	this->x = copy.x;
	this->y = copy.y;
	this->z = copy.z;
	return *this;
}

vec3 &vec3::operator=(float copy)
{
	this->x = copy;
	this->y = copy;
	this->z = copy;
	return *this;
}

vec3 vec3::operator+(const vec3 &rhs) const
{
	return vec3(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z);
}

vec3 vec3::operator-(const vec3 &rhs) const
{
	return vec3(this->x - rhs.x, this->y - rhs.y, this->z - rhs.z);
}

vec3 vec3::operator*(const vec3 &rhs) const
{
	return vec3(this->x * rhs.x, this->y * rhs.y, this->z * rhs.z);
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
		throw std::out_of_range("Index out of range");
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
		throw std::out_of_range("Index out of range");
	}
}


// non-member function

vec3 operator*(float scalar, const vec3 &vector)
{
	return vec3(scalar * vector.x, scalar * vector.y, scalar * vector.z);
}

vec3 operator*(const vec3 &vector, float scalar)
{
	return vec3(scalar * vector.x, scalar * vector.y, scalar * vector.z);
}

float dot(const vec3 &vector1, const vec3 &vector2)
{
	float ret = 0.0f;
	for (int i = 0; i < 3; i++)
	{
		ret += vector1[i] * vector2[i];
	}
	return ret;
}

float length(const vec3 &vector)
{
	return std::sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

float *value_ptr(vec3 &vector)
{
	return &vector[0];
}

vec3 cross(const vec3 &vector1, const vec3 &vector2)
{
	return vec3(vector1.y * vector2.z - vector1.z * vector2.y, vector1.z * vector2.x - vector1.x * vector2.z,
				vector1.x * vector2.y - vector1.y * vector2.x);
}

vec3 normalize(const vec3 &vector)
{
	float magnitude = 0;
	for (int i = 0; i < 3; i++)
	{
		magnitude += pow(vector[i], 2);
	}
	magnitude = sqrt(magnitude);
	return vec3(vector.x / magnitude, vector.y / magnitude, vector.z / magnitude);
}
} // namespace alglm