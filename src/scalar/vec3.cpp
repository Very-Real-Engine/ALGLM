#include "../../include/alglm.h"

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
		return this->x;
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
		return this->x;
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

vec3 operator/(float scalar, const vec3 &vector)
{
	return vec3(scalar / vector.x, scalar / vector.y, scalar / vector.z);
}

vec3 operator/(const vec3 &vector, float scalar)
{
	return vec3(vector.x / scalar, vector.y / scalar, vector.z / scalar);
}

float dot(const vec3 &vector1, const vec3 &vector2)
{
	return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
}

float length(const vec3 &vector)
{
	return std::sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

float length2(const vec3 &vector)
{
	return vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;
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
	float magnitude = length(vector);
	return vec3(vector.x / magnitude, vector.y / magnitude, vector.z / magnitude);
}

vec3 min(const vec3 &v1, const vec3 &v2)
{
	vec3 ret;

	ret.x = v1.x <= v2.x ? v1.x : v2.x;
	ret.y = v1.y <= v2.y ? v1.y : v2.y;
	ret.z = v1.z <= v2.z ? v1.z : v2.z;

	return ret;
}

vec3 max(const vec3 &v1, const vec3 &v2)
{
	vec3 ret;

	ret.x = v1.x >= v2.x ? v1.x : v2.x;
	ret.y = v1.y >= v2.y ? v1.y : v2.y;
	ret.z = v1.z >= v2.z ? v1.z : v2.z;

	return ret;
}

vec3 &vec3::operator-=(const vec3 &rhs)
{
	x -= rhs.x;
	y -= rhs.y;
	z -= rhs.z;
    return *this;
}

vec3 &vec3::operator+=(const vec3 &rhs)
{
	x += rhs.x;
	y += rhs.y;
	z += rhs.z;
    return *this;
}

vec3 &vec3::operator*=(const vec3 &rhs)
{
	x *= rhs.x;
	y *= rhs.y;
	z *= rhs.z;
    return *this;
}

vec3 vec3::operator-() const
{
	return vec3(-x, -y, -z);
}


} // namespace alglm