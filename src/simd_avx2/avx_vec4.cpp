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
	this->x = copy.x;
	this->y = copy.y;
	this->z = copy.z;
	this->w = copy.w;
	return *this;
}

float length(const vec4 &vector)
{
	return sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z + vector.w * vector.w);
}

float length2(const vec4 &vector)
{
	return vector.x * vector.x + vector.y * vector.y + vector.z * vector.z + vector.w * vector.w;
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
	return vec4(scalar * vector.x, scalar * vector.y, scalar * vector.z, scalar * vector.w);
}

vec4 operator*(const vec4 &vector, float scalar)
{
	return vec4(vector.x * scalar, vector.y * scalar, vector.z * scalar, vector.w * scalar);
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
	return (vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z + vector1.w * vector2.w);
}

vec4 normalize(const vec4 &vector)
{
	float magnitude = length(vector);
	return vec4(vector.x / magnitude, vector.y / magnitude, vector.z / magnitude, vector.w / magnitude);
}

vec4 mix(const vec4 &x, const vec4 &y, float a)
{
	return x * (1.0f - a) + y * a;
}

// ivec4

ivec4::ivec4() : x(0.0f), y(0.0f), z(0.0f), w(0.0f) {};

ivec4::ivec4(int32_t x) : x(x), y(x), z(x), w(x) {};

ivec4::ivec4(int32_t x, int32_t y, int32_t z, int32_t w) : x(x), y(y), z(z), w(w) {};

ivec4::ivec4(const ivec4 &copy) : x(copy.x), y(copy.y), z(copy.z), w(copy.w) {};

ivec4 &ivec4::operator=(const ivec4 &copy)
{
	this->x = copy.x;
	this->y = copy.y;
	this->z = copy.z;
	this->w = copy.w;
	return *this;
}

ivec4 ivec4::operator+(const ivec4 &rhs) const
{
	return ivec4(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z, this->w + rhs.w);
}

ivec4 ivec4::operator-(const ivec4 &rhs) const
{
	return ivec4(this->x - rhs.x, this->y - rhs.y, this->z - rhs.z, this->w - rhs.w);
}

ivec4 ivec4::operator*(const ivec4 &rhs) const
{
	return ivec4(this->x * rhs.x, this->y * rhs.y, this->z * rhs.z, this->w * rhs.w);
}

} // namespace alglm
