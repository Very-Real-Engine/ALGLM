#include "alglm.h"

namespace alglm
{
// member function

vec2::vec2() : x(0.0f), y(0.0f) {};

vec2::vec2(float x) : x(x), y(x) {};

vec2::vec2(float x, float y) : x(x), y(y) {};

vec2::vec2(const vec2 &copy) : x(copy.x), y(copy.y) {};

vec2 &vec2::operator=(const vec2 &copy)
{
	this->x = copy.x;
	this->y = copy.y;
	return *this;
}

vec2 &vec2::operator=(float copy)
{
	this->x = copy;
	this->y = copy;
	return *this;
}

vec2 vec2::operator+(const vec2 &rhs) const
{
	return vec2(this->x + rhs.x, this->y + rhs.y);
}

vec2 vec2::operator-(const vec2 &rhs) const
{
	return vec2(this->x - rhs.x, this->y - rhs.y);
}

vec2 vec2::operator*(const vec2 &rhs) const
{
	return vec2(this->x * rhs.x, this->y * rhs.y);
}

float &vec2::operator[](int idx)
{
	switch (idx)
	{
	case 0:
		return this->x;
	case 1:
		return this->y;
	default:
		return this->x;
	}
}

float vec2::operator[](int idx) const
{
	switch (idx)
	{
	case 0:
		return this->x;
	case 1:
		return this->y;
	default:
		return this->x;
	}
}


// non-member function

vec2 operator*(float scalar, const vec2 &vector)
{
	return vec2(scalar * vector.x, scalar * vector.y);
}

vec2 operator*(const vec2 &vector, float scalar)
{
	return vec2(scalar * vector.x, scalar * vector.y);
}

float dot(const vec2 &vector1, const vec2 &vector2)
{
	return vector1.x * vector2.x + vector1.y * vector2.y;
}

float length(const vec2 &vector)
{
	return std::sqrt(vector.x * vector.x + vector.y * vector.y);
}

float length2(const vec2 &vector)
{
	return vector.x * vector.x + vector.y * vector.y;
}

float *value_ptr(vec2 &vector)
{
	return &vector[0];
}

vec2 normalize(const vec2 &vector)
{
	float magnitude = length(vector);
	return vec2(vector.x / magnitude, vector.y / magnitude);
}
} // namespace alglm