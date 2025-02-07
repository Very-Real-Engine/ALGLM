#pragma once

#include <cmath>
#include <stdexcept>

namespace alglm
{

class vec3;

class quat
{
  public:
	float x;
	float y;
	float z;
	float w;

	quat();
	quat(float x);
	quat(float x, float y, float z, float w);
	quat(const vec3 &axis, float angle);
	quat(const vec3 &eulerAngle);

	quat operator*(const quat &rhs) const;
};

} // namespace alglm