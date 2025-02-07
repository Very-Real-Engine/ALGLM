#include "include/math/quat.h"
#include "include/math/math_utils.h"
#include "include/math/vec3.h"

namespace alglm
{
// member function

quat::quat() : x(0), y(0), z(0), w(1) {};

quat::quat(float x) : x(x), y(x), z(x), w(x) {};

quat::quat(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {};

quat::quat(const vec3 &axis, float angle)
{
	vec3 normalAxis = normalize(axis);
	float halfAngle = angle / 2.0f;
	float s = sin(halfAngle);
	w = cos(halfAngle);
	x = normalAxis.x * s;
	y = normalAxis.y * s;
	z = normalAxis.z * s;
}

quat::quat(const vec3 &eulerAngle)
{
	quat xQuat(vec3(1.0f, 0.0f, 0.0f), radians(eulerAngle.x));
	quat yQuat(vec3(0.0f, 1.0f, 0.0f), radians(eulerAngle.y));
	quat zQuat(vec3(0.0f, 0.0f, 1.0f), radians(eulerAngle.z));
	quat mulQuat = xQuat * yQuat * zQuat;

	x = mulQuat.x;
	y = mulQuat.y;
	z = mulQuat.z;
	w = mulQuat.w;
}

quat quat::operator*(const quat &rhs) const
{
	return quat(rhs.w * x + rhs.x * w + rhs.y * z - rhs.z * y, rhs.w * y - rhs.x * z + rhs.y * w + rhs.z * x,
				rhs.w * z + rhs.x * y - rhs.y * x + rhs.z * w, rhs.w * w - rhs.x * x - rhs.y * y - rhs.z * z);
}

} // namespace alglm
