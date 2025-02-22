#include "alglm.h"

namespace alglm
{

vec3 radians(const vec3 &v)
{
	return vec3(radians(v.x), radians(v.y), radians(v.z));
}

vec3 degrees(const vec3 &v)
{
	return vec3(degrees(v.x), degrees(v.y), degrees(v.z));
}

} // namespace alglm
