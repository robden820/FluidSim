#pragma once

#include "glm/vec3.hpp"

class Particle
{
	public:
		virtual ~Particle() = default;

		float GetMass() const { return mMass; };
		float GetRadius() const { return mRadius; };

	protected:

		float mMass;
		float mRadius;
};