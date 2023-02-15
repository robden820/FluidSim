#pragma once

#include "glm/vec3.hpp"

class Particle
{
	public:
		virtual ~Particle() = default;

		double GetMass() const { return mMass; };
		double GetRadius() const { return mRadius; };

	protected:

		double mMass;
		double mRadius;
};