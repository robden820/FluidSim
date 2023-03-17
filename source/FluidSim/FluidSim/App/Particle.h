#pragma once

#include "glm/vec3.hpp"

class Particle
{
	public:
		virtual ~Particle() = default;

		double GetMass() const { return mMass; };

	protected:

		double mMass;
};