#pragma once

#include "glm/vec3.hpp"

class Particle
{
	public:
		virtual ~Particle() = default;

		virtual void StepParticle(float deltaTime) = 0;

		virtual float GetMass() const { return mMass; }
		virtual float GetRadius() const{ return mRadius; }

	private:

		float mMass;
		float mRadius;
};