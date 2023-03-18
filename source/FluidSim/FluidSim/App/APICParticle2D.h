#pragma once

#include "Particle2D.h"

class APICParticle2D : public Particle2D
{
	public:

		APICParticle2D() = default;
		~APICParticle2D() = default;

		APICParticle2D(const glm::dvec2& inPosition, const glm::dvec2& inVelocity = glm::dvec2(0.0, 0.0), const glm::dvec3& inAngularVelocity = glm::dvec3(0.0, 0.0, 0.0), double inMass = 1.0);

		const glm::dvec3& GetAngularVelocity() const { return mAngularVelocity; }
		void SetAngularVelocity(const glm::dvec3& inAngVel) { mAngularVelocity = inAngVel; }

	private:
		// Vec 3 even in 2D fluid as angular momentum will point in the z direction.
		glm::dvec3 mAngularVelocity;
};