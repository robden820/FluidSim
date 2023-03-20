#pragma once

#include "Particle2D.h"

class APICParticle2D : public Particle2D
{
	public:

		APICParticle2D() = default;
		~APICParticle2D() = default;

		APICParticle2D(const glm::dvec2& inPosition, const glm::dvec2& inVelocity = glm::dvec2(0.0), double inMass = 1.0);

		const glm::dvec2& GetAffineState() { return mAffineState; }
		void SetAffineState(const glm::dvec2& inAff) { mAffineState = inAff; }

	private:
		glm::dvec2 mAffineState;
};