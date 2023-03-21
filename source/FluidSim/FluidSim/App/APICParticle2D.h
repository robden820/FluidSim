#pragma once

#include "Particle2D.h"

class APICParticle2D : public Particle2D
{
	public:

		APICParticle2D() = default;
		~APICParticle2D() = default;

		APICParticle2D(const glm::dvec2& inPosition, const glm::dvec2& inVelocity = glm::dvec2(0.0), double inMass = 1.0);

		const glm::dmat2& GetAffineState() { return mAffineState; }
		void SetAffineState(const glm::dmat2& inAffineState) { mAffineState = inAffineState; }

		const glm::dmat2& GetInertiaAnalog() { return mInertiaAnalog; }
		void SetInertiaAnalog(const glm::dmat2& inInertiaAnalog) { mInertiaAnalog = inInertiaAnalog; }

		bool HasValidInertiaAnalog() { return mValidInertiaAnalog; }
		void SetValidInertiaAnalog(bool isValid) { mValidInertiaAnalog = isValid; }

	private:

		glm::dmat2 mAffineState;
		glm::dmat2 mInertiaAnalog;

		bool mValidInertiaAnalog;
};