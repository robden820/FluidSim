#pragma once

#include "Particle2D.h"

class APICParticle2D : public Particle2D
{
	public:

		APICParticle2D() = default;
		~APICParticle2D() = default;

		APICParticle2D(const glm::dvec2& inPosition, const glm::dvec2& inVelocity = glm::dvec2(0.0, 0.0), double inAngularMomentum = 0.0, double inMass = 1.0);

		const double GetAngularMomentum() const { return mAngularMomentum; }
		void SetAngularMomentum(const double inAngMomentum) { mAngularMomentum = inAngMomentum; }

	private:
		double mAngularMomentum;
};