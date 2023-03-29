#pragma once

#include "Particle2D.h"

class RPICParticle2D : public Particle2D
{
public:
	~RPICParticle2D() = default;

	RPICParticle2D(const glm::dvec2& inPosition = glm::dvec2(0.0), const glm::dvec2& inVelocity = glm::dvec2(0.0), double inMass = 1.0);

	const glm::dvec3& GetAngularVelocity() { return mAngularVelocity; }
	void SetAngularVelocity(const glm::dvec3& inAngularVelocity) { mAngularVelocity = inAngularVelocity; }

private:

	glm::dvec3 mAngularVelocity;
};
