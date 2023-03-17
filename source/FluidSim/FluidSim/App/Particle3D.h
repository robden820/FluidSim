#pragma once

#include "Particle.h"
#include "glm/vec3.hpp"

class Particle3D : public Particle
{
public:
	Particle3D() = default;
	~Particle3D() = default;

	Particle3D(const glm::dvec3& inPosition, const glm::dvec3& inVelocity = glm::dvec3(0.0, 0.0, 0.0), double inMass = 1.0);

	void StepParticle(double deltaTime);

	const glm::dvec3& GetPosition() const { return mPosition; }
	void SetPosition(glm::dvec3 inPos) { mPosition = inPos; }

	const glm::dvec3& GetVelocity() const { return mVelocity; }
	void SetVelocity(const glm::dvec3& inVel) { mVelocity = inVel; }

	const glm::dvec3& GetAcceleration() const { return mAcceleration; }
	void SetAcceleration(const glm::dvec3& inAcc) { mAcceleration = inAcc; }

	void ApplyForce(const glm::dvec3& inForce) { mForceAccumulator += inForce; }

private:
	glm::dvec3 mPosition;
	glm::dvec3 mVelocity;
	glm::dvec3 mAcceleration;

	glm::dvec3 mForceAccumulator;
};

