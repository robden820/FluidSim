#pragma once

#include "glm\glm.hpp"
#include "Particle.h"

class Particle2D : public Particle
{
public:
	Particle2D() = default;
	~Particle2D() = default;

	Particle2D(const glm::dvec2& inPosition, const glm::dvec2& inVelocity = glm::dvec2(0.0, 0.0), double inMass = 1.0, double inRadius = 0.1);

	void StepRK3(double deltaTime, const glm::dvec2& K1, const glm::dvec2& K2, const glm::dvec2& K3);

	const glm::dvec2& GetPosition() const { return mPosition; }
	void SetPosition(const glm::dvec2& inPos) { mPosition = inPos; }

	const glm::dvec2& GetVelocity() const { return mVelocity; }
	void SetVelocity(const glm::dvec2& inVel) { mVelocity = inVel; }
	void SetXVelocity(double inXVel) { mVelocity.x = inXVel; }
	void SetYVelocity(double inYVel) { mVelocity.y = inYVel; }

	const glm::dvec2& GetPreviousVelocity() const { return mPrevVelocity; }
	void SetPreviousVelocity() { mPrevVelocity = mVelocity; }

	const glm::dvec2& GetAcceleration() const { return mAcceleration; }
	void SetAcceleration(const glm::dvec2& inAcc) { mAcceleration = inAcc; }

	void ApplyForce(const glm::dvec2& inForce) { mForceAccumulator += inForce; }

private:
	glm::dvec2 mPosition;
	glm::dvec2 mVelocity;
	glm::dvec2 mPrevVelocity;
	glm::dvec2 mAcceleration;

	glm::dvec2 mForceAccumulator; // For applying non-gravitational forces.
};