#pragma once

#include "glm\glm.hpp"
#include "Particle.h"

class Particle2D : public Particle
{
public:
	Particle2D() = default;
	~Particle2D() = default;

	Particle2D(const glm::vec2& inPosition, const glm::vec2& inVelocity = glm::vec2(0.0f, 0.0f), float inMass = 1.0f, float inRadius = 0.1f);

	void StepRK3(float deltaTime, const glm::vec2& K1, const glm::vec2& K2, const glm::vec2& K3);

	const glm::vec2& GetPosition() const { return mPosition; }
	void SetPosition(const glm::vec2& inPos) { mPosition = inPos; }

	const glm::vec2& GetVelocity() const { return mVelocity; }
	void SetVelocity(const glm::vec2& inVel) { mVelocity = inVel; }
	void SetXVelocity(float inXVel) { mVelocity.x = inXVel; }
	void SetYVelocity(float inYVel) { mVelocity.y = inYVel; }

	const glm::vec2& GetAcceleration() const { return mAcceleration; }
	void SetAcceleration(const glm::vec2& inAcc) { mAcceleration = inAcc; }

	void ApplyForce(const glm::vec2& inForce) { mForceAccumulator += inForce; }

private:
	glm::vec2 mPosition;
	glm::vec2 mVelocity;
	glm::vec2 mAcceleration;

	glm::vec2 mForceAccumulator; // For applying non-gravitational forces.
};