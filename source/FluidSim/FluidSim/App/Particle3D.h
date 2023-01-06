#pragma once

#include "Particle.h"
#include "glm/vec3.hpp"

class Particle3D : public Particle
{
public:
	Particle3D() = default;
	~Particle3D() = default;

	Particle3D(const glm::vec3& inPosition, const glm::vec3& inVelocity = glm::vec3(0.0f, 0.0f, 0.0f), float inMass = 1.0f, float inRadius = 0.1f);

	void StepParticle(float deltaTime) override;

	const glm::vec3& GetPosition() const { return mPosition; }
	void SetPosition(glm::vec3 inPos) { mPosition = inPos; }

	const glm::vec3& GetVelocity() const { return mVelocity; }
	void SetVelocity(const glm::vec3& inVel) { mVelocity = inVel; }

	const glm::vec3& GetAcceleration() const { return mAcceleration; }
	void SetAcceleration(const glm::vec3& inAcc) { mAcceleration = inAcc; }

	void ApplyForce(const glm::vec3& inForce) { mForceAccumulator += inForce; }

private:
	glm::vec3 mPosition;
	glm::vec3 mVelocity;
	glm::vec3 mAcceleration;

	glm::vec3 mForceAccumulator;
};

