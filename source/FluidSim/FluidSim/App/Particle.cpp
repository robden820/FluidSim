#include "Particle.h"

Particle::Particle(glm::vec3 inPosition, glm::vec3 inVelocity, float inMass)
{
	mPosition = inPosition;
	mVelocity = inVelocity;
	mMass = inMass;

	float g = -9.8f * mMass; // Downwards force due to gravity;

	mForceAccumulator = { 0.0f, g, 0.0f }; 
}

void Particle::StepParticle(float deltaTime)
{
	glm::vec3 acceleration = mForceAccumulator / mMass;

	mPosition += mVelocity * deltaTime + acceleration * deltaTime * deltaTime;
	mVelocity += acceleration * deltaTime;

	float g = -9.8f * mMass; // Downwards force due to gravity;

	mForceAccumulator = { 0.0f, g, 0.0f };
}