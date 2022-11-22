#include "Particle.h"

Particle::Particle(glm::vec3 inPosition, glm::vec3 inVelocity, float inMass, float inRadius)
{
	mPosition = inPosition;
	mVelocity = inVelocity;
	mMass = inMass;
	mRadius = inRadius;

	mAcceleration = { 0.0f, -9.8f, 0.0f }; // Acceleration due to gravity
}

void Particle::StepParticle(float deltaTime)
{

	mPosition += mVelocity * deltaTime + mAcceleration * deltaTime * deltaTime;
	mVelocity += mAcceleration * deltaTime;

	float g = -9.8f * mMass; // Downwards force due to gravity;

	mForceAccumulator = { 0.0f, g, 0.0f };
}