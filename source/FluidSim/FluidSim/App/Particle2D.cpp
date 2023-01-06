#include "Particle2D.h"

Particle2D::Particle2D(const glm::vec2& inPosition, const glm::vec2& inVelocity, float inMass, float inRadius)
{
	mPosition = inPosition;
	mVelocity = inVelocity;
	mMass = inMass;
	mRadius = inRadius;

	mAcceleration = { 0.0f, -9.8f}; // Acceleration due to gravity

	float g = -9.8f * mMass; // Downwards force due to gravity;
	mForceAccumulator = { 0.0f, g};
}

void Particle2D::StepParticle(float deltaTime)
{
	mPosition += mVelocity * deltaTime + mAcceleration * deltaTime * deltaTime;
	mVelocity += mAcceleration * deltaTime;

	float g = -9.8f * mMass; // Downwards force due to gravity;

	mForceAccumulator = { 0.0f, g};
}