#include "Particle3D.h"

Particle3D::Particle3D(const glm::vec3& inPosition, const glm::vec3& inVelocity, float inMass, float inRadius)
{
	mPosition = inPosition;
	mVelocity = inVelocity;
	mMass = inMass;
	mRadius = inRadius;

	mAcceleration = { 0.0f, -9.8f, 0.0f }; // Acceleration due to gravity

	float g = -9.8f * mMass; // Downwards force due to gravity;
	mForceAccumulator = { 0.0f, g, 0.0f };
}

void Particle3D::StepParticle(float deltaTime)
{

	mPosition += mVelocity * deltaTime + mAcceleration * deltaTime * deltaTime;
	mVelocity += mAcceleration * deltaTime;

	float g = -9.8f * mMass; // Downwards force due to gravity;

	mForceAccumulator = { 0.0f, g, 0.0f };
}