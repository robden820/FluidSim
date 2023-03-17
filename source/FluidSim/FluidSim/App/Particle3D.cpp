#include "Particle3D.h"

Particle3D::Particle3D(const glm::dvec3& inPosition, const glm::dvec3& inVelocity, double inMass)
{
	mPosition = inPosition;
	mVelocity = inVelocity;
	mMass = inMass;

	mAcceleration = { 0.0, -9.8f, 0.0 }; // Acceleration due to gravity

	double g = -9.8f * mMass; // Downwards force due to gravity;
	mForceAccumulator = { 0.0, g, 0.0 };
}

void Particle3D::StepParticle(double deltaTime)
{

	mPosition += mVelocity * deltaTime + mAcceleration * deltaTime * deltaTime;
	mVelocity += mAcceleration * deltaTime;

	double g = -9.8f * mMass; // Downwards force due to gravity;

	mForceAccumulator = { 0.0, g, 0.0 };
}