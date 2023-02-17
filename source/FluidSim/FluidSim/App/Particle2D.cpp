#include "Particle2D.h"

Particle2D::Particle2D(const glm::dvec2& inPosition, const glm::dvec2& inVelocity, double inMass, double inRadius)
{
	mPosition = inPosition;
	mVelocity = inVelocity;
	mPrevVelocity = inVelocity;
	mMass = inMass;
	mRadius = inRadius;

	mAcceleration = { 0.0, -9.8f}; // Acceleration due to gravity

	double g = -9.8f * mMass; // Downwards force due to gravity;
	mForceAccumulator = { 0.0, g};
}

void Particle2D::StepRK3(double deltaTime, const glm::dvec2& K1, const glm::dvec2& K2, const glm::dvec2& K3)
{
	double scaledDt = deltaTime / 9.0;

	mPosition += ((2.0 * K1) + (3.0 * K2) + (4.0 * K3)) * scaledDt;
}