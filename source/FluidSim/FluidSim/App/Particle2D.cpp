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

void Particle2D::StepRK3(float deltaTime, const glm::vec2& K1, const glm::vec2& K2, const glm::vec2& K3)
{
	float scaledDt = deltaTime / 9.0f;

	mPosition += (2.0f * scaledDt * K1) + (3.0f * scaledDt * K2) + (4.0f * scaledDt * K3);
}