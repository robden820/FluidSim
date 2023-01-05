#include "Particle2D.h"

Particle2D::Particle2D(const float inXPos, const float inYPos, const float inXVel, const float inYVel, float inMass, float inRadius)
{
	mXPos = inXPos;
	mYPos = inYPos;

	mXVelocity = inXVel;
	mYVelocity = inYVel;

	mMass = inMass;
	mRadius = inRadius;

	mXAcceleration = 0.0f;
	mYAcceleration = -9.8f;

	float g = -9.8f * mMass; // Downwards force due to gravity;
	mYForce += g;
}

void Particle2D::StepParticle(float deltaTime)
{
	mXPos += mXVelocity * deltaTime + mXAcceleration * deltaTime * deltaTime;
	mYPos += mYVelocity * deltaTime + mYAcceleration * deltaTime * deltaTime;

	mXVelocity += mXAcceleration * deltaTime;
	mYVelocity += mYAcceleration * deltaTime;
}