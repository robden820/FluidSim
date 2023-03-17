#include "APICParticle2D.h"

APICParticle2D::APICParticle2D(const glm::dvec2& inPosition, const glm::dvec2& inVelocity, double inAngularMomentum, double inMass) 
{
	mPosition = inPosition;
	mVelocity = inVelocity;
	mAngularMomentum = inAngularMomentum;
	mMass = inMass;
	mAcceleration = { 0.0, -9.8f }; // Acceleration due to gravity

	double g = -9.8f * mMass; // Downwards force due to gravity;
	mForceAccumulator = { 0.0, g };
}