#pragma once

#include <tuple>
#include "glm/vec3.hpp"

class Particle2D
{
public:
	Particle2D() = default;
	~Particle2D() = default;

	Particle2D(const float inXPos, const float inYPos,  const float inXVel = 0.0f, const float inYVel = 0.0f, float inMass = 1.0f, float inRadius = 0.1f);

	void StepParticle(float deltaTime);

	const std::tuple<float, float> GetXPosition() const { return std::tuple<float, float>(mXPos, mYPos); }
	void SetPosition(float inXPos, float inYPos) { mXPos = inXPos; mYPos = inYPos; }

	const std::tuple<float, float> GetVelocity() const { return std::tuple<float, float>(mXVelocity, mYVelocity); }
	void SetVelocity(float inXVel, float inYVel) { mXVelocity = inXVel; mYVelocity = inYVel; }

	const std::tuple<float, float> GetAcceleration() const { return std::tuple<float, float>(mXAcceleration, mYAcceleration); }
	void SetAcceleration(float inXAcc, float inYAcc) { mXAcceleration = inXAcc; mYAcceleration = inYAcc; }

	const float GetMass() const { return mMass; }
	const float GetRadius() const { return mRadius; }

	void ApplyForce(float inXForce, float inYForce) { mXForce += inXForce; mYForce += inYForce; }

private:
	float mXPos;
	float mYPos;
	
	float mXVelocity;
	float mYVelocity;

	float mXAcceleration;
	float mYAcceleration;

	float mXForce;
	float mYForce;

	float mMass;
	float mRadius;
};