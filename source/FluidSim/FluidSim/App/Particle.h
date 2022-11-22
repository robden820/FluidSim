#pragma once

#include "glm/vec3.hpp"

class Particle
{
	public:
		Particle() = default;
		~Particle() = default;

		Particle(glm::vec3 inPosition, glm::vec3 inVelocity = glm::vec3(0.0f, 0.0f, 0.0f), float inMass = 1.0f, float inRadius = 0.1f);

		void StepParticle(float deltaTime);

		glm::vec3 GetPosition() { return mPosition; }
		void SetPosition(glm::vec3 inPos) { mPosition = inPos; }

		glm::vec3 GetVelocity() { return mVelocity; }
		void SetVelocity(glm::vec3 inVel) { mVelocity = inVel; }

		glm::vec3 GetAcceleration() { return mAcceleration; }
		void SetAcceleration(glm::vec3 inAcc) { mAcceleration = inAcc; }

		float GetMass() { return mMass; }
		float GetRadius() { return mRadius; }

		void ApplyForce(glm::vec3 inForce) { mForceAccumulator += inForce; }

	private:
		glm::vec3 mPosition;
		glm::vec3 mVelocity;
		glm::vec3 mAcceleration;

		glm::vec3 mForceAccumulator;

		float mMass;
		float mRadius;
};