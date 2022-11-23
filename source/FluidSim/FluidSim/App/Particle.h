#pragma once

#include "glm/vec3.hpp"

class Particle
{
	public:
		Particle() = default;
		~Particle() = default;

		Particle(const glm::vec3& inPosition, const glm::vec3& inVelocity = glm::vec3(0.0f, 0.0f, 0.0f), float inMass = 1.0f, float inRadius = 0.1f);

		void StepParticle(float deltaTime);

		const glm::vec3& GetPosition() const { return mPosition; }
		void SetPosition(glm::vec3 inPos) { mPosition = inPos; }

		const glm::vec3& GetVelocity() const { return mVelocity; }
		void SetVelocity(const glm::vec3& inVel) { mVelocity = inVel; }

		const glm::vec3& GetAcceleration() const { return mAcceleration; }
		void SetAcceleration(const glm::vec3& inAcc) { mAcceleration = inAcc; }

		const float GetMass() const { return mMass; }
		const float GetRadius() const{ return mRadius; }

		void ApplyForce(const glm::vec3& inForce) { mForceAccumulator += inForce; }

	private:
		glm::vec3 mPosition;
		glm::vec3 mVelocity;
		glm::vec3 mAcceleration;

		glm::vec3 mForceAccumulator;

		float mMass;
		float mRadius;
};