#pragma once

#include "glm/vec3.hpp"

class Particle
{
	public:
		Particle() = default;
		~Particle() = default;

		Particle(glm::vec3 inPosition, glm::vec3 inVelocity = glm::vec3(0.0f, 0.0f, 0.0f));

		glm::vec3 GetPosition() { return mPosition; }
		void SetPosition(glm::vec3 inPos) { mPosition = inPos; }

		glm::vec3 GetVelocity() { return mVelocity; }
		void SetVelocity(glm::vec3 inVel) { mVelocity = inVel; }

	private:
		glm::vec3 mPosition;
		glm::vec3 mVelocity;
};