#pragma once

#include "glm/vec3.hpp"

class GridNode
{
	public:
		GridNode() = default;
		~GridNode() = default;

		GridNode(glm::vec3 inPosition);

		const glm::vec3& GetPosition() { return mPosition; }

		const glm::vec3& GetVelocity() { return mVelocity; }
		void SetVelocity(const glm::vec3& inVelocity) { mVelocity = inVelocity; }

	private:
		glm::vec3 mPosition;
		glm::vec3 mVelocity;
};