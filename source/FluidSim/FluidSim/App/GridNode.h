#pragma once

#include "glm/vec3.hpp"

class GridNode
{
	public:
		GridNode() = default;
		~GridNode() = default;

	private:
		glm::vec3 mPosition;
		glm::vec3 mVelocity;
};