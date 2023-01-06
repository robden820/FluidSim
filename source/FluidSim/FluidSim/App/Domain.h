#pragma once

#include "glm/vec3.hpp"

class Domain
{
	public:
		virtual ~Domain() = default;

	protected:
		virtual void CalculateSides() = 0;
};