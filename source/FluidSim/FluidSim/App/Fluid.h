#pragma once

#include <vector>
#include <memory>

#include "glm/glm.hpp"

#include "Particle.h"
#include "Domain.h"
#include "MACGrid.h"

class Fluid
{
	public:
		virtual ~Fluid() = default;

		virtual void StepSimulation(float deltaTime) = 0;
	
	protected:

		virtual void InterpolateToGrid() = 0;
		virtual void InterpolateFromGrid() = 0;

		int mMACGridResolution;
};