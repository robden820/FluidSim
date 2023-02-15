#pragma once

#include <vector>
#include <memory>

#include "glm/glm.hpp"

#include "Particle.h"
#include "MACGrid.h"

#include "ApplicationData.h"

class Fluid
{
	public:

		enum class SimulationType : bool
		{
			ePIC = false,
			eFLIP = true,
		};

		virtual ~Fluid() = default;

		virtual void Update(ApplicationData& inOutData) = 0;
	
	protected:

		virtual void InterpolateToGrid() = 0;
		virtual void InterpolateFromGrid() = 0;
};