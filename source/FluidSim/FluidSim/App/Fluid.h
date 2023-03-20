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

		virtual ~Fluid() = default;
};