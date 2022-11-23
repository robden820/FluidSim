#pragma once

#include <memory>

#include "Fluid.h"
#include "Shader.h"

#include "glm/vec3.hpp"
#include "glm/matrix.hpp"

enum class DrawMode {
	Lines, Loop, Strip, Points
};

class DrawFluid
{
	public:
		DrawFluid() = default;
		~DrawFluid() = default;

		DrawFluid(const Fluid& inFluid);

		void FromFluid(const Fluid& inFluid);

		void UpdateOpenGLBuffers();
		
		std::vector<glm::vec3> mParticlePoints;
};