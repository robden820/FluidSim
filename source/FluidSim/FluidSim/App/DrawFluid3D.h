#pragma once

#include <memory>

#include "Fluid3D.h"
#include "Shader.h"

#include "glm/vec3.hpp"
#include "glm/matrix.hpp"

enum class DrawMode {
	Lines, Loop, Strip, Points
};

class DrawFluid3D
{
	public:
		DrawFluid3D() = default;
		~DrawFluid3D() = default;

		DrawFluid3D(const Fluid3D& inFluid);

		void FromFluid(const Fluid3D& inFluid);

		void UpdateOpenGLBuffers();
		
		std::vector<glm::vec3> mParticlePoints;
};