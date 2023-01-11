#pragma once

#include <memory>

#include "DrawFluid.h"

#include "Fluid3D.h"
#include "Shader.h"

#include "glm/vec3.hpp"
#include "glm/matrix.hpp"

class DrawFluid3D : public DrawFluid
{
	public:
		DrawFluid3D() = default;
		~DrawFluid3D() = default;

		DrawFluid3D(const Fluid3D& inFluid);

		void FromFluid(const Fluid3D& inFluid);

		void UpdateOpenGLBuffers() override;
		
		std::vector<glm::vec3> mParticlePoints;
};