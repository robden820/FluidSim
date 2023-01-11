#pragma once

#include <memory>

#include "DrawFluid.h"

#include "Fluid2D.h"
#include "Shader.h"

#include "glm/vec2.hpp"
#include "glm/matrix.hpp"

class DrawFluid2D : public DrawFluid
{
public:
	DrawFluid2D() = default;
	~DrawFluid2D() = default;

	DrawFluid2D(const Fluid2D& inFluid);

	void FromFluid(const Fluid2D& inFluid);

	void UpdateOpenGLBuffers() override;

	std::vector<glm::vec2> mParticlePoints;
};