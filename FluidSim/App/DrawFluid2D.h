#pragma once

#include <memory>

#include "Fluid.h"
#include "Shader.h"

#include "glm/vec3.hpp"
#include "glm/matrix.hpp"

enum class DrawMode {
	Lines, Loop, Strip, Points
};

class DrawFluid2D
{
public:
	DrawFluid2D() = default;
	~DrawFluid2D() = default;

	DrawFluid2D(const Fluid2D& inFluid);

	void FromFluid(const Fluid2D& inFluid);

	void UpdateOpenGLBuffers();

	std::vector<glm::vec2> mParticlePoints;
};