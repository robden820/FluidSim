#pragma once

#include <memory>

#include "DrawFluid.h"

#include "Shader.h"

enum class DrawMode {
	Lines, Loop, Strip, Points
};


class DrawFluid
{
public:
	virtual ~DrawFluid() = default;

	virtual void UpdateOpenGLBuffers() = 0;
};
