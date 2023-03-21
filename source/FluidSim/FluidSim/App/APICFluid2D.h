#pragma once

#pragma once

#include <vector>
#include <memory>

#include "glm/glm.hpp"

#include "Fluid2D.h"

#include "APICParticle2D.h"
#include "MACGrid2D.h"

class APICFluid2D : public Fluid2D<APICParticle2D>
{
public:
	APICFluid2D() = default;
	~APICFluid2D() = default;

	APICFluid2D(const ApplicationData& inOutData);

	void StepParticles(double deltaTime, const MACGrid2D& inMACGrid) override; // RK3
	void StepParticlesEuler(double deltaTime, const MACGrid2D& inMACGrid) override;

	void UpdateApplicationData(ApplicationData& inOutData) override;

	void InterpolateToGrid(MACGrid2D& inMACGrid) override;

	void InterpolateFromGrid(const MACGrid2D& inMACGrid) override;
	void InterpolateFromGridBSpline(const MACGrid2D& inMACGrid) override;

private:

	// Interpolate from current cell velocities.
	glm::dvec2 InterpolateVelocityFromGridCell(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex);
	glm::dvec2 InterpolateVelocityFromGridCell(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex);

	glm::dvec2 InterpolateVelocityFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex);
	glm::dvec2 InterpolateVelocityFromGridCellBSpline(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex);

	void InterpolateAffineFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex);
};