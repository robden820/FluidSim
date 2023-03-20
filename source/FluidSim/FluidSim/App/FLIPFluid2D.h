#pragma once

#include <vector>
#include <memory>

#include "glm/glm.hpp"

#include "Fluid2D.h"

#include "Particle2D.h"
#include "MACGrid2D.h"

class FLIPFluid2D : public Fluid2D<Particle2D>
{
	enum class SimulationType : bool
	{
		ePIC = false,
		eFLIP = true,
	};

public:
	FLIPFluid2D() = default;
	~FLIPFluid2D() = default;

	FLIPFluid2D(const ApplicationData& inOutData);

	void StepParticles(double deltaTime, const MACGrid2D& inMACGrid) override; // RK3
	void StepParticlesEuler(double deltaTime, const MACGrid2D& inMACGrid) override;

	void UpdateApplicationData(ApplicationData& inOutData) override;

	void InterpolateToGrid(MACGrid2D& inMACGrid) override;
	
	void InterpolateFromGrid(const MACGrid2D& inMACGrid) override;
	void InterpolateFromGridBSpline(const MACGrid2D& inMACGrid) override;

private:

	// Interpolate from current cell velocities.
	glm::dvec2 InterpolateFromGridCell(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex, SimulationType simType);
	glm::dvec2 InterpolateFromGridCell(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex, SimulationType simType);

	glm::dvec2 InterpolateFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex, SimulationType simType);
	glm::dvec2 InterpolateFromGridCellBSpline(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex, SimulationType simType);

	// Interpolate from previous cell velocities
	glm::dvec2 InterpolateFromGridCellBSplinePrev(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex);
	glm::dvec2 InterpolateFromGridCellBSplinePrev(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex);

	double mFLIPBlend; // % FLIP Blend, clamped between [0, 1]
};

