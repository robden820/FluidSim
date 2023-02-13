#pragma once

#include <vector>
#include <memory>

#include "glm/glm.hpp"

#include "Fluid.h"

#include "Particle2D.h"
#include "MACGrid2D.h"

class Fluid2D : public Fluid
{
public:
	Fluid2D() = default;
	~Fluid2D() = default;

	Fluid2D(const ApplicationData& inOutData);

	void StepParticles(float deltaTime, const MACGrid2D& inMACGrid); // RK3
	void StepParticlesFLIP(float deltaTime, const MACGrid2D& inMACGrid);
	void UpdateApplicationData(ApplicationData& inOutData);

	void InterpolateToGrid(MACGrid2D& inMACGrid);

	void InterpolateFromGrid(const MACGrid2D& inMACGrid);
	void InterpolateFromGridBSpline(const MACGrid2D& inMACGrid);

	void InterpolateFromGridFLIP(const MACGrid2D& inMACGrid);
	void InterpolateFromGridBSplineFLIP(const MACGrid2D& inMACGrid);

	const std::vector<Particle2D>& GetParticles() const { return mParticles; }
	const Particle2D& GetParticle(int index) const { return mParticles[index]; }
	size_t GetNumParticles() const { return mParticles.size(); }

	void DeleteBadParticles(const MACGrid2D& inMACGrid);

	// TO DO: remove the use of these functions. Requires updating Fluid class. Will do this when updating Fluid3D.
	void Update(ApplicationData& inOutData) override {};
	void InterpolateToGrid() override {};
	void InterpolateFromGrid() override {};

private:

	void SeedParticles(const ApplicationData& inOutData);

	float InterpolateSupport(const glm::vec2& diff, float invCellSize);
	float BSpline(float input);

	int ClosestCellToParticle(const MACGrid2D& inMACGrid, const Particle2D& particle);

	glm::dvec2 InterpolateFromGridCell(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex);
	glm::dvec2 InterpolateFromGridCell(const MACGrid2D& inMACGrid, const glm::vec2& particlePosition, int cellIndex);

	glm::dvec2 InterpolateFromGridCellFLIP(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex);
	glm::dvec2 InterpolateFromGridCellFLIP(const MACGrid2D& inMACGrid, const glm::vec2& particlePosition, int cellIndex);

	glm::dvec2 InterpolateFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex);
	glm::dvec2 InterpolateFromGridCellBSpline(const MACGrid2D& inMACGrid, const glm::vec2& particlePosition, int cellIndex);

	void ProjectParticleToFluid(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex);

	std::vector<Particle2D> mParticles;
	std::vector<glm::vec2> mParticlePositions;
};

