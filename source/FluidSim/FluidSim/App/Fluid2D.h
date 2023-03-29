#pragma once

#include "Fluid.h"

#include "glm/glm.hpp"

#include "MACGrid2D.h"

#include "Particle2D.h"
#include "APICParticle2D.h"
#include "RPICParticle2D.h"

class IFluid2D : public Fluid
{
public:
	virtual ~IFluid2D() = default;

	virtual void StepParticles(double deltaTime, const MACGrid2D& inMACGrid) = 0;
	virtual void StepParticlesEuler(double deltaTime, const MACGrid2D& inMACGrid) = 0;

	virtual void UpdateApplicationData(ApplicationData& inOutData) = 0;

	virtual void InterpolateToGrid(MACGrid2D& inMACGrid) = 0;

	// FLIP/PIC blend.
	virtual void InterpolateFromGrid(const MACGrid2D& inMACGrid) = 0;
	virtual void InterpolateFromGridBSpline(const MACGrid2D& inMACGrid) = 0;

	virtual void DeleteBadParticles(const MACGrid2D& inMACGrid) = 0;

	virtual void CalculateSystemEnergy() = 0;
};

template <typename Tparticle>
class Fluid2D : public IFluid2D
{
public:
	virtual ~Fluid2D() = default;

	void SeedParticles(const ApplicationData& inOutData);

	int ClosestCellToParticle(const MACGrid2D& inMACGrid, const Tparticle& particle);
	void ProjectParticleToFluid(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex);

	void DeleteBadParticles(const MACGrid2D& inMACGrid) override;

	const std::vector<Tparticle>& GetParticles() const { return mParticles; }
	const Tparticle& GetParticle(int index) const { return mParticles[index]; }
	size_t GetNumParticles() const { return mParticles.size(); }

	void CalculateSystemEnergy() override;

protected:
	double InterpolateSupport(const glm::dvec2& diff, double invCellSize);
	double BSpline(double input);

	std::vector<Tparticle> mParticles;
	std::vector<glm::dvec2> mParticlePositions;
};

template class Fluid2D<Particle2D>;
template class Fluid2D<APICParticle2D>;
template class Fluid2D<RPICParticle2D>;