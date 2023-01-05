#pragma once

#include <vector>
#include <memory>

#include "glm/glm.hpp"

#include "Particle2D.h"
#include "Domain2D.h"
#include "MACGrid2D.h"

class Fluid2D
{
public:
	Fluid2D() = default;
	~Fluid2D() = default;

	Fluid2D(int numParticles);

	void StepSimulation(float deltaTime);

	const std::vector<Particle2D>& GetParticles() const { return mParticles; }
	const Particle2D& GetParticle(int index) const { return mParticles[index]; }
	int GetNumParticles() const { return mParticles.size(); }

	const Domain& GetDomain() const { return mDomain; }
	const MACGrid& GetMACGrid() const { return mMACGrid; }

	int GetMACGridResolution() const { return mMACGridResolution; }

	void ClampParticleToDomain(Particle2D& particle);

private:

	void InterpolateToGrid();
	void InterpolateFromGrid();

	int ClosestCellToParticle(const Particle2D& particle);

	std::vector<Particle2D> mParticles;
	std::vector<glm::vec3> mParticlePositions;

	MACGrid mMACGrid;
	int mMACGridResolution;

	Domain mDomain;
};

