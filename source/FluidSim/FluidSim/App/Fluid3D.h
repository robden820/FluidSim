#pragma once

#include <vector>
#include <memory>

#include "glm/glm.hpp"

#include "Fluid.h"

#include "Particle3D.h"
#include "Domain3D.h"
#include "MACGrid3D.h"

class Fluid3D : public Fluid
{
public:
	Fluid3D() = default;
	~Fluid3D() = default;

	Fluid3D(int numParticles);

	void StepSimulation(float deltaTime) override;

	const std::vector<Particle3D>& GetParticles() const { return mParticles; }
	const Particle3D& GetParticle(int index) const { return mParticles[index]; }
	int GetNumParticles() const { return mParticles.size(); }

	const Domain3D& GetDomain() const { return mDomain; }
	const MACGrid3D& GetMACGrid() const { return mMACGrid; }

	int GetMACGridResolution() const { return mMACGridResolution; }

	void ClampParticleToDomain(Particle3D& particle);

private:

	void InterpolateToGrid() override;
	void InterpolateFromGrid() override;

	int ClosestCellToParticle(const Particle3D& particle);

	std::vector<Particle3D> mParticles;
	std::vector<glm::vec3> mParticlePositions;

	MACGrid3D mMACGrid;
	Domain3D mDomain;
};
