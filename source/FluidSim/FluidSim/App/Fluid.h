#pragma once

#include <vector>
#include <memory>

#include "glm/glm.hpp"

#include "Particle.h"
#include "Domain.h"
#include "MACGrid.h"

class Fluid
{
	public:
		Fluid() = default;
		~Fluid() = default;

		Fluid(int numParticles);

		void StepSimulation(float deltaTime);

		const std::vector<Particle>& GetParticles() const { return mParticles; }
		const Particle& GetParticle(int index) const { return mParticles[index]; }
		int GetNumParticles() const { return mParticles.size(); }

		const Domain& GetDomain() const { return mDomain; }
		const MACGrid& GetMACGrid() const { return mMACGrid; }

		void ClampParticleToDomain(Particle& particle);
	
	private:

		void InterpolateToGrid();
		void InterpolateFromGrid();
		
		const MACGridCell& ClosestCellToParticle(const Particle& particle) const;

		Particle& ClosestParticleToCell(const glm::vec3& cellCenter);

		std::vector<Particle> mParticles;

		MACGrid mMACGrid;
		
		Domain mDomain;
};