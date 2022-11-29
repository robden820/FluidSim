#pragma once

#include <vector>
#include <memory>

#include "glm/glm.hpp"

#include "Particle.h"
#include "GridNode.h"
#include "Domain.h"
#include "MACGrid.h"
#include "Grid.h"

class Fluid
{
	public:
		Fluid() = default;
		~Fluid() = default;

		Fluid(float numParticles);

		void StepSimulation(float deltaTime);

		const std::vector<Particle>& GetParticles() const { return mParticles; }
		const Particle& GetParticle(int index) const { return mParticles[index]; }
		int GetNumParticles() const { return mParticles.size(); }

		const Domain& GetDomain() const { return mDomain; }
//		const Grid& GetGrid() const { return mGrid; }
		const MACGrid& GetMACGrid() const { return mMACGrid; }

		void ClampParticleToDomain(Particle& particle);
	
	private:

		void InterpolateToGrid();
		void InterpolateFromGrid();
		
		GridNode& ClosestNodeToParticle(Particle& particle);
		const MACGridCell& ClosestCellToParticle(const Particle& particle) const;

		Particle& ClosestParticleToNode(GridNode& node);
		Particle& ClosestParticleToCell(const glm::vec3& cellCenter);

		std::vector<Particle> mParticles;
		Grid mGrid;

		MACGrid mMACGrid;
		
		Domain mDomain;
};