#pragma once

#include <vector>
#include <memory>

#include "glm/glm.hpp"

#include "Particle.h"
#include "GridNode.h"
#include "Domain.h"
#include "Grid.h"

class Fluid
{
	public:
		Fluid() = default;
		~Fluid() = default;

		Fluid(float numParticles);

		void StepSimulation(float deltaTime);

		const std::vector<Particle>& GetParticles() const { return mParticles; }
		int GetNumParticles() const { return mParticles.size(); }

		void ClampParticleToDomain(Particle& particle);
	
	private:

		void InterpolateToGrid();
		void InterpolateFromGrid();
		
		GridNode& ClosestNodeToParticle(Particle& particle);
		Particle& ClosestParticleToNode(GridNode& node);

		std::vector<Particle> mParticles;
		Grid mGrid;
		
		Domain mDomain;
};