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

		std::vector<Particle>& GetParticles() { return mParticles; }
		void SetParticles(std::vector<Particle> inParticles) { mParticles = inParticles; }
		int GetNumParticles() { return mParticles.size(); }

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