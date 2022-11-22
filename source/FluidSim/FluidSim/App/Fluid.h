#pragma once

#include <vector>
#include <memory>

#include "Particle.h"
#include "GridNode.h"
#include "Domain.h"

class Fluid
{
	public:
		Fluid() = default;
		~Fluid() = default;

		Fluid(float numParticles, float gridResolution);
		Fluid(std::vector<Particle> inParticles, std::vector<GridNode> inGridNodes);

		void StepSimulation(float deltaTime);

		void SetParticles(std::vector<Particle> inParticles) { mParticles = inParticles; }
		void SetGridNodes(std::vector<GridNode> inGridNodes) { mGridNodes = inGridNodes; }

		int GetNumParticles() { return mParticles.size(); }

		void ClampParticleToDomain(Particle& particle);
	
		std::vector<Particle> mParticles;
		std::vector<GridNode> mGridNodes;

		Domain mDomain;
};