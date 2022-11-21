#include "Fluid.h"

Fluid::Fluid(float numParticles, float gridResolution)
{
	mParticles.reserve(numParticles);

	for(int p = 0; p < numParticles; p++)
	{
		glm::vec3 position(0.0f, 0.0f, 0.0f);

		Particle particle(position);

		mParticles.push_back(particle);
	}
}

Fluid::Fluid(std::vector<Particle> inParticles, std::vector<GridNode> inGridNodes)
{
	mParticles = inParticles;
	mGridNodes = inGridNodes;
}