#include "Fluid.h"

Fluid::Fluid(float numParticles, float gridResolution)
{
	mParticles.reserve(numParticles);

	for (int x = -5; x < 5; x++)
	{
		for (int y = 0; y < 10; y++)
		{
			for (int z = 0; z < 10; z++)
			{
				glm::vec3 position(x, y, z);
				glm::vec3 velocity(x, y, z);

				Particle particle(position, velocity);

				mParticles.push_back(particle);
			}
		}
	}

	Domain d(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(10.0f, 10.0f, 10.0f));
	mDomain = d;
}

Fluid::Fluid(std::vector<Particle> inParticles, std::vector<GridNode> inGridNodes)
{
	mParticles = inParticles;
	mGridNodes = inGridNodes;
}

void Fluid::StepSimulation(float deltaTime)
{
	for (int p = 0; p < GetNumParticles(); p++)
	{
		mParticles[p].StepParticle(deltaTime);

		if (!mDomain.IsPointInDomain(mParticles[p].GetPosition()))
		{
			ClampParticleToDomain(mParticles[p]);
		}
	}
}

void Fluid::ClampParticleToDomain(Particle& particle)
{
	glm::vec3 particlePos = particle.GetPosition();
	glm::vec3 particleVel = particle.GetVelocity();

	if (particlePos.x < mDomain.GetLeft())
	{
		particlePos.x = mDomain.GetLeft();
		particleVel.x *= -1.0f;
	}
	else if (particlePos.x > mDomain.GetRight())
	{
		particlePos.x = mDomain.GetRight();
		particleVel.x *= -1.0f;
	}

	if (particlePos.y < mDomain.GetBottom())
	{
		particlePos.y = mDomain.GetBottom();
		particleVel.y *= -1.0f;
	}
	else if (particlePos.y > mDomain.GetTop())
	{
		particlePos.y = mDomain.GetTop();
		particleVel.y *= -1.0f;
	}

	if (particlePos.z < mDomain.GetBack())
	{
		particlePos.z = mDomain.GetBack();
		particleVel.z *= -1.0f;
	}
	else if (particlePos.z > mDomain.GetFront())
	{
		particlePos.z = mDomain.GetFront();
		particleVel.z *= -1.0f;
	}

	particle.SetPosition(particlePos);
	particle.SetVelocity(particleVel);
}