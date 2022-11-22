#include "Fluid.h"
#include <iostream>

Fluid::Fluid(float numParticles, float gridResolution)
{
	mParticles.reserve(numParticles);

	for (int x = 0; x < 10; x++)
	{
		for (int y = 0; y < 10; y++)
		{
			for (int z = 0; z < 10; z++)
			{
				glm::vec3 position(rand() % 20 - 10, rand() % 20, rand() % 20 - 10);
				glm::vec3 velocity(rand() % 100 * 0.01f, rand() % 100 * 0.01f, rand() % 100 * 0.01f);

				Particle particle(position, velocity);

				mParticles.push_back(particle);
			}
		}
	}

	Domain d(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(10.0f, 10.0f, 10.0f));
	mDomain = d;

	ContactResolver resolver(100);
	mContactResolver = resolver;
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

	CheckParticleCollisions();
	if (mContacts.size() > 0)
	{
		mContactResolver.ResolveContacts(mContacts.data(), mContacts.size(), deltaTime);
		mContacts.clear();
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

void Fluid::CheckParticleCollisions()
{
	for (int p1 = 0; p1 < GetNumParticles(); p1++)
	{
		for (int p2 = p1 + 1; p2 < GetNumParticles(); p2++)
		{
			glm::vec3 dist = mParticles[p1].GetPosition() - mParticles[p2].GetPosition();
			
			if (glm::length(dist) < mParticles[p1].GetRadius() + mParticles[p2].GetRadius())
			{
				Contact contact(mParticles[p1], mParticles[p2], 1.0f, mParticles[p1].GetRadius() + mParticles[p2].GetRadius() - glm::length(dist) , glm::normalize(dist));

				mContacts.push_back(contact);
			}
		}
	}
}