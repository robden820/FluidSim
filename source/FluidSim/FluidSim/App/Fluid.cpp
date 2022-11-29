#include "Fluid.h"
#include <iostream>

Fluid::Fluid(float numParticles)
{
	// Initialize particles
	mParticles.reserve(numParticles);

	for (int x = 0; x < 10; x++)
	{
		for (int y = 0; y < 10; y++)
		{
			for (int z = 0; z < 10; z++)
			{
				glm::vec3 position((x + 2.5) * 0.5f, (y + 5) * 0.5f, (z + 2.5) * 0.5f);

				Particle particle(position);

				mParticles.push_back(particle);
			}
		}
	}
	
	// Initialize simulation Domain.
	Domain d(glm::vec3(5.0f, 5.0f, 5.0f), glm::vec3(5.0f, 5.0f, 5.0f));
	mDomain = d;

	// Initialize Grid.
	//Grid g(mDomain, 0.5f);
	//mGrid = g;

	MACGrid g(mDomain, 0.5f);
	mMACGrid = g;


//	ContactResolver resolver(100);
//	mContactResolver = resolver;
}

void Fluid::StepSimulation(float deltaTime)
{

	// Transfer particle velocities to grid.
	InterpolateToGrid();

	// Update grid velocities.
	//mGrid.StepGrid(deltaTime);
	mMACGrid.Update(deltaTime);

	// Interpolate velocities back to particles.
	InterpolateFromGrid();

	// Advect particles.
	for (int p = 0; p < GetNumParticles(); p++)
	{
		
		mParticles[p].StepParticle(deltaTime);
	}
	
	
	// Ensure particles stay inside the simulation domain.
	for (int p = 0; p < GetNumParticles(); p++)
	{
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
		particleVel.x = 0.0f;
	}
	else if (particlePos.x > mDomain.GetRight())
	{
		particlePos.x = mDomain.GetRight();
		particleVel.x = 0.0f;
	}

	if (particlePos.y < mDomain.GetBottom())
	{
		particlePos.y = mDomain.GetBottom();
		particleVel.y = 0.0f;
	}
	else if (particlePos.y > mDomain.GetTop())
	{
		particlePos.y = mDomain.GetTop();
		particleVel.y = 0.0f;
	}

	if (particlePos.z < mDomain.GetBack())
	{
		particlePos.z = mDomain.GetBack();
		particleVel.z = 0.0f;
	}
	else if (particlePos.z > mDomain.GetFront())
	{
		particlePos.z = mDomain.GetFront();
		particleVel.z = 0.0f;
	}

	particle.SetPosition(particlePos);
	particle.SetVelocity(particleVel);
}

void Fluid::InterpolateToGrid()
{
	/*
	for (int n = 0; n < mGrid.GetNumGridNodes(); n++)
	{
		Particle& particle = ClosestParticleToNode(mGrid.GetGridNode(n));

		mGrid.GetGridNode(n).SetVelocity(particle.GetVelocity());
	}
	*/

	for (int n = 0; n < mMACGrid.GetNumCells(); n++)
	{
		Particle& particle = ClosestParticleToCell(mMACGrid.GetCellCenter(n));

		mMACGrid.SetGridCellVelocity(n, particle.GetVelocity());
	}
}

void Fluid::InterpolateFromGrid()
{
	/*
	for (int p = 0; p < GetNumParticles(); p++)
	{
		GridNode& node = ClosestNodeToParticle(mParticles[p]);

		mParticles[p].SetVelocity(node.GetVelocity());
	}
	*/

	for (int p = 0; p < GetNumParticles(); p++)
	{
		const MACGridCell& cell = ClosestCellToParticle(mParticles[p]);

		glm::vec3 pVelocity(0.f, 0.f, 0.f);
		cell.GetCellVelocity(pVelocity);

		mParticles[p].SetVelocity(pVelocity);
	}
}

GridNode& Fluid::ClosestNodeToParticle(Particle& particle)
{
	glm::vec3 particlePos = particle.GetPosition();

	float shortestDist = 100.0f;

	GridNode& closestNode = mGrid.GetGridNode(0);

	for (int n = 0; n < mGrid.GetNumGridNodes(); n++)
	{
		glm::vec3 pToG = mGrid.GetGridNode(n).GetPosition() - particlePos;

		float distSqr = (pToG.x * pToG.x) + (pToG.y * pToG.y) + (pToG.z * pToG.z);

		if (shortestDist < distSqr)
		{
			shortestDist = distSqr;
			closestNode = mGrid.GetGridNode(n);
		}
	}

	return closestNode;
}

const MACGridCell& Fluid::ClosestCellToParticle(const Particle& particle) const
{
	glm::vec3 particlePos = particle.GetPosition();

	float shortestDist = 100.0f;

	int closestCell = 0;

	for (int n = 0; n < mMACGrid.GetNumCells(); n++)
	{
		glm::vec3 pToG = mMACGrid.GetCellCenter(n) - particlePos;

		float distSqr = (pToG.x * pToG.x) + (pToG.y * pToG.y) + (pToG.z * pToG.z);

		if (shortestDist < distSqr)
		{
			shortestDist = distSqr;
			closestCell = n;
		}
	}

	return mMACGrid.GetGridCell(closestCell);
}

Particle& Fluid::ClosestParticleToNode(GridNode& node)
{
	glm::vec3 nodePos = node.GetPosition();

	float shortestDist = 100.0f;

	Particle& closestParticle = mParticles[0];

	for (int p = 0; p < GetNumParticles(); p++)
	{
		glm::vec3 GToP = mParticles[p].GetPosition() - nodePos;

		float distSqr = (GToP.x * GToP.x) + (GToP.y * GToP.y) + (GToP.z * GToP.z);

		if (shortestDist < distSqr)
		{
			shortestDist = distSqr;
			closestParticle = mParticles[p];
		}
	}

	return closestParticle;
}

Particle& Fluid::ClosestParticleToCell(const glm::vec3& cellCenter)
{
	float shortestDist = 100.0f;

	Particle& closestParticle = mParticles[0];

	for (int p = 0; p < GetNumParticles(); p++)
	{
		glm::vec3 GToP = mParticles[p].GetPosition() - cellCenter;

		float distSqr = (GToP.x * GToP.x) + (GToP.y * GToP.y) + (GToP.z * GToP.z);

		if (shortestDist < distSqr)
		{
			shortestDist = distSqr;
			closestParticle = mParticles[p];
		}
	}

	return closestParticle;
}