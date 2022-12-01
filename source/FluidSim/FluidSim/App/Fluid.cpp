#include "Fluid.h"
#include <iostream>

Fluid::Fluid(int numParticles)
{
	// Initialize particles
	mParticles.reserve(numParticles);
	mParticlePositions.reserve(numParticles);

	for (int x = 0; x < 10; x++)
	{
		for (int y = 0; y < 10; y++)
		{
			for (int z = 0; z < 10; z++)
			{
				glm::vec3 position((x + 2.5) * 0.5f, (y + 5) * 0.5f, (z + 2.5) * 0.5f);

				Particle particle(position);

				mParticles.push_back(particle);
				mParticlePositions.push_back(position);
			}
		}
	}

	// Initialize simulation Domain.
	Domain d(glm::vec3(5.0f, 5.0f, 5.0f), glm::vec3(5.0f, 5.0f, 5.0f));
	mDomain = d;

	// Initialize Grid.
	MACGrid g(mDomain, mParticlePositions, 20.f);
	mMACGrid = g;
}

void Fluid::StepSimulation(float deltaTime)
{

	// Transfer particle velocities to grid.
	InterpolateToGrid();

	// Update grid velocities.
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
	for (int c = 0; c < mMACGrid.GetNumCells(); c++)
	{
		if (mMACGrid.GetCellType(c) != MACGrid::CellType::eSOLID)
		{
			mMACGrid.SetCellType(c, MACGrid::CellType::eAIR);
		}
	}

	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(mParticles[p]);

		mMACGrid.SetCellType(cellIndex, MACGrid::CellType::eFLUID);

		mMACGrid.SetCellVelocity(p, mParticles[p].GetVelocity());
	}
}

void Fluid::InterpolateFromGrid()
{
	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(mParticles[p]);

		int x, y, z;
		std::tie(x, y, z) = mMACGrid.GetXYZFromIndex(cellIndex);

		int neighbourLeft = mMACGrid.GetIndexFromXYZ(x - 1, y, z);
		int neighbourBottom = mMACGrid.GetIndexFromXYZ(x, y - 1, z);
		int neighbourBack = mMACGrid.GetIndexFromXYZ(x, y, z - 1);

		float velocityX = (mMACGrid.GetCellXVelocity(cellIndex) + mMACGrid.GetCellXVelocity(neighbourLeft)) * 0.5f;
		float velocityY = (mMACGrid.GetCellYVelocity(cellIndex) + mMACGrid.GetCellYVelocity(neighbourBottom)) * 0.5f;
		float velocityZ = (mMACGrid.GetCellZVelocity(cellIndex) + mMACGrid.GetCellZVelocity(neighbourBack)) * 0.5f;

		mParticles[p].SetVelocity(glm::vec3(velocityX, velocityY, velocityZ));
	}
}

int Fluid::ClosestCellToParticle(const Particle& particle)
{
	glm::vec3 particlePos = particle.GetPosition();

	float shortestDist = 100.0f;

	int closestCell = 0;

	for (int n = 0; n < mMACGrid.GetNumCells(); n++)
	{
		glm::vec3 pToG = mMACGrid.GetCellCenter(n) - particlePos;

		float distSqr = (pToG.x * pToG.x) + (pToG.y * pToG.y) + (pToG.z * pToG.z);

		if (distSqr < shortestDist)
		{
			shortestDist = distSqr;
			closestCell = n;
		}
	}

	return closestCell;
}

Particle& Fluid::ClosestParticleToCell(const glm::vec3& cellCenter)
{
	float shortestDist = 100.0f;

	Particle& closestParticle = mParticles[0];

	for (int p = 0; p < GetNumParticles(); p++)
	{
		glm::vec3 GToP = mParticles[p].GetPosition() - cellCenter;

		float distSqr = (GToP.x * GToP.x) + (GToP.y * GToP.y) + (GToP.z * GToP.z);

		if (distSqr < shortestDist)
		{
			shortestDist = distSqr;
			closestParticle = mParticles[p];
		}
	}

	return closestParticle;
}