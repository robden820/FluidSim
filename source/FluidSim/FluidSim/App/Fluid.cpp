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

	std::vector<int> numParticlesToCell;
	std::vector<float> interpXVelocities;
	std::vector<float> interpYVelocities;
	std::vector<float> interpZVelocities;
	numParticlesToCell.assign(mMACGrid.GetNumCells(), 0);
	interpXVelocities.assign(mMACGrid.GetNumCells(), 0.f);
	interpYVelocities.assign(mMACGrid.GetNumCells(), 0.f);
	interpZVelocities.assign(mMACGrid.GetNumCells(), 0.f);

	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(mParticles[p]);

		if (mMACGrid.GetCellType(cellIndex) != MACGrid::CellType::eSOLID)
		{
			mMACGrid.SetCellType(cellIndex, MACGrid::CellType::eFLUID);

			interpXVelocities[cellIndex] += mParticles[p].GetVelocity().x;
			interpYVelocities[cellIndex] += mParticles[p].GetVelocity().y;
			interpZVelocities[cellIndex] += mParticles[p].GetVelocity().z;
			++numParticlesToCell[cellIndex];

			int x, y, z;
			std::tie(x, y, z) = mMACGrid.GetXYZFromIndex(cellIndex);

			if (x < 0)
			{
				int neighbourLeft = mMACGrid.GetIndexFromXYZ(x - 1, y, z);

				if (mMACGrid.GetCellType(cellIndex) != MACGrid::CellType::eSOLID)
				{
					interpXVelocities[neighbourLeft] += mParticles[p].GetVelocity().x;
					interpYVelocities[neighbourLeft] += mParticles[p].GetVelocity().y;
					interpZVelocities[neighbourLeft] += mParticles[p].GetVelocity().z;
					++numParticlesToCell[neighbourLeft];
				}
			}
			if (y < 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXYZ(x, y - 1, z);

				if (mMACGrid.GetCellType(cellIndex) != MACGrid::CellType::eSOLID)
				{
					interpXVelocities[neighbourBottom] += mParticles[p].GetVelocity().x;
					interpYVelocities[neighbourBottom] += mParticles[p].GetVelocity().y;
					interpZVelocities[neighbourBottom] += mParticles[p].GetVelocity().z;
					++numParticlesToCell[neighbourBottom];
				}
			}
			if (z < 0)
			{
				int neighbourBack = mMACGrid.GetIndexFromXYZ(x, y, z - 1);

				if (mMACGrid.GetCellType(cellIndex) != MACGrid::CellType::eSOLID)
				{
					interpXVelocities[neighbourBack] += mParticles[p].GetVelocity().x;
					interpYVelocities[neighbourBack] += mParticles[p].GetVelocity().y;
					interpZVelocities[neighbourBack] += mParticles[p].GetVelocity().z;
					++numParticlesToCell[neighbourBack];
				}
			}
		}
	}

	for (int c = 0; c < mMACGrid.GetNumCells(); c++)
	{
		if (numParticlesToCell[c] != 0.f)
		{
			interpXVelocities[c] *= 1.0f / numParticlesToCell[c];
			interpYVelocities[c] *= 1.0f / numParticlesToCell[c];
			interpZVelocities[c] *= 1.0f / numParticlesToCell[c];
		}
		
		if (mMACGrid.GetCellType(c) != MACGrid::CellType::eSOLID)
		{
			mMACGrid.SetCellXVelocity(c, interpXVelocities[c]);
			mMACGrid.SetCellYVelocity(c, interpYVelocities[c]);
			mMACGrid.SetCellZVelocity(c, interpZVelocities[c]);
		}
	}
}

void Fluid::InterpolateFromGrid()
{
	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(mParticles[p]);
		
		if (mMACGrid.GetCellType(cellIndex) != MACGrid::CellType::eSOLID)
		{
			int x, y, z;
			std::tie(x, y, z) = mMACGrid.GetXYZFromIndex(cellIndex);

			float velocityX = mMACGrid.GetCellXVelocity(cellIndex);
			float velocityY = mMACGrid.GetCellYVelocity(cellIndex);
			float velocityZ = mMACGrid.GetCellZVelocity(cellIndex);

			if (x > 0)
			{
				int neighbourLeft = mMACGrid.GetIndexFromXYZ(x - 1, y, z);

				if (mMACGrid.GetCellType(neighbourLeft) != MACGrid::CellType::eSOLID)
				{
					velocityX += mMACGrid.GetCellXVelocity(neighbourLeft);
					velocityX *= 0.5f;
				}
			}
			if (y > 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXYZ(x, y - 1, z);
				
				if (mMACGrid.GetCellType(neighbourBottom) != MACGrid::CellType::eSOLID)
				{
					velocityY += mMACGrid.GetCellYVelocity(neighbourBottom);
					velocityY *= 0.5f;
				}
			}
			if (z > 0)
			{
				int neighbourBack = mMACGrid.GetIndexFromXYZ(x, y, z - 1);
				
				if (mMACGrid.GetCellType(neighbourBack) != MACGrid::CellType::eSOLID)
				{
					velocityZ += mMACGrid.GetCellXVelocity(neighbourBack);
					velocityZ *= 0.5f;
				}
			}

			mParticles[p].SetVelocity(glm::vec3(velocityX, velocityY, velocityZ));
		}
		else
		{
			std::cout << "WARNING - particle may have penetrated solid boundary. \n";
			std::cout << "Fluid.pp - InterpolateFromGrid function. \n";
		}
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