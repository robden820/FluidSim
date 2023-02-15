#include "Fluid3D.h"
#include <iostream>

#include "oneapi/tbb.h"
#include "GLFW/glfw3.h"

Fluid3D::Fluid3D(const ApplicationData& inData)
{
	double start = glfwGetTime();
	std::cout << "Initializing particles: ";
	// Initialize particles
	mParticles.reserve(inData.GetNumParticles());
	mParticlePositions.reserve(inData.GetNumParticles());

	// TO DO: fix particle initialisation.
	for (int x = 0; x < inData.GetNumGridCellsWidth(); x++)
	{
		for (int y = 0; y < inData.GetNumGridCellsHeight(); y++)
		{
			for (int z = 0; z < inData.GetNumGridCellsLength(); z++)
			{
				glm::dvec3 position((x + 10) * 0.5f, (y + 10) * 0.5f, (z + 10) * 0.5f);

				Particle3D particle(position);

				mParticles.push_back(particle);
				mParticlePositions.push_back(position);
			}
		}
	}

	std::cout << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	std::cout << "------------------------ \n";
	std::cout << "Initializing MAC Grid: \n";

	// Initialize Grid.
	MACGrid3D grid(inData);
	mMACGrid = grid;

	std::cout <<"Total : " << glfwGetTime() - start << "\n";
	std::cout << "------------------------ \n";
}

void Fluid3D::Update(ApplicationData& inOutData)
{
	double deltaTime = inOutData.GetDeltaTime();

	// Transfer particle velocities to grid.
	InterpolateToGrid();

	// Update grid velocities.
	mMACGrid.Update(inOutData);

	// Interpolate velocities back to particles.
	InterpolateFromGrid();

	// Advect particles.
	for (int p = 0; p < GetNumParticles(); p++)
	{
		mParticles[p].StepParticle(deltaTime);
	}

	/* TO DO: shouldn't be necessary
	// Ensure particles stay inside the simulation domain.
	for (int p = 0; p < GetNumParticles(); p++)
	{
		if (!mDomain.IsPointInDomain(mParticles[p].GetPosition()))
		{
			ClampParticleToDomain(mParticles[p]);
		}
	}
	*/

	inOutData.Set3DParticlePositions(mParticlePositions);
}

/*
void Fluid3D::ClampParticleToDomain(Particle3D& particle)
{
	glm::dvec3 particlePos = particle.GetPosition();
	glm::dvec3 particleVel = particle.GetVelocity();

	if (particlePos.x < mDomain.GetLeft())
	{
		particlePos.x = mDomain.GetLeft();
		particleVel.x = 0.0;
	}
	else if (particlePos.x > mDomain.GetRight())
	{
		particlePos.x = mDomain.GetRight();
		particleVel.x = 0.0;
	}

	if (particlePos.y < mDomain.GetBottom())
	{
		particlePos.y = mDomain.GetBottom();
		particleVel.y = 0.0;
	}
	else if (particlePos.y > mDomain.GetTop())
	{
		particlePos.y = mDomain.GetTop();
		particleVel.y = 0.0;
	}

	if (particlePos.z < mDomain.GetBack())
	{
		particlePos.z = mDomain.GetBack();
		particleVel.z = 0.0;
	}
	else if (particlePos.z > mDomain.GetFront())
	{
		particlePos.z = mDomain.GetFront();
		particleVel.z = 0.0;
	}

	particle.SetPosition(particlePos);
	particle.SetVelocity(particleVel);
}*/

void Fluid3D::InterpolateToGrid()
{
	for (int c = 0; c < mMACGrid.GetNumCells(); c++)
	{
		if (mMACGrid.GetCellType(c) != CellType::eSOLID)
		{
			mMACGrid.SetCellType(c, CellType::eAIR);
		}
	}

	std::vector<double> contributedXWeights;
	std::vector<double> contributedYWeights;
	std::vector<double> contributedZWeights;

	std::vector<double> interpXVelocities;
	std::vector<double> interpYVelocities;
	std::vector<double> interpZVelocities;

	interpXVelocities.assign(mMACGrid.GetNumCells(), 0.f);
	interpYVelocities.assign(mMACGrid.GetNumCells(), 0.f);
	interpZVelocities.assign(mMACGrid.GetNumCells(), 0.f);

	contributedXWeights.assign(mMACGrid.GetNumCells(), 0.f);
	contributedYWeights.assign(mMACGrid.GetNumCells(), 0.f);
	contributedZWeights.assign(mMACGrid.GetNumCells(), 0.f);

	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(mParticles[p]);

		glm::dvec3 particlePos = mParticles[p].GetPosition();
		glm::dvec3 cellPos = mMACGrid.GetCellCenter(cellIndex);

		double xDiff = particlePos.x - cellPos.x;
		double yDiff = particlePos.y - cellPos.y;
		double zDiff = particlePos.z - cellPos.z;

		double xWeight = (xDiff / mMACGrid.GetCellSize()) + 0.5f;
		double yWeight = (yDiff / mMACGrid.GetCellSize()) + 0.5f;
		double zWeight = (zDiff / mMACGrid.GetCellSize()) + 0.5f;

		double velocityX = mMACGrid.GetCellXVelocity(cellIndex) * xWeight;
		double velocityY = mMACGrid.GetCellYVelocity(cellIndex) * yWeight;
		double velocityZ = mMACGrid.GetCellZVelocity(cellIndex) * zWeight;

		if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
		{
			mMACGrid.SetCellType(cellIndex, CellType::eFLUID);

			interpXVelocities[cellIndex] += mParticles[p].GetVelocity().x * xWeight;
			interpYVelocities[cellIndex] += mParticles[p].GetVelocity().y * yWeight;
			interpZVelocities[cellIndex] += mParticles[p].GetVelocity().z * zWeight;
			
			contributedXWeights[cellIndex] += xWeight;
			contributedYWeights[cellIndex] += yWeight;
			contributedZWeights[cellIndex] += zWeight;

			int x, y, z;
			std::tie(x, y, z) = mMACGrid.GetXYZFromIndex(cellIndex);

			if (x < 0)
			{
				int neighbourLeft = mMACGrid.GetIndexFromXYZ(x - 1, y, z);

				if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
				{
					interpXVelocities[neighbourLeft] += mParticles[p].GetVelocity().x * (1 - xWeight);
					interpYVelocities[neighbourLeft] += mParticles[p].GetVelocity().y * (1 - yWeight);
					interpZVelocities[neighbourLeft] += mParticles[p].GetVelocity().z * (1 - zWeight);
					
					contributedXWeights[neighbourLeft] += (1 - xWeight);
					contributedYWeights[neighbourLeft] += (1 - yWeight);
					contributedZWeights[neighbourLeft] += (1 - zWeight);
				}
			}
			if (y < 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXYZ(x, y - 1, z);

				if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
				{
					interpXVelocities[neighbourBottom] += mParticles[p].GetVelocity().x * (1 - xWeight);
					interpYVelocities[neighbourBottom] += mParticles[p].GetVelocity().y * (1 - yWeight);
					interpZVelocities[neighbourBottom] += mParticles[p].GetVelocity().z * (1 - zWeight);
					
					contributedXWeights[neighbourBottom] += (1 - xWeight);
					contributedYWeights[neighbourBottom] += (1 - yWeight);
					contributedZWeights[neighbourBottom] += (1 - zWeight);
				}
			}
			if (z < 0)
			{
				int neighbourBack = mMACGrid.GetIndexFromXYZ(x, y, z - 1);

				if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
				{
					interpXVelocities[neighbourBack] += mParticles[p].GetVelocity().x * (1 - xWeight);
					interpYVelocities[neighbourBack] += mParticles[p].GetVelocity().y * (1 - yWeight);
					interpZVelocities[neighbourBack] += mParticles[p].GetVelocity().z * (1 - zWeight);
					
					contributedXWeights[neighbourBack] += (1 - xWeight);
					contributedYWeights[neighbourBack] += (1 - yWeight);
					contributedZWeights[neighbourBack] += (1 - zWeight);
				}
			}
		}
	}

	for (int c = 0; c < mMACGrid.GetNumCells(); c++)
	{
		if (contributedXWeights[c] != 0.f)
		{
			interpXVelocities[c] *= 1.0 / contributedXWeights[c];
		}
		if (contributedYWeights[c] != 0.f)
		{
			interpYVelocities[c] *= 1.0 / contributedYWeights[c];
		}
		if (contributedZWeights[c] != 0.f)
		{
			interpZVelocities[c] *= 1.0 / contributedZWeights[c];
		}
		
		if (mMACGrid.GetCellType(c) != CellType::eSOLID)
		{
			mMACGrid.SetCellXVelocity(c, interpXVelocities[c]);
			mMACGrid.SetCellYVelocity(c, interpYVelocities[c]);
			mMACGrid.SetCellZVelocity(c, interpZVelocities[c]);
		}
	}
}

void Fluid3D::InterpolateFromGrid()
{
	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(mParticles[p]);
		
		if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
		{
			int x, y, z;
			std::tie(x, y, z) = mMACGrid.GetXYZFromIndex(cellIndex);

			glm::dvec3 particlePos = mParticles[p].GetPosition();
			glm::dvec3 cellPos = mMACGrid.GetCellCenter(cellIndex);

			double xDiff = particlePos.x - cellPos.x;
			double yDiff = particlePos.y - cellPos.y;
			double zDiff = particlePos.z - cellPos.z;

			double xWeight = (xDiff / mMACGrid.GetCellSize()) + 0.5f;
			double yWeight = (yDiff / mMACGrid.GetCellSize()) + 0.5f;
			double zWeight = (zDiff / mMACGrid.GetCellSize()) + 0.5f;

			double velocityX = mMACGrid.GetCellXVelocity(cellIndex) * xWeight;
			double velocityY = mMACGrid.GetCellYVelocity(cellIndex) * yWeight;
			double velocityZ = mMACGrid.GetCellZVelocity(cellIndex) * zWeight;

			if (x > 0)
			{
				int neighbourLeft = mMACGrid.GetIndexFromXYZ(x - 1, y, z);

				velocityX += mMACGrid.GetCellXVelocity(neighbourLeft);
				velocityX *= (1 - xWeight);
			}
			if (y > 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXYZ(x, y - 1, z);

				velocityY += mMACGrid.GetCellYVelocity(neighbourBottom);
				velocityY *= (1 - yWeight);
			}
			if (z > 0)
			{
				int neighbourBack = mMACGrid.GetIndexFromXYZ(x, y, z - 1);
				
				velocityZ += mMACGrid.GetCellZVelocity(neighbourBack);
				velocityZ *= (1 - xWeight);
			}

			mParticles[p].SetVelocity(glm::dvec3(velocityX, velocityY, velocityZ));
		}
		else
		{
			std::cout << "WARNING - particle may have penetrated solid boundary. \n";
			std::cout << "Fluid3D.pp - InterpolateFromGrid function. \n";
		}
	}
}

int Fluid3D::ClosestCellToParticle(const Particle3D& particle)
{
	glm::dvec3 particlePos = particle.GetPosition();

	return mMACGrid.GetClosestCell(particlePos);
}