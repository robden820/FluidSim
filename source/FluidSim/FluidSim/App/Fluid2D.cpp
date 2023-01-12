#include "Fluid2D.h"
#include <iostream>

#include "oneapi/tbb.h"
#include "GLFW/glfw3.h"

Fluid2D::Fluid2D(int numParticles)
{
	float start = glfwGetTime();
	std::cout << "Initializing particles: ";
	// Initialize particles
	mParticles.reserve(numParticles);
	mParticlePositions.reserve(numParticles);

	for (int x = 0; x < 10; x++)
	{
		for (int y = 0; y < 10; y++)
		{
				glm::vec2 position((x + 10) * 0.5f, (y + 10) * 0.5f);

				Particle2D particle(position);

				mParticles.push_back(particle);
				mParticlePositions.push_back(position);
		}
	}

	std::cout << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	std::cout << "Initializing domain: ";

	// Initialize simulation Domain.
	Domain2D d(glm::vec2(10.0f, 10.0f), glm::vec2(10.0f, 10.0f));
	mDomain = d;

	std::cout << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	std::cout << "------------------------ \n";
	std::cout << "Initializing MAC Grid: \n";

	// Initialize Grid.
	mMACGridResolution = 50;
	MACGrid2D g(mDomain, mParticlePositions, mMACGridResolution);
	mMACGrid = g;

	std::cout << "Total : " << glfwGetTime() - start << "\n";
	std::cout << "------------------------ \n";
}

void Fluid2D::Update(ApplicationData& inOutData)
{
	float deltaTime = inOutData.GetDeltaTime();
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
		mParticlePositions[p] = mParticles[p].GetPosition();
	}

	// Ensure particles stay inside the simulation domain.
	for (int p = 0; p < GetNumParticles(); p++)
	{
		if (!mDomain.IsPointInDomain(mParticles[p].GetPosition()))
		{
			ClampParticleToDomain(mParticles[p]);
		}
	}

	inOutData.Set2DParticlePositions(mParticlePositions);
}

void Fluid2D::ClampParticleToDomain(Particle2D& particle)
{
	glm::vec2 particlePos = particle.GetPosition();
	glm::vec2 particleVel = particle.GetVelocity();

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

	particle.SetPosition(particlePos);
	particle.SetVelocity(particleVel);
}

void Fluid2D::InterpolateToGrid()
{
	for (int c = 0; c < mMACGrid.GetNumCells(); c++)
	{
		if (mMACGrid.GetCellType(c) != MACGrid::CellType::eSOLID)
		{
			mMACGrid.SetCellType(c, MACGrid::CellType::eAIR);
		}
	}

	std::vector<float> contributedXWeights;
	std::vector<float> contributedYWeights;

	std::vector<float> interpXVelocities;
	std::vector<float> interpYVelocities;

	interpXVelocities.assign(mMACGrid.GetNumCells(), 0.f);
	interpYVelocities.assign(mMACGrid.GetNumCells(), 0.f);

	contributedXWeights.assign(mMACGrid.GetNumCells(), 0.f);
	contributedYWeights.assign(mMACGrid.GetNumCells(), 0.f);

	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(mParticles[p]);

		glm::vec2 particlePos = mParticles[p].GetPosition();
		glm::vec2 cellPos = mMACGrid.GetCellCenter(cellIndex);

		float xDiff = particlePos.x - cellPos.x;
		float yDiff = particlePos.y - cellPos.y;

		float xWeight = (xDiff / mMACGrid.GetCellSize()) + 0.5f;
		float yWeight = (yDiff / mMACGrid.GetCellSize()) + 0.5f;

		float velocityX = mMACGrid.GetCellXVelocity(cellIndex) * xWeight;
		float velocityY = mMACGrid.GetCellYVelocity(cellIndex) * yWeight;

		if (mMACGrid.GetCellType(cellIndex) != MACGrid::CellType::eSOLID)
		{
			mMACGrid.SetCellType(cellIndex, MACGrid::CellType::eFLUID);

			interpXVelocities[cellIndex] += mParticles[p].GetVelocity().x * xWeight;
			interpYVelocities[cellIndex] += mParticles[p].GetVelocity().y * yWeight;

			contributedXWeights[cellIndex] += xWeight;
			contributedYWeights[cellIndex] += yWeight;

			int x, y;
			std::tie(x, y) = mMACGrid.GetXYFromIndex(cellIndex);

			if (x < 0)
			{
				int neighbourLeft = mMACGrid.GetIndexFromXY(x - 1, y);

				if (mMACGrid.GetCellType(cellIndex) != MACGrid::CellType::eSOLID)
				{
					interpXVelocities[neighbourLeft] += mParticles[p].GetVelocity().x * (1 - xWeight);
					interpYVelocities[neighbourLeft] += mParticles[p].GetVelocity().y * (1 - yWeight);

					contributedXWeights[neighbourLeft] += (1 - xWeight);
					contributedYWeights[neighbourLeft] += (1 - yWeight);
				}
			}
			if (y < 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXY(x, y - 1);

				if (mMACGrid.GetCellType(cellIndex) != MACGrid::CellType::eSOLID)
				{
					interpXVelocities[neighbourBottom] += mParticles[p].GetVelocity().x * (1 - xWeight);
					interpYVelocities[neighbourBottom] += mParticles[p].GetVelocity().y * (1 - yWeight);

					contributedXWeights[neighbourBottom] += (1 - xWeight);
					contributedYWeights[neighbourBottom] += (1 - yWeight);
				}
			}
		}
	}

	for (int c = 0; c < mMACGrid.GetNumCells(); c++)
	{
		if (contributedXWeights[c] != 0.f)
		{
			interpXVelocities[c] *= 1.0f / contributedXWeights[c];
		}
		if (contributedYWeights[c] != 0.f)
		{
			interpYVelocities[c] *= 1.0f / contributedYWeights[c];
		}

		if (mMACGrid.GetCellType(c) != MACGrid::CellType::eSOLID)
		{
			mMACGrid.SetCellXVelocity(c, interpXVelocities[c]);
			mMACGrid.SetCellYVelocity(c, interpYVelocities[c]);
		}
	}
}

void Fluid2D::InterpolateFromGrid()
{
	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(mParticles[p]);

		if (mMACGrid.GetCellType(cellIndex) != MACGrid::CellType::eSOLID)
		{
			int x, y;
			std::tie(x, y) = mMACGrid.GetXYFromIndex(cellIndex);

			glm::vec2 particlePos = mParticles[p].GetPosition();
			glm::vec2 cellPos = mMACGrid.GetCellCenter(cellIndex);

			float xDiff = particlePos.x - cellPos.x;
			float yDiff = particlePos.y - cellPos.y;

			float xWeight = (xDiff / mMACGrid.GetCellSize()) + 0.5f;
			float yWeight = (yDiff / mMACGrid.GetCellSize()) + 0.5f;

			float velocityX = mMACGrid.GetCellXVelocity(cellIndex) * xWeight;
			float velocityY = mMACGrid.GetCellYVelocity(cellIndex) * yWeight;

			if (x > 0)
			{
				int neighbourLeft = mMACGrid.GetIndexFromXY(x - 1, y);

				velocityX += mMACGrid.GetCellXVelocity(neighbourLeft);
				velocityX *= (1 - xWeight);
			}
			if (y > 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXY(x, y - 1);

				velocityY += mMACGrid.GetCellYVelocity(neighbourBottom);
				velocityY *= (1 - yWeight);
			}

			mParticles[p].SetVelocity(glm::vec2(velocityX, velocityY));
		}
		else
		{
			std::cout << "WARNING - particle may have penetrated solid boundary. \n";
			std::cout << "Fluid2D.pp - InterpolateFromGrid function. \n";
		}
	}
}

int Fluid2D::ClosestCellToParticle(const Particle2D& particle)
{
	glm::vec2 particlePos = particle.GetPosition();

	return mMACGrid.GetClosestCell(particlePos);
}