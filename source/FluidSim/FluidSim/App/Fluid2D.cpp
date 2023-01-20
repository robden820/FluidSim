#include "Fluid2D.h"
#include <iostream>

#include "oneapi/tbb.h"
#include "GLFW/glfw3.h"

Fluid2D::Fluid2D(ApplicationData& inOutData)
{
	float start = glfwGetTime();
	std::cout << "Initializing particles: ";
	// Initialize particles
	mParticles.reserve(inOutData.GetNumParticles());
	mParticlePositions.reserve(inOutData.GetNumParticles());

	//TO DO : fix this initialization to not rely on literals.
	for (int x = 0; x < 30; x++)
	{
		for (int y = 0; y < 30; y++)
		{
				glm::vec2 position((x - 5) * 0.1f, (y - 5) * 0.1f);

				Particle2D particle(position);

				mParticles.push_back(particle);
				mParticlePositions.push_back(position);
		}
	}

	inOutData.Set2DParticlePositions(mParticlePositions);

	std::cout << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	std::cout << "------------------------ \n";
	std::cout << "Initializing MAC Grid: \n";

	// Initialize Grid.
	MACGrid2D grid(inOutData);
	mMACGrid = grid;

	std::cout << "Total : " << glfwGetTime() - start << "\n";
	std::cout << "------------------------ \n";
}

void Fluid2D::Update(ApplicationData& inOutData)
{
	float deltaTime = inOutData.GetDeltaTime();
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
		mParticlePositions[p] = mParticles[p].GetPosition();
	}

	inOutData.Set2DParticlePositions(mParticlePositions);
}

void Fluid2D::InterpolateToGrid()
{
	for (int c = 0; c < mMACGrid.GetNumCells(); c++)
	{
		if (mMACGrid.GetCellType(c) != CellType::eSOLID)
		{
			mMACGrid.SetCellType(c, CellType::eAIR);
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

		if (cellIndex < 0)
		{
			continue;
		}

		glm::vec2 particlePos = mParticles[p].GetPosition();
		glm::vec2 cellPos = mMACGrid.GetCellCenter(cellIndex);

		float xDiff = particlePos.x - cellPos.x;
		float yDiff = particlePos.y - cellPos.y;

		float xWeight = (xDiff / mMACGrid.GetCellSize()) + 0.5f;
		float yWeight = (yDiff / mMACGrid.GetCellSize()) + 0.5f;

		if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
		{
			mMACGrid.SetCellType(cellIndex, CellType::eFLUID);

			interpXVelocities[cellIndex] += mParticles[p].GetVelocity().x * xWeight;
			interpYVelocities[cellIndex] += mParticles[p].GetVelocity().y * yWeight;

			contributedXWeights[cellIndex] += xWeight;
			contributedYWeights[cellIndex] += yWeight;

			int x, y;
			std::tie(x, y) = mMACGrid.GetXYFromIndex(cellIndex);

			if (x < 0)
			{
				int neighbourLeft = mMACGrid.GetIndexFromXY(x - 1, y);

				if (mMACGrid.GetCellType(neighbourLeft) != CellType::eSOLID)
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

				if (mMACGrid.GetCellType(neighbourBottom) != CellType::eSOLID)
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

		mMACGrid.SetCellXVelocity(c, interpXVelocities[c]);
		mMACGrid.SetCellYVelocity(c, interpYVelocities[c]);
	}
}

void Fluid2D::InterpolateFromGrid()
{
	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(mParticles[p]);

		if (cellIndex < 0)
		{
			continue;
		}

		if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
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

				velocityX += mMACGrid.GetCellXVelocity(neighbourLeft) * (1 - xWeight);
			}
			if (y > 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXY(x, y - 1);

				velocityY += mMACGrid.GetCellYVelocity(neighbourBottom) * (1 - yWeight);
			}

			// Handle bottom solid boundary
			if (y > 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXY(x, y - 1);

				if (mMACGrid.GetCellType(neighbourBottom) == CellType::eSOLID)
				{
					mParticles[p].ApplyForce(glm::vec2{ 0.f, 9.8f } * mParticles[p].GetMass());
				}
			}

			// Handle left and right solid boundaries
			if (x > 0)
			{
				int neighbourLeft = mMACGrid.GetIndexFromXY(x - 1, y);

				if (mMACGrid.GetCellType(neighbourLeft) == CellType::eSOLID)
				{
					velocityX = 0.f;
				}
			}
			if (x < mMACGrid.GetNumCellsWidth() - 1)
			{
				int neighbourRight = mMACGrid.GetIndexFromXY(x + 1, y);

				if (mMACGrid.GetCellType(neighbourRight) == CellType::eSOLID)
				{
					velocityX = 0.f;
				}
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