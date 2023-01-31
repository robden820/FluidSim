#include "Fluid2D.h"
#include <iostream>

#include "oneapi/tbb.h"
#include "GLFW/glfw3.h"

Fluid2D::Fluid2D(ApplicationData& inOutData)
{
	double start = glfwGetTime();
	std::cout << "Initializing particles: ";
	// Initialize particles
	mParticles.reserve(inOutData.GetNumParticles());
	mParticlePositions.reserve(inOutData.GetNumParticles());

	//TO DO : fix this initialization to not rely on literals.
	for (int x = 0; x < 30; x++)
	{
		for (int y = 0; y < 30; y++)
		{
				glm::vec2 position((x - 15) * 0.2f, (y - 15) * 0.2f);

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
	// Reset all non-solid cell types as air.
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
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
		{
			// Update cell type as it contains a fluid particle
			mMACGrid.SetCellType(cellIndex, CellType::eFLUID);

			// Interpolate the particle velocities to the grid cell and right hand neighbour.
			glm::vec2 particlePos = mParticles[p].GetPosition();
			glm::vec2 cellPos = mMACGrid.GetCellCenter(cellIndex);

			glm::vec2 diff = cellPos - particlePos;
			glm::vec2 weight = diff * mMACGrid.GetInverseCellSize() + 0.5f;

			interpXVelocities[cellIndex] += mParticles[p].GetVelocity().x * weight.x;
			interpYVelocities[cellIndex] += mParticles[p].GetVelocity().y * weight.y;

			contributedXWeights[cellIndex] += weight.x;
			contributedYWeights[cellIndex] += weight.y;

			int x, y;
			std::tie(x, y) = mMACGrid.GetXYFromIndex(cellIndex);

			if (x < mMACGrid.GetNumCellsWidth() - 1)
			{
				int neighbourRight = mMACGrid.GetIndexFromXY(x + 1, y);

				if (mMACGrid.GetCellType(neighbourRight) != CellType::eSOLID)
				{
					interpXVelocities[neighbourRight] += mParticles[p].GetVelocity().x * (1 - weight.x);

					contributedXWeights[neighbourRight] += (1 - weight.x);
				}
			}
			if (y < mMACGrid.GetNumCellsHeight() - 1)
			{
				int neighbourTop = mMACGrid.GetIndexFromXY(x, y + 1);

				if (mMACGrid.GetCellType(neighbourTop) != CellType::eSOLID)
				{
					interpYVelocities[neighbourTop] += mParticles[p].GetVelocity().y * (1 - weight.y);

					contributedYWeights[neighbourTop] += (1 - weight.y);
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
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
		{
			glm::vec2 particlePos = mParticles[p].GetPosition();
			glm::vec2 cellPos = mMACGrid.GetCellCenter(cellIndex);

			glm::vec2 diff = cellPos - particlePos;

			glm::dvec2 weight = diff * mMACGrid.GetInverseCellSize() + 0.5f;

			double velocityX = mMACGrid.GetCellXVelocity(cellIndex) * weight.x;
			double velocityY = mMACGrid.GetCellYVelocity(cellIndex) * weight.y;

			int x, y;
			std::tie(x, y) = mMACGrid.GetXYFromIndex(cellIndex);

			if (x < mMACGrid.GetNumCellsWidth() - 1)
			{
				int neighbourRight = mMACGrid.GetIndexFromXY(x + 1, y);

				velocityX += mMACGrid.GetCellXVelocity(neighbourRight) * (1 - weight.x);
			}
			else
			{
				// Normalize the velocity if weights don't sum to 1.
				velocityX *= 1 / weight.x;
			}
			if (y < mMACGrid.GetNumCellsHeight() - 1)
			{
				int neighbourTop = mMACGrid.GetIndexFromXY(x, y + 1);
				velocityY += mMACGrid.GetCellYVelocity(neighbourTop) * (1 - weight.y);
			}
			else
			{
				// Normalize the velocity if weights don't sum to 1.
				velocityY *= 1 / weight.y;
			}

			// Handle bottom solid boundary
			if (y > 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXY(x, y - 1);

				if (mMACGrid.GetCellType(neighbourBottom) == CellType::eSOLID)
				{
					mParticles[p].ApplyForce(glm::vec2{ 0.f, 9.8f } *mParticles[p].GetMass());
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
	}
}


int Fluid2D::ClosestCellToParticle(const Particle2D& particle)
{
	glm::vec2 particlePos = particle.GetPosition();

	return mMACGrid.GetClosestCell(particlePos);
}