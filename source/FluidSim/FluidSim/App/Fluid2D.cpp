#include "Fluid2D.h"
#include <iostream>

#include "oneapi/tbb.h"
#include "GLFW/glfw3.h"

Fluid2D::Fluid2D(const ApplicationData& inOutData)
{
	double start = glfwGetTime();
	std::cout << "Initializing particles: ";
	// Initialize particles
	mParticles.reserve(inOutData.GetNumParticles());
	mParticlePositions = inOutData.Get2DParticlePositions();

	for (int pIndex = 0; pIndex < inOutData.GetNumParticles(); pIndex++)
	{
		Particle2D particle(mParticlePositions[pIndex]);
		mParticles.push_back(particle);
	}

	std::cout << "Total : " << glfwGetTime() - start << "\n";
	std::cout << "------------------------ \n";
}

void Fluid2D::UpdateApplicationData(ApplicationData& inOutData)
{
	inOutData.Set2DParticlePositions(mParticlePositions);
}

void Fluid2D::StepParticles(float deltaTime)
{
	for (int p = 0; p < GetNumParticles(); p++)
	{
		mParticles[p].StepParticle(deltaTime);
		mParticlePositions[p] = mParticles[p].GetPosition();
	}
}

void Fluid2D::InterpolateToGrid(MACGrid2D& inMACGrid)
{
	std::vector<float> contributedXWeights;
	std::vector<float> contributedYWeights;

	std::vector<float> interpXVelocities;
	std::vector<float> interpYVelocities;

	interpXVelocities.assign(inMACGrid.GetNumCells(), 0.f);
	interpYVelocities.assign(inMACGrid.GetNumCells(), 0.f);

	contributedXWeights.assign(inMACGrid.GetNumCells(), 0.f);
	contributedYWeights.assign(inMACGrid.GetNumCells(), 0.f);

	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[p]);

		if (cellIndex < 0)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		if (inMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
		{
			// Interpolate the particle velocities to the grid cell and right hand neighbour.
			glm::vec2 particlePos = mParticles[p].GetPosition();
			glm::vec2 cellPos = inMACGrid.GetCellCenter(cellIndex);

			glm::vec2 diff = cellPos - particlePos;
			glm::vec2 weight = diff * inMACGrid.GetInverseCellSize() + 0.5f;

			interpXVelocities[cellIndex] += mParticles[p].GetVelocity().x * weight.x;
			interpYVelocities[cellIndex] += mParticles[p].GetVelocity().y * weight.y;

			contributedXWeights[cellIndex] += weight.x;
			contributedYWeights[cellIndex] += weight.y;

			int x, y;
			std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

			if (x < inMACGrid.GetNumCellsWidth() - 1)
			{
				int neighbourRight = inMACGrid.GetIndexFromXY(x + 1, y);

				if (inMACGrid.GetCellType(neighbourRight) != CellType::eSOLID)
				{
					interpXVelocities[neighbourRight] += mParticles[p].GetVelocity().x * (1 - weight.x);

					contributedXWeights[neighbourRight] += (1 - weight.x);
				}
			}
			if (y < inMACGrid.GetNumCellsHeight() - 1)
			{
				int neighbourTop = inMACGrid.GetIndexFromXY(x, y + 1);

				if (inMACGrid.GetCellType(neighbourTop) != CellType::eSOLID)
				{
					interpYVelocities[neighbourTop] += mParticles[p].GetVelocity().y * (1 - weight.y);

					contributedYWeights[neighbourTop] += (1 - weight.y);
				}
			}

		}
	}

	for (int c = 0; c < inMACGrid.GetNumCells(); c++)
	{
		if (contributedXWeights[c] != 0.f)
		{
			interpXVelocities[c] *= 1.0f / contributedXWeights[c];
		}
		if (contributedYWeights[c] != 0.f)
		{
			interpYVelocities[c] *= 1.0f / contributedYWeights[c];
		}

		inMACGrid.SetCellXVelocity(c, interpXVelocities[c]);
		inMACGrid.SetCellYVelocity(c, interpYVelocities[c]);
	}
}

void Fluid2D::InterpolateFromGrid(const MACGrid2D& inMACGrid)
{
	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[p]);

		if (cellIndex < 0)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		if (inMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
		{
			glm::vec2 particlePos = mParticles[p].GetPosition();
			glm::vec2 cellPos = inMACGrid.GetCellCenter(cellIndex);

			glm::vec2 diff = cellPos - particlePos;

			glm::dvec2 weight = diff * inMACGrid.GetInverseCellSize() + 0.5f;

			double velocityX = inMACGrid.GetCellXVelocity(cellIndex) * weight.x;
			double velocityY = inMACGrid.GetCellYVelocity(cellIndex) * weight.y;

			int x, y;
			std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

			if (x < inMACGrid.GetNumCellsWidth() - 1)
			{
				int neighbourRight = inMACGrid.GetIndexFromXY(x + 1, y);

				velocityX += inMACGrid.GetCellXVelocity(neighbourRight) * (1 - weight.x);
			}
			else
			{
				// Normalize the velocity if weights don't sum to 1.
				velocityX *= 1 / weight.x;
			}
			if (y < inMACGrid.GetNumCellsHeight() - 1)
			{
				int neighbourTop = inMACGrid.GetIndexFromXY(x, y + 1);
				velocityY += inMACGrid.GetCellYVelocity(neighbourTop) * (1 - weight.y);
			}
			else
			{
				// Normalize the velocity if weights don't sum to 1.
				velocityY *= 1 / weight.y;
			}

			// Handle bottom solid boundary
			if (y > 0)
			{
				int neighbourBottom = inMACGrid.GetIndexFromXY(x, y - 1);

				if (inMACGrid.GetCellType(neighbourBottom) == CellType::eSOLID)
				{
					mParticles[p].ApplyForce(glm::vec2{ 0.f, 9.8f } *mParticles[p].GetMass());
				}
			}

			// Handle left and right solid boundaries
			if (x > 0)
			{
				int neighbourLeft = inMACGrid.GetIndexFromXY(x - 1, y);

				if (inMACGrid.GetCellType(neighbourLeft) == CellType::eSOLID)
				{
					velocityX = 0.0;
				}
			}
			if (x < inMACGrid.GetNumCellsWidth() - 1)
			{
				int neighbourRight = inMACGrid.GetIndexFromXY(x + 1, y);

				if (inMACGrid.GetCellType(neighbourRight) == CellType::eSOLID)
				{
					velocityX = 0.0;
				}
			}

			mParticles[p].SetVelocity(glm::vec2(velocityX, velocityY));
		}
	}
}


int Fluid2D::ClosestCellToParticle(const MACGrid2D& inMACGrid, const Particle2D& particle)
{
	glm::vec2 particlePos = particle.GetPosition();

	return inMACGrid.GetClosestCell(particlePos);
}