#include "Fluid2D.h"
#include <iostream>

#include "oneapi/tbb.h"
#include "GLFW/glfw3.h"

Fluid2D::Fluid2D(const ApplicationData& inOutData)
{
	double start = glfwGetTime();
	std::cout << "Initializing particles: ";

	SeedParticles(inOutData);

	std::cout << "Total : " << glfwGetTime() - start << "\n";
	std::cout << "------------------------ \n";
}

void Fluid2D::SeedParticles(const ApplicationData& inOutData)
{
	float gridCellSize = inOutData.GetGridCellSize();
	float halfGridCell = gridCellSize * 0.5f;

	// To do: reserve memory space for particles array.

	for (int cellIndex = 0; cellIndex < inOutData.GetNumGridCells(); cellIndex++)
	{
		if (inOutData.GetCellType(cellIndex) == CellType::eFLUID)
		{
			// Seed 4 particles at each fluid cell.
			// To do : allow a user to set this value.

			for (int i = 0; i < 4; i++)
			{
				glm::vec2 particlePos = inOutData.GetCellCenter2D(cellIndex);

				particlePos.x += (rand() % 100) * 0.01f * gridCellSize + halfGridCell;
				particlePos.y += (rand() % 100) * 0.01f * gridCellSize + halfGridCell;

				Particle2D particle(particlePos);

				mParticles.push_back(particle);
				mParticlePositions.push_back(particlePos);
			}
		}
	}
}

void Fluid2D::UpdateApplicationData(ApplicationData& inOutData)
{
	inOutData.Set2DParticlePositions(mParticlePositions);
	inOutData.SetNumParticles(mParticles.size());
}

void Fluid2D::StepParticles(float deltaTime)
{
	for (int p = 0; p < GetNumParticles(); p++)
	{
		mParticles[p].StepParticle(deltaTime);
		mParticlePositions[p] = mParticles[p].GetPosition();
	}
}

void Fluid2D::StepParticlesRK3(float deltaTime, const MACGrid2D& inMACGrid)
{
	for (int particleIndex = 0; particleIndex < GetNumParticles(); particleIndex++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		if (cellIndex < 0)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		// Calculate K1 value.
		glm::vec2 K1 = InterpolateFromGridCell(inMACGrid, particleIndex, cellIndex);

		// Calculate K2 value.
		glm::vec2 K2Pos = mParticles[particleIndex].GetPosition() + deltaTime * 0.5f * K1;
		int K2CellIndex = inMACGrid.GetClosestCell(K2Pos);

		glm::vec2 K2 = InterpolateFromGridCell(inMACGrid, K2Pos, K2CellIndex);

		// Calculate K3 value.
		glm::vec K3Pos = mParticles[particleIndex].GetPosition() + deltaTime * 0.75f * K2;
		int K3CellIndex = inMACGrid.GetClosestCell(K3Pos);

		glm::vec2 K3 = InterpolateFromGridCell(inMACGrid, K3Pos, K3CellIndex);

		// Step particle
		mParticles[particleIndex].StepRK3(deltaTime, K1, K2, K3);

		mParticlePositions[particleIndex] = mParticles[particleIndex].GetPosition();
	}
}

void Fluid2D::InterpolateToGrid(MACGrid2D& inMACGrid)
{
	std::vector<double> contributedXWeights;
	std::vector<double> contributedYWeights;

	std::vector<double> interpXVelocities;
	std::vector<double> interpYVelocities;

	int numCells = inMACGrid.GetNumCells();

	interpXVelocities.assign(numCells, 0.f);
	interpYVelocities.assign(numCells, 0.f);

	contributedXWeights.assign(numCells, 0.f);
	contributedYWeights.assign(numCells, 0.f);

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

			glm::dvec2 diff = cellPos - particlePos;
			glm::dvec2 weight = diff * (double)inMACGrid.GetInverseCellSize() + 0.5;

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

	for (int c = 0; c < numCells; c++)
	{
		if (contributedXWeights[c] != 0.f)
		{
			interpXVelocities[c] *= 1.0f / contributedXWeights[c];
		}
		if (contributedYWeights[c] != 0.f)
		{
			interpYVelocities[c] *= 1.0f / contributedYWeights[c];
		}

		inMACGrid.SetIntXVelocity(c, interpXVelocities[c]);
		inMACGrid.SetIntYVelocity(c, interpYVelocities[c]);
	}
}

void Fluid2D::InterpolateToGridBSpline(MACGrid2D& inMACGrid)
{
	std::vector<double> contributedXWeights;
	std::vector<double> contributedYWeights;

	std::vector<double> interpXVelocities;
	std::vector<double> interpYVelocities;

	int numCells = inMACGrid.GetNumCells();

	interpXVelocities.assign(numCells, 0.f);
	interpYVelocities.assign(numCells, 0.f);

	contributedXWeights.assign(numCells, 0.f);
	contributedYWeights.assign(numCells, 0.f);

	for (int particleIndex = 0; particleIndex < GetNumParticles(); particleIndex++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		if (cellIndex < 0)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		int x, y;
		std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

		for (int i = x - 3; i <= x + 3; i++)
		{
			// If there isn't a cell where we are looking, continue.
			if (i < 0 || i > inMACGrid.GetNumCellsWidth() - 1)
			{
				continue;
			}

			for (int j = y - 3; j <= y + 3; j++)
			{
				// If there isn't a cell where we are looking, continue.
				if (j < 0 || j > inMACGrid.GetNumCellsHeight() - 1)
				{
					continue;
				}

				int nearbyCellIndex = inMACGrid.GetIndexFromXY(i, j);

				if (inMACGrid.GetCellType(nearbyCellIndex) != CellType::eSOLID)
				{
					glm::vec2 particlePos = mParticles[particleIndex].GetPosition();
					glm::vec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);

					glm::vec2 diff = particlePos - nearbyCellPos;

					float k = InterpolateToGridSupport(diff, inMACGrid.GetInverseCellSize());

					glm::vec2 newVelocity = mParticles[particleIndex].GetVelocity() * k;

					interpXVelocities[nearbyCellIndex] += newVelocity.x;
					interpYVelocities[nearbyCellIndex] += newVelocity.y;

					contributedXWeights[nearbyCellIndex] += k;
					contributedYWeights[nearbyCellIndex] += k;
				}
			}
		}
	}

	for (int c = 0; c < numCells; c++)
	{
		if (contributedXWeights[c] != 0.f)
		{
			interpXVelocities[c] *= 1.0f / contributedXWeights[c];
		}
		if (contributedYWeights[c] != 0.f)
		{
			interpYVelocities[c] *= 1.0f / contributedYWeights[c];
		}

		inMACGrid.SetIntXVelocity(c, interpXVelocities[c]);
		inMACGrid.SetIntYVelocity(c, interpYVelocities[c]);
	}
}

float Fluid2D::InterpolateToGridSupport(const glm::vec2& diff, float invCellSize)
{
	glm::vec2 scaled = diff * invCellSize;

	return BSpline(scaled.x) * BSpline(scaled.y);
}

float Fluid2D::BSpline(float input)
{
	float output = 0.0f;

	// If input is in [-1.5, 0.5)
	if (input >= -1.5f && input < -0.5f)
	{
		float a = input + 1.5f;
		output = 0.5f * a * a;
	}
	// if input is in [-0.5, 0.5)
	else if (input >= -0.5f && input < 0.5f)
	{
		output = 0.75f - input * input;
	}
	// If input is in [0.5, 1.5)
	else if (input >= 0.5f && input < 1.5f)
	{
		float a = 1.5f - input;
		output = 0.5f * a * a;
	}

	return output;
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
			glm::vec2 velocity = InterpolateFromGridCell(inMACGrid, p, cellIndex);
			
			mParticles[p].SetVelocity(velocity);
		}
	}
}

glm::dvec2 Fluid2D::InterpolateFromGridCell(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::vec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateFromGridCell(inMACGrid, particlePos, cellIndex);
}

glm::dvec2 Fluid2D::InterpolateFromGridCell(const MACGrid2D& inMACGrid, const glm::vec2& particlePosition, int cellIndex)
{
	glm::vec2 cellPos = inMACGrid.GetCellCenter(cellIndex);

	glm::vec2 diff = cellPos - particlePosition;

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
	if (y < inMACGrid.GetNumCellsHeight() - 1)
	{
		int neighbourTop = inMACGrid.GetIndexFromXY(x, y + 1);
		velocityY += inMACGrid.GetCellYVelocity(neighbourTop) * (1 - weight.y);
	}

	return glm::dvec2(velocityX, velocityY);
}

int Fluid2D::ClosestCellToParticle(const MACGrid2D& inMACGrid, const Particle2D& particle)
{
	glm::vec2 particlePos = particle.GetPosition();

	return inMACGrid.GetClosestCell(particlePos);
}