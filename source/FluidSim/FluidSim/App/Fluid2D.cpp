#include "Fluid2D.h"
#include <iostream>

#include "oneapi/tbb.h"

#include "Interpolation.h"

Fluid2D::Fluid2D(const ApplicationData& inOutData)
{
	SeedParticles(inOutData);
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

void Fluid2D::StepParticles(float deltaTime, const MACGrid2D& inMACGrid, float blend)
{
	// Make sure our blend value is between 0 and 1.
	std::clamp(blend, 0.0f, 1.0f);

	// RK3
	for (int particleIndex = 0; particleIndex < GetNumParticles(); particleIndex++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		if (cellIndex < 0)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		glm::vec2 particlePosition = mParticles[particleIndex].GetPosition();
		glm::dvec2 particleVelocity = mParticles[particleIndex].GetVelocity();

		// Calculate K1 value.
		glm::vec2 K1PIC = InterpolateFromGridCellPIC(inMACGrid, particleIndex, cellIndex);
		glm::vec2 K1FLIP = particleVelocity + InterpolateFromGridCellFLIP(inMACGrid, particleIndex, cellIndex);

		glm::vec2 K1 = (K1FLIP * blend) + (K1PIC * (1.0f - blend));

		// Calculate K2 value.
		glm::vec2 K2Pos = particlePosition + deltaTime * 0.5f * K1;
		int K2CellIndex = inMACGrid.GetClosestCell(K2Pos);

		glm::vec2 K2PIC = InterpolateFromGridCellPIC(inMACGrid, K2Pos, K2CellIndex);
		glm::vec2 K2FLIP = particleVelocity + InterpolateFromGridCellFLIP(inMACGrid, K2Pos, K2CellIndex);

		glm::vec2 K2 = (K2FLIP * blend) + (K2PIC * (1.0f - blend));

		// Calculate K3 value.
		glm::vec K3Pos = particlePosition + deltaTime * 0.75f * K2;
		int K3CellIndex = inMACGrid.GetClosestCell(K3Pos);

		glm::vec2 K3PIC = InterpolateFromGridCellPIC(inMACGrid, K3Pos, K3CellIndex);
		glm::vec2 K3FLIP = particleVelocity + InterpolateFromGridCellFLIP(inMACGrid, K3Pos, K3CellIndex);

		glm::vec2 K3 = (K3FLIP * blend) + (K3PIC * (1.0f - blend));

		// Step particle
		mParticles[particleIndex].StepRK3(deltaTime, K1, K2, K3);

		mParticlePositions[particleIndex] = mParticles[particleIndex].GetPosition();

		int finalCellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		if (finalCellIndex >= 0 && inMACGrid.GetCellType(finalCellIndex) == CellType::eSOLID)
		{
			ProjectParticleToFluid(inMACGrid, particleIndex, finalCellIndex);
		}
	}
}

void Fluid2D::ProjectParticleToFluid(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	// TO DO: project out using directions of level set.
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	if (x < 1)
	{
		int right = inMACGrid.GetIndexFromXY(x + 1, y);
		x++;
		mParticlePositions[particleIndex].x += inMACGrid.GetCellSize();
	}
	else if (x > inMACGrid.GetNumCellsWidth() - 2)
	{
		int left = inMACGrid.GetIndexFromXY(x - 1, y);
		x--;
		mParticlePositions[particleIndex].x -= inMACGrid.GetCellSize();
	}

	if (y < 1)
	{
		int top = inMACGrid.GetIndexFromXY(x, y + 1);
		y++;
		mParticlePositions[particleIndex].y += inMACGrid.GetCellSize();
	}
	else if (y > inMACGrid.GetNumCellsHeight() - 2)
	{
		int bottom = inMACGrid.GetIndexFromXY(x, y - 1);
		y--;
		mParticlePositions[particleIndex].y -= inMACGrid.GetCellSize();
	}

	mParticles[particleIndex].SetVelocity(glm::vec2(0.0f, 0.0f));
	mParticles[particleIndex].SetPosition(mParticlePositions[particleIndex]);
}

void Fluid2D::InterpolateToGrid(MACGrid2D& inMACGrid)
{
	std::vector<double> contributedXWeights;
	std::vector<double> contributedYWeights;

	std::vector<double> interpXVelocities;
	std::vector<double> interpYVelocities;

	int numCells = inMACGrid.GetNumCells();

	interpXVelocities.assign(numCells, 0.0);
	interpYVelocities.assign(numCells, 0.0);

	contributedXWeights.assign(numCells, 0.0);
	contributedYWeights.assign(numCells, 0.0);

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
		
		// Interpolate to cells that may be close enough to particle.
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

					glm::vec2 diff = particlePos - nearbyCellPos + inMACGrid.GetCellSize() * 0.5f;

					float k = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

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

float Fluid2D::InterpolateSupport(const glm::vec2& diff, float invCellSize)
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

void Fluid2D::InterpolateFromGrid(const MACGrid2D& inMACGrid, float blend)
{
	// Make sure our blend value is between 0 and 1.
	std::clamp(blend, 0.0f, 1.0f);

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
			// Interpolate new velocity for PIC.
			glm::vec2 velocityPIC = InterpolateFromGridCellPIC(inMACGrid, p, cellIndex);

			// Interpolate difference for FLIP
			glm::vec2 velocityDiffFLIP = InterpolateFromGridCellFLIP(inMACGrid, p, cellIndex);
			glm::vec2 velocityFLIP = mParticles[p].GetVelocity() + velocityDiffFLIP;

			// Blend together to get final velocity
			glm::dvec2 finalVelocity = (velocityFLIP * blend) + (velocityPIC * (1.0f - blend));
			mParticles[p].SetVelocity(finalVelocity);
		}
	}
}

void Fluid2D::InterpolateFromGridBSpline(const MACGrid2D& inMACGrid, float blend)
{
	// Make sure our blend value is between 0 and 1.
	std::clamp(blend, 0.0f, 1.0f);

	for (int particleIndex = 0; particleIndex < GetNumParticles(); particleIndex++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		if (cellIndex < 0 || inMACGrid.GetCellType(cellIndex) == CellType::eSOLID)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		// Interpolate new velocity for PIC.
		glm::vec2 velocityPIC = InterpolateFromGridCellBSplinePIC(inMACGrid, particleIndex, cellIndex);

		// Interpolate difference for FLIP
		glm::vec2 velocityDiffFLIP = InterpolateFromGridCellBSplineFLIP(inMACGrid, particleIndex, cellIndex);
		glm::vec2 velocityFLIP = mParticles[particleIndex].GetVelocity() + velocityDiffFLIP;

		// Blend together to get final velocity
		glm::dvec2 finalVelocity = (velocityFLIP * blend) + (velocityPIC * (1.0f - blend));
		mParticles[particleIndex].SetVelocity(finalVelocity);
	}
}

glm::dvec2 Fluid2D::InterpolateFromGridCellPIC(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::vec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateFromGridCellPIC(inMACGrid, particlePos, cellIndex);
}

glm::dvec2 Fluid2D::InterpolateFromGridCellPIC(const MACGrid2D& inMACGrid, const glm::vec2& particlePosition, int cellIndex)
{
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	double cellSize = inMACGrid.GetCellSize();

	double xVelocities[4][4];
	double yVelocities[4][4];

	// Initialise all velocities with z zero value in case we are close to the edge of the grid.
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			xVelocities[i][j] = 0.0;
			yVelocities[i][j] = 0.0;
		}
	}

	// Find velocities for cells nearby where we are testing.
	for (int i = 0; i < 4; i++)
	{
		int indexX = x + i - 1;

		if (indexX < 0 || indexX > inMACGrid.GetNumCellsWidth() - 1)
		{
			continue;
		}

		for (int j = 0; j < 4; j++)
		{
			int indexY = y + j - 1;

			if (indexY < 0 || indexY > inMACGrid.GetNumCellsHeight() - 1)
			{
				continue;
			}

			int nearbyCellIndex = inMACGrid.GetIndexFromXY(indexX, indexY);

			xVelocities[j][i] = inMACGrid.GetCellXVelocity(nearbyCellIndex);
			yVelocities[j][i] = inMACGrid.GetCellYVelocity(nearbyCellIndex);
		}
	}

	// Put our particle position within the interval of a single grid cell.
	glm::vec2 t = particlePosition - inMACGrid.GetCellCenter(cellIndex);
	t += cellSize * 0.5;

	// Interpolate to find new velocities.
	double xVel = Interpolation::BicubicInterpolate(xVelocities, cellSize, t);
	double yVel = Interpolation::BicubicInterpolate(yVelocities, cellSize, t);

	return glm::dvec2(xVel, yVel);
}

glm::dvec2 Fluid2D::InterpolateFromGridCellFLIP(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::vec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateFromGridCellFLIP(inMACGrid, particlePos, cellIndex);
}

glm::dvec2 Fluid2D::InterpolateFromGridCellFLIP(const MACGrid2D& inMACGrid, const glm::vec2& particlePosition, int cellIndex)
{
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	double cellSize = inMACGrid.GetCellSize();

	double xVelocitiesDiff[4][4];
	double yVelocitiesDiff[4][4];

	// Initialise all velocities with z zero value in case we are close to the edge of the grid.
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			xVelocitiesDiff[i][j] = 0.0;
			yVelocitiesDiff[i][j] = 0.0;
		}
	}

	// Find velocities for cells nearby where we are testing.
	for (int i = 0; i < 4; i++)
	{
		int indexX = x + i - 1;

		if (indexX < 0 || indexX > inMACGrid.GetNumCellsWidth() - 1)
		{
			continue;
		}

		for (int j = 0; j < 4; j++)
		{
			int indexY = y + j - 1;

			if (indexY < 0 || indexY > inMACGrid.GetNumCellsHeight() - 1)
			{
				continue;
			}

			int nearbyCellIndex = inMACGrid.GetIndexFromXY(indexX, indexY);

			xVelocitiesDiff[j][i] = inMACGrid.GetCellXVelocityDiff(nearbyCellIndex);
			yVelocitiesDiff[j][i] = inMACGrid.GetCellYVelocityDiff(nearbyCellIndex);
		}
	}

	// Put our particle position within the interval of a single grid cell.
	glm::vec2 t = particlePosition - inMACGrid.GetCellCenter(cellIndex);
	t += cellSize * 0.5;

	// Interpolate to find new velocities.
	double xVel = Interpolation::BicubicInterpolate(xVelocitiesDiff, cellSize, t);
	double yVel = Interpolation::BicubicInterpolate(yVelocitiesDiff, cellSize, t);

	return glm::dvec2(xVel, yVel);
}

glm::dvec2 Fluid2D::InterpolateFromGridCellBSplinePIC(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::vec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateFromGridCellBSplinePIC(inMACGrid, particlePos, cellIndex);
}

glm::dvec2 Fluid2D::InterpolateFromGridCellBSplinePIC(const MACGrid2D& inMACGrid, const glm::vec2& particlePosition, int cellIndex)
{
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	glm::dvec2 newVelocitySum(0.0, 0.0);
	float totalWeight = 0.0f;

	// Interpolate to cells that may be close enough to particle.
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
				glm::vec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);

				glm::vec2 diff = particlePosition - nearbyCellPos + inMACGrid.GetCellSize() * 0.5f;

				float k = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

				newVelocitySum.x += inMACGrid.GetCellXVelocity(nearbyCellIndex) * k;
				newVelocitySum.y += inMACGrid.GetCellYVelocity(nearbyCellIndex) * k;

				totalWeight += k;
			}
		}
	}

	return (newVelocitySum / (double)totalWeight);
}

glm::dvec2 Fluid2D::InterpolateFromGridCellBSplineFLIP(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::vec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateFromGridCellBSplineFLIP(inMACGrid, particlePos, cellIndex);
}

glm::dvec2 Fluid2D::InterpolateFromGridCellBSplineFLIP(const MACGrid2D& inMACGrid, const glm::vec2& particlePosition, int cellIndex)
{
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	glm::dvec2 diffVelocitySum(0.0, 0.0);
	float totalWeight = 0.0f;

	// Interpolate to cells that may be close enough to particle.
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
				glm::vec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);

				glm::vec2 diff = particlePosition - nearbyCellPos + inMACGrid.GetCellSize() * 0.5f;

				float k = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

				diffVelocitySum.x += inMACGrid.GetCellXVelocityDiff(nearbyCellIndex) * k;
				diffVelocitySum.y += inMACGrid.GetCellYVelocityDiff(nearbyCellIndex) * k;

				totalWeight += k;
			}
		}
	}

	glm::vec2 diffVelocity = diffVelocitySum / (double)totalWeight;

	return diffVelocity;
}

int Fluid2D::ClosestCellToParticle(const MACGrid2D& inMACGrid, const Particle2D& particle)
{
	glm::vec2 particlePos = particle.GetPosition();

	return inMACGrid.GetClosestCell(particlePos);
}

void Fluid2D::DeleteBadParticles(const MACGrid2D& inMACGrid)
{
	std::vector<Particle2D> safeParticles;
	std::vector<glm::vec2> safeParticlePositions;
	safeParticles.reserve(mParticles.size());
	safeParticlePositions.reserve(mParticlePositions.size());

	for (int pIndex = 0; pIndex < mParticles.size(); pIndex++)
	{
		if (ClosestCellToParticle(inMACGrid, mParticles[pIndex]) >= 0)
		{
			safeParticles.push_back(mParticles[pIndex]);
			safeParticlePositions.push_back(mParticlePositions[pIndex]);
		}
	}

	safeParticles.shrink_to_fit();
	safeParticlePositions.shrink_to_fit();

	mParticles = safeParticles;
	mParticlePositions = safeParticlePositions;
}