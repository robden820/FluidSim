#include "APICFluid2D.h"


#include "oneapi/tbb.h"

#include "Interpolation.h"

APICFluid2D::APICFluid2D(const ApplicationData& inOutData)
{
	SeedParticles(inOutData);
}

void APICFluid2D::UpdateApplicationData(ApplicationData& inOutData)
{
	inOutData.Set2DParticlePositions(mParticlePositions);
	inOutData.SetNumParticles(mParticles.size());
}

void APICFluid2D::StepParticlesEuler(double deltaTime, const MACGrid2D& inMACGrid)
{
	for (int particleIndex = 0; particleIndex < GetNumParticles(); particleIndex++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		if (cellIndex < 0)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		glm::dvec2 particlePosition = mParticles[particleIndex].GetPosition();

		glm::dvec2 vel = mParticles[particleIndex].GetVelocity();

		glm::dvec2 newPos = particlePosition + vel * deltaTime;

		mParticles[particleIndex].SetPosition(newPos);
		mParticlePositions[particleIndex] = mParticles[particleIndex].GetPosition();

		int finalCellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		// If the final cell is still valid but is a solid cell, project the particle outside of the solid.
		if (finalCellIndex >= 0 && inMACGrid.GetCellType(finalCellIndex) == CellType::eSOLID)
		{
			ProjectParticleToFluid(inMACGrid, particleIndex, finalCellIndex);
		}
	}
}

void APICFluid2D::StepParticles(double deltaTime, const MACGrid2D& inMACGrid)
{
	// RK3
	for (int particleIndex = 0; particleIndex < GetNumParticles(); particleIndex++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		if (cellIndex < 0)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		glm::dvec2 particlePosition = mParticles[particleIndex].GetPosition();

		// Calculate K1 value.
		glm::dvec2 K1 = InterpolateVelocityFromGridCellBSpline(inMACGrid, particleIndex, cellIndex);

		// Calculate K2 value.
		glm::dvec2 K2Pos = particlePosition + deltaTime * 0.5 * K1;
		int K2CellIndex = inMACGrid.GetClosestCell(K2Pos);

		glm::dvec2 K2 = InterpolateVelocityFromGridCellBSpline(inMACGrid, K2Pos, K2CellIndex);

		// Calculate K3 value.
		glm::dvec2 K3Pos = particlePosition + deltaTime * 0.75 * K2;
		int K3CellIndex = inMACGrid.GetClosestCell(K3Pos);

		glm::dvec2 K3 = InterpolateVelocityFromGridCellBSpline(inMACGrid, K3Pos, K3CellIndex);

		// Step particle
		mParticles[particleIndex].StepRK3(deltaTime, K1, K2, K3);

		mParticlePositions[particleIndex] = mParticles[particleIndex].GetPosition();

		int finalCellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		// If the final cell is still valid but is a solid cell, project the particle outside of the solid.
		if (finalCellIndex >= 0 && inMACGrid.GetCellType(finalCellIndex) == CellType::eSOLID)
		{
			ProjectParticleToFluid(inMACGrid, particleIndex, finalCellIndex);
		}
	}
}

void APICFluid2D::InterpolateToGrid(MACGrid2D& inMACGrid)
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

				if (inMACGrid.GetCellType(nearbyCellIndex) == CellType::eSOLID)
				{
					continue;
				}

				glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();
				glm::dvec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);

				glm::dvec2 diff = particlePos - nearbyCellPos + inMACGrid.GetCellSize() * 0.5;

				double k = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

				glm::dvec2 newVelocity = mParticles[particleIndex].GetVelocity() * k;

				interpXVelocities[nearbyCellIndex] += newVelocity.x;
				interpYVelocities[nearbyCellIndex] += newVelocity.y;

				contributedXWeights[nearbyCellIndex] += k;
				contributedYWeights[nearbyCellIndex] += k;
			}
		}
	}

	for (int c = 0; c < numCells; c++)
	{
		if (contributedXWeights[c] != 0.0)
		{
			interpXVelocities[c] *= 1.0 / contributedXWeights[c];
		}
		if (contributedYWeights[c] != 0.0)
		{
			interpYVelocities[c] *= 1.0 / contributedYWeights[c];
		}

		inMACGrid.SetIntXVelocity(c, interpXVelocities[c]);
		inMACGrid.SetIntYVelocity(c, interpYVelocities[c]);
	}
}

void APICFluid2D::InterpolateFromGrid(const MACGrid2D& inMACGrid)
{
	for (int particleIndex = 0; particleIndex < GetNumParticles(); particleIndex++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		if (cellIndex < 0 || inMACGrid.GetCellType(cellIndex) == CellType::eSOLID)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		glm::dvec2 newVelocity = InterpolateVelocityFromGridCell(inMACGrid, particleIndex, cellIndex);

		mParticles[particleIndex].SetVelocity(newVelocity);
	}
}

void APICFluid2D::InterpolateFromGridBSpline(const MACGrid2D& inMACGrid)
{
	for (int particleIndex = 0; particleIndex < GetNumParticles(); particleIndex++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		if (cellIndex < 0 || inMACGrid.GetCellType(cellIndex) == CellType::eSOLID)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		glm::dvec2 newVelocity = InterpolateVelocityFromGridCellBSpline(inMACGrid, particleIndex, cellIndex);
		double newAngularMomentum = InterpolateAngularFromGridCellBSpline(inMACGrid, particleIndex, cellIndex);

		mParticles[particleIndex].SetVelocity(newVelocity);
		mParticles[particleIndex].SetAngularMomentum(newAngularMomentum);
	}
}

glm::dvec2 APICFluid2D::InterpolateVelocityFromGridCell(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateVelocityFromGridCell(inMACGrid, particlePos, cellIndex);
}

glm::dvec2 APICFluid2D::InterpolateVelocityFromGridCell(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex)
{
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	double cellSize = inMACGrid.GetCellSize();

	// Initialise all velocities with z zero value in case we are close to the edge of the grid.
	double xVelocities[4][4] = { {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0} };
	double yVelocities[4][4] = { {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0} };

	// Find velocities for cells nearby where we are testing.
	for (int i = x - 1; i <= x + 2; i++)
	{
		if (i < 0 || i > inMACGrid.GetNumCellsWidth() - 1)
		{
			continue;
		}

		for (int j = y - 1; j <= y + 2; j++)
		{

			if (j < 0 || j > inMACGrid.GetNumCellsHeight() - 1)
			{
				continue;
			}

			int nearbyCellIndex = inMACGrid.GetIndexFromXY(i, j);

				xVelocities[j][i] = inMACGrid.GetCellXVelocity(nearbyCellIndex);
				yVelocities[j][i] = inMACGrid.GetCellYVelocity(nearbyCellIndex);
		}
	}

	// Put our particle position within the interval of a single grid cell.
	glm::dvec2 t = particlePosition - inMACGrid.GetCellCenter(cellIndex);
	t += cellSize * 0.5;

	// Interpolate to find new velocities.
	double xVel = Interpolation::BicubicInterpolate(xVelocities, cellSize, t);
	double yVel = Interpolation::BicubicInterpolate(yVelocities, cellSize, t);

	return glm::dvec2(xVel, yVel);
}

glm::dvec2 APICFluid2D::InterpolateVelocityFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateVelocityFromGridCellBSpline(inMACGrid, particlePos, cellIndex);
}

glm::dvec2 APICFluid2D::InterpolateVelocityFromGridCellBSpline(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex)
{
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	glm::dvec2 velocitySum(0.0, 0.0);
	double totalWeight = 0.0;

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

			if (inMACGrid.GetCellType(nearbyCellIndex) == CellType::eSOLID)
			{
				continue;
			}

			glm::dvec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);

			glm::dvec2 diff = particlePosition - nearbyCellPos + inMACGrid.GetCellSize() * 0.5f;

			double k = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

			velocitySum.x += inMACGrid.GetCellXVelocity(nearbyCellIndex) * k;
			velocitySum.y += inMACGrid.GetCellYVelocity(nearbyCellIndex) * k;

			totalWeight += k;
		}
	}

	return (velocitySum / (double)totalWeight);
}


double APICFluid2D::InterpolateAngularFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();
	double particleMass = mParticles[particleIndex].GetMass();

	return InterpolateAngularFromGridCellBSpline(inMACGrid, particlePos, particleMass, cellIndex);
}

double APICFluid2D::InterpolateAngularFromGridCellBSpline(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, const double particleMass, int cellIndex)
{
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	double angularSum = 0.0;
	double totalWeight = 0.0;

	double halfCell = inMACGrid.GetCellSize() * 0.5f;

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

			/*
			if (inMACGrid.GetCellType(nearbyCellIndex) == CellType::eSOLID)
			{
				continue;
			}
			*/

			glm::dvec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);

			// Add half the cell size because the velocities are offset. 
			glm::dvec2 diff = particlePosition - nearbyCellPos + halfCell;
			// TO DO: how to use the weight?
			double weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

			glm::dvec2 diffToX = particlePosition - nearbyCellPos + glm::dvec2(halfCell, 0.0);
			glm::dvec2 diffToY = particlePosition - nearbyCellPos + glm::dvec2(0.0, halfCell);

			double diffToXLen = glm::length(diffToX);
			double diffToYLen = glm::length(diffToY);
			
			// Calculat vectors of equal length to diff but perpendicular.
			glm::dvec2 perpX(1.0, -1.0 * diffToX.x / diffToX.y);
			glm::normalize(perpX);
			perpX *= diffToXLen * weight;

			glm::dvec2 perpY(1.0, -1.0 * diffToY.x / diffToY.y);
			glm::normalize(perpY);
			perpY *= diffToXLen * weight;

			glm::dvec2 velocity(inMACGrid.GetCellXVelocity(nearbyCellIndex), inMACGrid.GetCellYVelocity(nearbyCellIndex));
			glm::dvec2 momentum = particleMass * velocity;

			double angularMomentumX = glm::dot(perpX, momentum) / diffToXLen;
			double angularMomentumY = glm::dot(perpY, momentum) / diffToYLen;

			double signX = diffToX.x * inMACGrid.GetCellXVelocity(nearbyCellIndex);
			signX /= abs(signX);

			double signY = diffToY.y * inMACGrid.GetCellYVelocity(nearbyCellIndex);
			signY /= abs(signY);

			angularSum += angularMomentumX * signX;
			angularSum += angularMomentumY * signY;


			totalWeight += weight;
		}
	}

	return (angularSum / (double)totalWeight);
}