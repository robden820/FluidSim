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
		glm::dvec2 cellPos = inMACGrid.GetCellCenter(cellIndex);

		glm::dvec2 velocity = mParticles[particleIndex].GetVelocity();

		glm::dvec2 newPos = particlePosition + velocity * deltaTime;

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
	std::vector<double> interpXMass;
	std::vector<double> interpYMass;

	int numCells = inMACGrid.GetNumCells();

	interpXVelocities.assign(numCells, 0.0);
	interpYVelocities.assign(numCells, 0.0);

	interpXMass.assign(numCells, 0.0);
	interpYMass.assign(numCells, 0.0);

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
				double cellSize = inMACGrid.GetCellSize();

				glm::dvec2 diffToX = (nearbyCellPos - (cellSize * 0.5, 0.0)) - particlePos;
				glm::dvec2 diffToY = (nearbyCellPos - (0.0, cellSize * 0.5)) - particlePos;
				
				// TO DO: different mass contribution weights for each face. SEE trilinear interpolating kernel.
				double weightX = InterpolateSupport(diffToX, inMACGrid.GetInverseCellSize());
				double weightY = InterpolateSupport(diffToY, inMACGrid.GetInverseCellSize());

				// Calculate contribution to cell mass
				double particleMass = mParticles[particleIndex].GetMass();
				interpXMass[nearbyCellIndex] += particleMass * weightX;
				interpYMass[nearbyCellIndex] += particleMass * weightY;

				// Calculate contribution to cell velocity.
				glm::dvec2 velocity = mParticles[particleIndex].GetVelocity();
				glm::dvec2 affineState = mParticles[particleIndex].GetAffineState();

				double velocityContributionX = particleMass * weightX * (velocity.x + affineState.x * diffToX.x);
				double velocityContributionY = particleMass * weightY * (velocity.y + affineState.y * diffToY.y);

				interpXVelocities[nearbyCellIndex] += velocityContributionX;
				interpYVelocities[nearbyCellIndex] += velocityContributionY;

				// Update contributed weights.
				contributedXWeights[nearbyCellIndex] += weightX;
				contributedYWeights[nearbyCellIndex] += weightY;
			}
		}
	}

	for (int c = 0; c < numCells; c++)
	{
		// Normalize all values by dividing by weight.
		if (contributedXWeights[c] != 0.0)
		{
			interpXVelocities[c] *= 1.0 / contributedXWeights[c];
			interpXMass[c] *= 1.0 / contributedXWeights[c];
		}
		if (contributedYWeights[c] != 0.0)
		{
			interpYVelocities[c] *= 1.0 / contributedYWeights[c];
			interpYMass[c] *= 1.0 / contributedYWeights[c];
		}

		// Update MACGrid values.
		if (interpXMass[c] != 0.0)
		{
			inMACGrid.SetIntXVelocity(c, interpXVelocities[c] / interpXMass[c]);
		}
		if (interpYMass[c] != 0.0)
		{
			inMACGrid.SetIntYVelocity(c, interpYVelocities[c] / interpYMass[c]);
		}
		// TO DO: else set velocity to zero?
		
		inMACGrid.SetCellMassX(c, interpXMass[c]);
		inMACGrid.SetCellMassY(c, interpYMass[c]);

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

		glm::dvec2 newAffineState = InterpolateAffineFromGridCellBSpline(inMACGrid, particleIndex, cellIndex);


		mParticles[particleIndex].SetVelocity(newVelocity);
		mParticles[particleIndex].SetAffineState(newAffineState);
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


glm::dvec2 APICFluid2D::InterpolateAffineFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();
	double particleMass = mParticles[particleIndex].GetMass();

	return InterpolateAffineFromGridCellBSpline(inMACGrid, particlePos, particleMass, cellIndex);
}

glm::dvec2 APICFluid2D::InterpolateAffineFromGridCellBSpline(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, const double particleMass, int cellIndex)
{
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	double affineSumX(0.0);
	double affineSumY(0.0);
	double totalWeight = 0.0;

	std::vector<std::vector<double>> weights;
	std::vector<double> zeros;
	zeros.assign(5, 0.0);
	weights.assign(5, zeros);


	// Interpolate to cells that may be close enough to particle.
	for (int i = x - 1; i <= x + 2; i++)
	{
		// If there isn't a cell where we are looking, continue.
		if (i < 0 || i > inMACGrid.GetNumCellsWidth() - 1)
		{
			continue;
		}

		for (int j = y - 1; j <= y + 2; j++)
		{
			// If there isn't a cell where we are looking, continue.
			if (j < 0 || j > inMACGrid.GetNumCellsHeight() - 1)
			{
				continue;
			}

			int nearbyCellIndex = inMACGrid.GetIndexFromXY(i, j);

			glm::dvec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);

			glm::dvec2 diff = particlePosition - nearbyCellPos;

			double weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

			int a = i - (x - 1);
			int b = j - (y - 1);

			weights[a][b] = weight;

			totalWeight += weight;
		}
	}

	// Interpolate to cells that may be close enough to particle.
	for (int i = x - 1; i <= x + 2; i++)
	{
		// If there isn't a cell where we are looking, continue.
		if (i < 0 || i > inMACGrid.GetNumCellsWidth() - 1)
		{
			continue;
		}

		for (int j = y - 1; j <= y + 2; j++)
		{
			// If there isn't a cell where we are looking, continue.
			if (j < 0 || j > inMACGrid.GetNumCellsHeight() - 1)
			{
				continue;
			}

			int nearbyCellIndex = inMACGrid.GetIndexFromXY(i, j);

			int a = i - (x - 1);
			int b = j - (y - 1);

			double deltaWeightX = weights[a + 1][b] - weights[a][b];
			double deltaWeightY = weights[a][b + 1] - weights[a][b];

			double invCellSize = inMACGrid.GetInverseCellSize();

			glm::dvec2 weightGrad(deltaWeightX * invCellSize, deltaWeightY * invCellSize);

			affineSumX += weightGrad.x * inMACGrid.GetCellXVelocity(nearbyCellIndex);
			affineSumY += weightGrad.y * inMACGrid.GetCellYVelocity(nearbyCellIndex);
		}
	}

	double affineStateX = affineSumX / totalWeight;
	double affineStateY = affineSumY / totalWeight;

	return glm::dvec2(affineStateX, affineStateY);
}