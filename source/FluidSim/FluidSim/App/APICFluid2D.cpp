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

		glm::dvec2 diff = cellPos - particlePosition;

		// Is this correct?
		glm::dvec3 linearVelocity = glm::cross(mParticles[particleIndex].GetAngularVelocity(), glm::dvec3(diff.x, diff.y, 0.0));

		glm::dvec2 vel = mParticles[particleIndex].GetVelocity() + glm::dvec2(linearVelocity.x, linearVelocity.y);

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
	std::vector<double> contributedMassWeights;

	std::vector<double> interpXVelocities;
	std::vector<double> interpYVelocities;
	std::vector<double> interpMass;

	int numCells = inMACGrid.GetNumCells();

	interpXVelocities.assign(numCells, 0.0);
	interpYVelocities.assign(numCells, 0.0);

	interpMass.assign(numCells, 0.0);

	contributedXWeights.assign(numCells, 0.0);
	contributedYWeights.assign(numCells, 0.0);
	contributedMassWeights.assign(numCells, 0.0);

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

				glm::dvec2 diffToCenter = particlePos - nearbyCellPos;
				glm::dvec2 diffToVelocities = diffToCenter + inMACGrid.GetCellSize() * 0.5;
				
				double cellWeight = InterpolateSupport(diffToCenter, inMACGrid.GetInverseCellSize());
				double velocityWeight = InterpolateSupport(diffToVelocities, inMACGrid.GetInverseCellSize());

				// Calculate contribution to cell mass
				double particleMass = mParticles[particleIndex].GetMass();
				interpMass[nearbyCellIndex] += mParticles[particleIndex].GetMass() * cellWeight;

				// Calculate contribution to cell velocity.
				glm::dvec3 angularVelocity = mParticles[particleIndex].GetAngularVelocity();
				glm::dvec3 linearVelocity = glm::cross(angularVelocity, glm::dvec3(-diffToCenter.x, -diffToCenter.y, 0.0));

				glm::dvec2 particleVelocity = mParticles[particleIndex].GetVelocity();

				glm::dvec3 velocityContribution = velocityWeight * particleMass * (glm::dvec3(particleVelocity.x, particleVelocity.y, 0.0) + linearVelocity);

				interpXVelocities[nearbyCellIndex] += velocityContribution.x;
				interpYVelocities[nearbyCellIndex] += velocityContribution.y;

				// Update contributed weights.
				contributedXWeights[nearbyCellIndex] += velocityWeight;
				contributedYWeights[nearbyCellIndex] += velocityWeight;
				contributedMassWeights[nearbyCellIndex] += cellWeight;
			}
		}
	}

	for (int c = 0; c < numCells; c++)
	{
		// Normalize all values by dividing by weight.
		if (contributedXWeights[c] != 0.0)
		{
			interpXVelocities[c] *= 1.0 / contributedXWeights[c];
		}
		if (contributedYWeights[c] != 0.0)
		{
			interpYVelocities[c] *= 1.0 / contributedYWeights[c];
		}
		if (contributedMassWeights[c] != 0.0)
		{
			interpMass[c] *= 1.0 / contributedMassWeights[c];
		}

		// Update MACGrid values.
		inMACGrid.SetIntXVelocity(c, interpXVelocities[c]);
		inMACGrid.SetIntYVelocity(c, interpYVelocities[c]);
		inMACGrid.SetCellMass(c, interpMass[c]);
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

		glm::dvec3 newAngularVelocity = InterpolateAngularFromGridCellBSpline(inMACGrid, particleIndex, cellIndex);

		mParticles[particleIndex].SetVelocity(newVelocity);
		mParticles[particleIndex].SetAngularVelocity(newAngularVelocity);
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


glm::dvec3 APICFluid2D::InterpolateAngularFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();
	double particleMass = mParticles[particleIndex].GetMass();

	return InterpolateAngularFromGridCellBSpline(inMACGrid, particlePos, particleMass, cellIndex);
}

glm::dvec3 APICFluid2D::InterpolateAngularFromGridCellBSpline(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, const double particleMass, int cellIndex)
{
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	glm::dvec3 angularSum(0.0);
	glm::dmat3 inertiaSum(0.0);
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

			/*
			if (inMACGrid.GetCellType(nearbyCellIndex) == CellType::eSOLID)
			{
				continue;
			}
			*/

			glm::dvec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);

			glm::dvec2 diff = nearbyCellPos - particlePosition;
			double weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

			// Vector 3
			glm::dvec3 weightedDiff(diff.x * weight, diff.y * weight, 0.0);

			// Calculate the angular momentum contributed by the cell.
			glm::dvec2 cellVelocity = inMACGrid.GetCellVelocity(nearbyCellIndex);
			// Vector 3
			glm::dvec3 cellMomentum(cellVelocity.x * particleMass, cellVelocity.y * particleMass, 0.0);

			angularSum += glm::cross(weightedDiff, cellMomentum);
			totalWeight += weight;

			// Calculate the cells contribution to the inertia tensor.
			glm::dmat3 crossMatDiff(0.0, 0.0, diff.y,
									0.0, 0.0, diff.x,
									-diff.y, diff.x, 0.0);
			glm::dmat3 transpose = glm::transpose(crossMatDiff);

			glm::dmat3 cellInertiaTensor = crossMatDiff * transpose;
			cellInertiaTensor *= weight * particleMass;

			inertiaSum += cellInertiaTensor;
		}
	}

	glm::dvec3 angularMomentum = angularSum / totalWeight;
	glm::dmat3 inertiaTensor = inertiaSum / totalWeight;

	glm::dvec3 angularVelocity = glm::inverse(inertiaTensor) * angularMomentum;

	return angularVelocity;
}