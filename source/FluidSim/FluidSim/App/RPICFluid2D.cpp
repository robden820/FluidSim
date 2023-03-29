#include "RPICFluid2D.h"

#include "oneapi/tbb.h"

#include "Interpolation.h"

RPICFluid2D::RPICFluid2D(const ApplicationData& inOutData)
{
	SeedParticles(inOutData);
}

void RPICFluid2D::UpdateApplicationData(ApplicationData& inOutData)
{
	inOutData.Set2DParticlePositions(mParticlePositions);
	inOutData.SetNumParticles(mParticles.size());
}

void RPICFluid2D::StepParticlesEuler(double deltaTime, const MACGrid2D& inMACGrid)
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

void RPICFluid2D::StepParticles(double deltaTime, const MACGrid2D& inMACGrid)
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
		glm::dvec2 nearbyCellPos = inMACGrid.GetCellCenter(cellIndex);

		glm::dvec2 diff = nearbyCellPos - particlePosition;
		glm::dvec3 diff3(diff, 0.0);

		double weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

		// Calculate K1 value.
		glm::dvec2 cellVel = InterpolateVelocityFromGridCellBSpline(inMACGrid, particleIndex, cellIndex);
		glm::dvec3 vel3(cellVel, 0.0);

		glm::dvec3 K1Vel = vel3 + glm::cross(mParticles[particleIndex].GetAngularVelocity(), diff3);
		glm::dvec2 K1(K1Vel.x, K1Vel.y);

		// Calculate K2 value.
		glm::dvec2 K2Pos = particlePosition + deltaTime * 0.5 * K1;
		int K2CellIndex = inMACGrid.GetClosestCell(K2Pos);
		glm::dvec2 K2(0.0);
		if (K2CellIndex >= 0)
		{
			diff = inMACGrid.GetCellCenter(K2CellIndex) - K2Pos;
			diff3 = glm::dvec3(diff, 0.0);

			weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

			cellVel = InterpolateVelocityFromGridCellBSpline(inMACGrid, K2Pos, K2CellIndex);
			vel3 = glm::dvec3(cellVel, 0.0);

			glm::dvec3 K2Vel = vel3 + glm::cross(mParticles[particleIndex].GetAngularVelocity(), diff3);
			K2 = glm::dvec2(K2Vel.x, K2Vel.y);
		}

		// Calculate K3 value.
		glm::dvec2 K3Pos = particlePosition + deltaTime * 0.75 * K2;
		int K3CellIndex = inMACGrid.GetClosestCell(K3Pos);
		glm::dvec2 K3(0.0);
		if (K3CellIndex >= 0)
		{
			diff = inMACGrid.GetCellCenter(K3CellIndex) - K3Pos;
			diff3 = glm::dvec3(diff, 0.0);

			weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

			cellVel = InterpolateVelocityFromGridCellBSpline(inMACGrid, K3Pos, K3CellIndex);
			vel3 = glm::dvec3(cellVel, 0.0);

			glm::dvec3 K3Vel = vel3 + glm::cross(mParticles[particleIndex].GetAngularVelocity(), diff3);
			K3 = glm::dvec2(K3Vel.x, K3Vel.y);
		}
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

void RPICFluid2D::InterpolateToGrid(MACGrid2D& inMACGrid)
{
	std::vector<double> contributedWeights;

	std::vector<double> interpXVelocities;
	std::vector<double> interpYVelocities;
	std::vector<double> interpMass;

	int numCells = inMACGrid.GetNumCells();

	interpXVelocities.assign(numCells, 0.0);
	interpYVelocities.assign(numCells, 0.0);

	interpMass.assign(numCells, 0.0);

	contributedWeights.assign(numCells, 0.0);

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

				if (inMACGrid.GetCellType(nearbyCellIndex) == CellType::eSOLID)
				{
					continue;
				}

				glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();
				glm::dvec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);

				glm::dvec2 diff = nearbyCellPos - particlePos;
				glm::dvec3 diff3(diff, 0.0);

				double weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

				// Calculate contribution to cell mass
				double particleMass = mParticles[particleIndex].GetMass();
				interpMass[nearbyCellIndex] += particleMass * weight;

				// Calculate contribution to cell velocity.
				glm::dvec2 particleVel = mParticles[particleIndex].GetVelocity();
				glm::dvec3 vel3(particleVel, 0.0);

				glm::dvec3 velocityContribution(0.0);

				velocityContribution = weight * particleMass * (vel3 + glm::cross(mParticles[particleIndex].GetAngularVelocity(), diff3));

				interpXVelocities[nearbyCellIndex] += velocityContribution.x;
				interpYVelocities[nearbyCellIndex] += velocityContribution.y;

				// Update contributed weights.
				contributedWeights[nearbyCellIndex] += weight;
			}
		}
	}

	for (int c = 0; c < numCells; c++)
	{
		// Normalize all values by dividing by weight.
		if (contributedWeights[c] != 0.0)
		{
			interpXVelocities[c] *= 1.0 / contributedWeights[c];
			interpYVelocities[c] *= 1.0 / contributedWeights[c];
			interpMass[c] *= 1.0 / contributedWeights[c];
		}

		// Update MACGrid values.
		if (interpMass[c] != 0.0)
		{
			inMACGrid.SetIntXVelocity(c, interpXVelocities[c] / interpMass[c]);
			inMACGrid.SetIntYVelocity(c, interpYVelocities[c] / interpMass[c]);
		}

		inMACGrid.SetCellMassX(c, interpMass[c]);
		inMACGrid.SetCellMassY(c, interpMass[c]);

	}
}

void RPICFluid2D::InterpolateFromGrid(const MACGrid2D& inMACGrid)
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

void RPICFluid2D::InterpolateFromGridBSpline(const MACGrid2D& inMACGrid)
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

glm::dvec2 RPICFluid2D::InterpolateVelocityFromGridCell(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateVelocityFromGridCell(inMACGrid, particlePos, cellIndex);
}

glm::dvec2 RPICFluid2D::InterpolateVelocityFromGridCell(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex)
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

glm::dvec2 RPICFluid2D::InterpolateVelocityFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateVelocityFromGridCellBSpline(inMACGrid, particlePos, cellIndex);
}

glm::dvec2 RPICFluid2D::InterpolateVelocityFromGridCellBSpline(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex)
{
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	glm::dvec2 velocitySum(0.0, 0.0);
	double totalWeight = 0.0;

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

	if (totalWeight > 0.0)
	{
		velocitySum /= totalWeight;
	}
	else
	{
		velocitySum = glm::dvec2(0.0);
	}

	return velocitySum;
}

glm::dvec3 RPICFluid2D::InterpolateAngularFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();
	double particleMass = mParticles[particleIndex].GetMass();

	return InterpolateAngularFromGridCellBSpline(inMACGrid, particlePos, particleMass, cellIndex);
}

glm::dvec3 RPICFluid2D::InterpolateAngularFromGridCellBSpline(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, double particleMass, int cellIndex)
{
	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	glm::dvec3 angularMomentumSum(0.0);
	glm::dmat3 inertiaSum(0.0);
	double totalWeight = 0.0;


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

			if (inMACGrid.GetCellType(nearbyCellIndex) == CellType::eSOLID)
			{
				continue;
			}

			glm::dvec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);
			glm::dvec2 nearbyCellVelocity(inMACGrid.GetCellXVelocity(nearbyCellIndex), inMACGrid.GetCellYVelocity(nearbyCellIndex));

			glm::dvec2 diff = nearbyCellPos - particlePosition;

			double weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

			// Vector 3 versions of diff and velocity.
			glm::dvec3 diff3(diff, 0.0);
			glm::dvec3 velocity3(nearbyCellVelocity, 0.0);

			// Calculate contributions to angular momentum.
			glm::dvec3 angularMomentum = weight * particleMass * glm::cross(diff3, velocity3);

			// Calculate contributions to inertia.
			const glm::dmat3 crossProductMat(0.0, 0.0, diff.y,
											 0.0, 0.0, -diff.x,
											-diff.y, diff.x, 0.0);
			const glm::dmat3 crossProductMatT(0.0, 0.0, -diff.y,
											  0.0, 0.0, diff.x,
											  diff.y, -diff.x, 0.0);

			glm::dmat3 inertia = weight * particleMass * crossProductMat * crossProductMatT;

			// Add to sum.
			inertiaSum += inertia;
			angularMomentumSum += angularMomentum;

			totalWeight += weight;
		}
	}

	// Normalize the sums incase weight doesn't add to 1.
	if (totalWeight > 0.0)
	{
		angularMomentumSum /= totalWeight;
		inertiaSum /= totalWeight;
	}

	glm::dvec3 angularVelocity = glm::inverse(inertiaSum) * angularMomentumSum;

	return angularVelocity;
}