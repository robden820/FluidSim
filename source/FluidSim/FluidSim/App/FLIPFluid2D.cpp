#include "FLIPFluid2D.h"
#include <iostream>

#include "oneapi/tbb.h"

#include "Interpolation.h"

FLIPFluid2D::FLIPFluid2D(const ApplicationData& inOutData)
{
	SeedParticles(inOutData);

	mFLIPBlend = inOutData.GetFLIPBlend();
}

void FLIPFluid2D::UpdateApplicationData(ApplicationData& inOutData)
{
	inOutData.Set2DParticlePositions(mParticlePositions);
	inOutData.SetNumParticles(mParticles.size());
}

void FLIPFluid2D::StepParticlesEuler(double deltaTime, const MACGrid2D& inMACGrid)
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

void FLIPFluid2D::StepParticles(double deltaTime, const MACGrid2D& inMACGrid)
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

		glm::dvec2 K1PrevVel = InterpolateFromGridCellBSplinePrev(inMACGrid, particlePosition, cellIndex);

		// Calculate K1 value.
		glm::dvec2 K1PIC = InterpolateFromGridCellBSpline(inMACGrid, particleIndex, cellIndex, SimulationType::ePIC);
		glm::dvec2 K1FLIP = K1PrevVel + InterpolateFromGridCellBSpline(inMACGrid, particleIndex, cellIndex, SimulationType::eFLIP);

		glm::dvec2 K1 = glm::mix(K1PIC, K1FLIP, mFLIPBlend);

		// Calculate K2 value.
		glm::dvec2 K2Pos = particlePosition + deltaTime * 0.5 * K1;
		int K2CellIndex = inMACGrid.GetClosestCell(K2Pos);

		glm::dvec2 K2PrevVel = InterpolateFromGridCellBSplinePrev(inMACGrid, K2Pos, K2CellIndex);

		glm::dvec2 K2PIC = InterpolateFromGridCellBSpline(inMACGrid, K2Pos, K2CellIndex, SimulationType::ePIC);
		glm::dvec2 K2FLIP = K2PrevVel + InterpolateFromGridCellBSpline(inMACGrid, K2Pos, K2CellIndex, SimulationType::eFLIP);

		glm::dvec2 K2 = glm::mix(K2PIC, K2FLIP, mFLIPBlend);

		// Calculate K3 value.
		glm::dvec2 K3Pos = particlePosition + deltaTime * 0.75 * K2;
		int K3CellIndex = inMACGrid.GetClosestCell(K3Pos);

		glm::dvec2 K3PrevVel = InterpolateFromGridCellBSplinePrev(inMACGrid, K3Pos, K3CellIndex);

		glm::dvec2 K3PIC = InterpolateFromGridCellBSpline(inMACGrid, K3Pos, K3CellIndex, SimulationType::ePIC);
		glm::dvec2 K3FLIP = K3PrevVel + InterpolateFromGridCellBSpline(inMACGrid, K3Pos, K3CellIndex, SimulationType::eFLIP);

		glm::dvec2 K3 = glm::mix(K3PIC, K3FLIP, mFLIPBlend);

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

void FLIPFluid2D::InterpolateToGrid(MACGrid2D& inMACGrid)
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

void FLIPFluid2D::InterpolateFromGrid(const MACGrid2D& inMACGrid)
{
	for (int particleIndex = 0; particleIndex < GetNumParticles(); particleIndex++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		if (cellIndex < 0 || inMACGrid.GetCellType(cellIndex) == CellType::eSOLID)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		// Interpolate new velocity for PIC.
		glm::dvec2 velocityPIC = InterpolateFromGridCell(inMACGrid, particleIndex, cellIndex, SimulationType::ePIC);

		// Interpolate difference for FLIP
		glm::dvec2 velocityDiffFLIP = InterpolateFromGridCell(inMACGrid, particleIndex, cellIndex, SimulationType::eFLIP);
		glm::dvec2 velocityFLIP = mParticles[particleIndex].GetVelocity() + velocityDiffFLIP;

		// Blend together to get final velocity
		glm::dvec2 finalVelocity = glm::mix(velocityPIC, velocityFLIP, mFLIPBlend);

		mParticles[particleIndex].SetVelocity(finalVelocity);
	}
}

void FLIPFluid2D::InterpolateFromGridBSpline(const MACGrid2D& inMACGrid)
{
	for (int particleIndex = 0; particleIndex < GetNumParticles(); particleIndex++)
	{
		int cellIndex = ClosestCellToParticle(inMACGrid, mParticles[particleIndex]);

		if (cellIndex < 0 || inMACGrid.GetCellType(cellIndex) == CellType::eSOLID)
		{
			// Invalid cell index returned if no closest cell found, i.e. particle has left grid.
			continue;
		}

		// Interpolate new velocity for PIC.
		glm::dvec2 velocityPIC = InterpolateFromGridCellBSpline(inMACGrid, particleIndex, cellIndex, SimulationType::ePIC);

		// Interpolate difference for FLIP
		glm::dvec2 velocityDiffFLIP = InterpolateFromGridCellBSpline(inMACGrid, particleIndex, cellIndex, SimulationType::eFLIP);
		glm::dvec2 velocityFLIP = mParticles[particleIndex].GetVelocity() + velocityDiffFLIP;

		// Blend together to get final velocity
		glm::dvec2 finalVelocity = glm::mix(velocityPIC, velocityFLIP, mFLIPBlend);

		mParticles[particleIndex].SetVelocity(finalVelocity);
	}
}

glm::dvec2 FLIPFluid2D::InterpolateFromGridCell(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex, SimulationType simType)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateFromGridCell(inMACGrid, particlePos, cellIndex, simType);
}

glm::dvec2 FLIPFluid2D::InterpolateFromGridCell(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex, SimulationType simType)
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

			if (simType == SimulationType::eFLIP) // If using FLIP interpolate the change in velocity.
			{
				xVelocities[j][i] = inMACGrid.GetCellXVelocityDiff(nearbyCellIndex);
				yVelocities[j][i] = inMACGrid.GetCellYVelocityDiff(nearbyCellIndex);
			}
			else // IF using PIC interpolate the new velocity.
			{
				xVelocities[j][i] = inMACGrid.GetCellXVelocity(nearbyCellIndex);
				yVelocities[j][i] = inMACGrid.GetCellYVelocity(nearbyCellIndex);
			}
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

glm::dvec2 FLIPFluid2D::InterpolateFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex, SimulationType simType)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateFromGridCellBSpline(inMACGrid, particlePos, cellIndex, simType);
}

glm::dvec2 FLIPFluid2D::InterpolateFromGridCellBSpline(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex, SimulationType simType)
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

			if (simType == SimulationType::eFLIP) // If using FLIP interpolate the change in velocity.
			{
				velocitySum.x += inMACGrid.GetCellXVelocityDiff(nearbyCellIndex) * k;
				velocitySum.y += inMACGrid.GetCellYVelocityDiff(nearbyCellIndex) * k;
			}
			else // If using PIC interpolate the new velocity.
			{
				velocitySum.x += inMACGrid.GetCellXVelocity(nearbyCellIndex) * k;
				velocitySum.y += inMACGrid.GetCellYVelocity(nearbyCellIndex) * k;
			}
				
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

glm::dvec2 FLIPFluid2D::InterpolateFromGridCellBSplinePrev(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();

	return InterpolateFromGridCellBSplinePrev(inMACGrid, particlePos, cellIndex);
}

glm::dvec2 FLIPFluid2D::InterpolateFromGridCellBSplinePrev(const MACGrid2D& inMACGrid, const glm::dvec2& particlePosition, int cellIndex)
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

			// Add half the cell size because the velocities are offset.
			glm::dvec2 diff = particlePosition - nearbyCellPos + inMACGrid.GetCellSize() * 0.5f;

			double k = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

			velocitySum.x += inMACGrid.GetCellXVelocityPrev(nearbyCellIndex) * k;
			velocitySum.y += inMACGrid.GetCellYVelocityPrev(nearbyCellIndex) * k;

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