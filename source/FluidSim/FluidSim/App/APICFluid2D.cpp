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

	std::vector<double> speeds;
	speeds.reserve(mParticles.size());

	for (int i = 0; i < mParticles.size(); i++)
	{
		speeds.push_back(glm::length(mParticles[i].GetVelocity()));
	}

	inOutData.SetParticleSpeeds(speeds);
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
		glm::dvec2 nearbyCellPos = inMACGrid.GetCellCenter(cellIndex);

		glm::dvec2 diff = nearbyCellPos - particlePosition;

		double weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

		// Calculate K1 value.
		glm::dvec2 K1(0.0);

		if (mParticles[particleIndex].HasValidInertiaAnalog())
		{
			K1 = InterpolateVelocityFromGridCellBSpline(inMACGrid, particleIndex, cellIndex) + (mParticles[particleIndex].GetAffineState() * glm::inverse(mParticles[particleIndex].GetInertiaAnalog())) * diff;
		}
		else
		{
			K1 = InterpolateVelocityFromGridCellBSpline(inMACGrid, particleIndex, cellIndex);
		}

		// Calculate K2 value.
		glm::dvec2 K2Pos = particlePosition + deltaTime * 0.5 * K1;
		int K2CellIndex = inMACGrid.GetClosestCell(K2Pos);
		glm::dvec2 K2(0.0);
		if (K2CellIndex >= 0)
		{
			diff = inMACGrid.GetCellCenter(K2CellIndex) - K2Pos;

			weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

			if (mParticles[particleIndex].HasValidInertiaAnalog())
			{
				K2 = InterpolateVelocityFromGridCellBSpline(inMACGrid, K2Pos, K2CellIndex) + (mParticles[particleIndex].GetAffineState() * glm::inverse(mParticles[particleIndex].GetInertiaAnalog())) * diff;
			}
			else
			{
				K2 = InterpolateVelocityFromGridCellBSpline(inMACGrid, K2Pos, K2CellIndex);
			}
		}

		// Calculate K3 value.
		glm::dvec2 K3Pos = particlePosition + deltaTime * 0.75 * K2;
		int K3CellIndex = inMACGrid.GetClosestCell(K3Pos);
		glm::dvec2 K3(0.0);
		if (K3CellIndex >= 0)
		{
			diff = inMACGrid.GetCellCenter(K3CellIndex) - K3Pos;

			weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

			if (mParticles[particleIndex].HasValidInertiaAnalog())
			{
				K3 = InterpolateVelocityFromGridCellBSpline(inMACGrid, K3Pos, K3CellIndex) + (mParticles[particleIndex].GetAffineState() * glm::inverse(mParticles[particleIndex].GetInertiaAnalog())) * diff;
			}
			else
			{
				K3 = InterpolateVelocityFromGridCellBSpline(inMACGrid, K3Pos, K3CellIndex);
			}
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

void APICFluid2D::InterpolateToGrid(MACGrid2D& inMACGrid)
{
	std::vector<double> contributedVelWeights;
	std::vector<double> contributedMassWeights;

	std::vector<double> interpXVelocities;
	std::vector<double> interpYVelocities;
	std::vector<double> interpMass;

	int numCells = inMACGrid.GetNumCells();

	interpXVelocities.assign(numCells, 0.0);
	interpYVelocities.assign(numCells, 0.0);

	interpMass.assign(numCells, 0.0);

	contributedVelWeights.assign(numCells, 0.0);
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

				if (inMACGrid.GetCellType(nearbyCellIndex) == CellType::eSOLID)
				{
					continue;
				}

				glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();
				glm::dvec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);

				glm::dvec2 diff = nearbyCellPos - particlePos;
				// TO DO: are these weights correct? Does the diff to each velocity face need to be considered?
				double weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

				int a = i - (x - 1);
				int b = j - (y - 1);

				weights[a][b] = weight;
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

				if (inMACGrid.GetCellType(nearbyCellIndex) == CellType::eSOLID)
				{
					continue;
				}

				glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();
				glm::dvec2 nearbyCellPos = inMACGrid.GetCellCenter(nearbyCellIndex);

				glm::dvec2 diff = nearbyCellPos - particlePos;

				double weight = InterpolateSupport(diff, inMACGrid.GetInverseCellSize());

				// Calculate contribution to cell mass
				double particleMass = mParticles[particleIndex].GetMass();
				interpMass[nearbyCellIndex] += particleMass * weight;

				// Calculate contribution to cell velocity.
				glm::dvec2 particleVel = mParticles[particleIndex].GetVelocity();

				glm::dvec2 velocityContribution(0.0);

				// Calculate separate weights for each face velocity.
				glm::dvec2 diffToX = diff - glm::dvec2(0.5 * inMACGrid.GetCellSize(), 0.0);
				double weightX = InterpolateSupport(diffToX, inMACGrid.GetInverseCellSize());

				glm::dvec2 diffToY = diff - glm::dvec2(0.0, 0.5 * inMACGrid.GetCellSize());
				double weightY = InterpolateSupport(diffToY, inMACGrid.GetInverseCellSize());

				// Include angular velocity if there is a valid inertia
				if (mParticles[particleIndex].HasValidInertiaAnalog())
				{
					glm::dmat2 angularMomentum = mParticles[particleIndex].GetAffineState() * glm::inverse(mParticles[particleIndex].GetInertiaAnalog());

					// Calculate velocity contribution for each face.
					velocityContribution += weightX * particleMass * (particleVel + angularMomentum * diffToX);
					velocityContribution += weightY * particleMass * (particleVel + angularMomentum * diffToY);
				}
				else
				{
					int a = i - (x - 1);
					int b = j - (y - 1);
					double deltaWeightX = weights[a + 1][b] - weights[a][b];
					double deltaWeightY = weights[a][b + 1] - weights[a][b];

					glm::dvec2 weightGradient(deltaWeightX, deltaWeightY);

					velocityContribution = particleMass * (particleVel + mParticles[particleIndex].GetAffineState() * weightGradient);
				}

				interpXVelocities[nearbyCellIndex] += velocityContribution.x;
				interpYVelocities[nearbyCellIndex] += velocityContribution.y;

				// Update contributed weights.
				contributedVelWeights[nearbyCellIndex] += weightX + weightY;
				contributedMassWeights[nearbyCellIndex] += weight;
			}
		}
	}

	for (int c = 0; c < numCells; c++)
	{
		// Normalize all values by dividing by weight.
		if (contributedVelWeights[c] != 0.0)
		{
			interpXVelocities[c] *= 1.0 / contributedVelWeights[c];
			interpYVelocities[c] *= 1.0 / contributedVelWeights[c];
		}
		if (contributedMassWeights[c] != 0.0)
		{
			interpMass[c] *= 1.0 / contributedMassWeights[c];
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
		mParticles[particleIndex].SetVelocity(newVelocity);

		InterpolateAffineFromGridCellBSpline(inMACGrid, particleIndex, cellIndex);
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


void APICFluid2D::InterpolateAffineFromGridCellBSpline(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
{
	glm::dvec2 particlePos = mParticles[particleIndex].GetPosition();
	double particleMass = mParticles[particleIndex].GetMass();

	int x, y;
	std::tie(x, y) = inMACGrid.GetXYFromIndex(cellIndex);

	glm::dmat2 affineSum(0.0);
	glm::dmat2 inertiaSum(0.0);
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

			glm::dvec2 diff = nearbyCellPos - particlePos;

			// Handle the X velocity.
			glm::dvec2 diffToX = diff - glm::dvec2(0.5 * inMACGrid.GetCellSize(), 0.0);
			glm::dvec2 velocityX(inMACGrid.GetCellXVelocity(nearbyCellIndex), 0.0);

			double weightX = InterpolateSupport(diffToX, inMACGrid.GetInverseCellSize());

			// Calculate contributions to inertia and affine state of cell.
			glm::dmat2 inertiaAnalogX = glm::outerProduct(diffToX, diffToX) * weightX;
			glm::dmat2 affineStateX = glm::outerProduct(velocityX, diffToX) * weightX;

			// Handle the Y velocity.
			glm::dvec2 diffToY = diff - glm::dvec2(0.0, 0.5 * inMACGrid.GetCellSize());
			glm::dvec2 velocityY(0.0, inMACGrid.GetCellYVelocity(nearbyCellIndex));

			double weightY = InterpolateSupport(diffToY, inMACGrid.GetInverseCellSize());

			// Calculate contributions to inertia and affine state of cell.
			glm::dmat2 inertiaAnalogY = glm::outerProduct(diffToY, diffToY) * weightY;
			glm::dmat2 affineStateY = glm::outerProduct(velocityY, diffToY) * weightY;

			// Add to sum.
			inertiaSum += inertiaAnalogX;
			inertiaSum += inertiaAnalogY;

			affineSum += affineStateX;
			affineSum += affineStateY;

			totalWeight += weightX;
			totalWeight += weightY;
		}
	}

	// Normalize the sums incase weight doesn't add to 1.
	if (totalWeight > 0.0)
	{
		inertiaSum /= totalWeight;
		affineSum /= totalWeight;
	}  
	else
	{
		inertiaSum = glm::dmat2(0.0);
		affineSum = glm::dmat2(0.0);
	}
	

	// Set particle values.
	mParticles[particleIndex].SetAffineState(affineSum);

	// If the inertia sum has an index, calculate and return the velocity derivatives matrix.
	if (glm::determinant(inertiaSum) != 0.0)
	{
		mParticles[particleIndex].SetInertiaAnalog(inertiaSum);
		mParticles[particleIndex].SetValidInertiaAnalog(true);
	}
	else
	{
		mParticles[particleIndex].SetInertiaAnalog(glm::dmat3(0.0));
		mParticles[particleIndex].SetValidInertiaAnalog(false);
	}
}