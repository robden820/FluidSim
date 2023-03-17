# include "Fluid2D.h"

template <typename Tparticle>
void Fluid2D<Tparticle>::SeedParticles(const ApplicationData& inOutData)
{
	double gridCellSize = inOutData.GetGridCellSize();
	double halfGridCell = gridCellSize * 0.5;

	// To do: reserve memory space for particles array.

	for (int cellIndex = 0; cellIndex < inOutData.GetNumGridCells(); cellIndex++)
	{
		if (inOutData.GetCellType(cellIndex) == CellType::eFLUID)
		{
			// Seed 4 particles at each fluid cell.
			// To do : allow a user to set this value.

			for (int i = 0; i < 4; i++)
			{
				glm::dvec2 particlePos = inOutData.GetCellCenter2D(cellIndex);

				particlePos.x += (rand() % 100) * 0.01 * gridCellSize + halfGridCell;
				particlePos.y += (rand() % 100) * 0.01 * gridCellSize + halfGridCell;

				Tparticle particle(particlePos);

				mParticles.push_back(particle);
				mParticlePositions.push_back(particlePos);
			}
		}
	}
}

template <typename Tparticle>
int Fluid2D<Tparticle>::ClosestCellToParticle(const MACGrid2D& inMACGrid, const Tparticle& particle)
{
	glm::dvec2 particlePos = particle.GetPosition();

	return inMACGrid.GetClosestCell(particlePos);
}

template <typename Tparticle>
void Fluid2D<Tparticle>::ProjectParticleToFluid(const MACGrid2D& inMACGrid, int particleIndex, int cellIndex)
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

	mParticles[particleIndex].SetVelocity(glm::dvec2(0.0, 0.0));
	mParticles[particleIndex].SetPosition(mParticlePositions[particleIndex]);
}

template <typename Tparticle>
void Fluid2D<Tparticle>::DeleteBadParticles(const MACGrid2D& inMACGrid)
{
	std::vector<Tparticle> safeParticles;
	std::vector<glm::dvec2> safeParticlePositions;
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

template <typename Tparticle>
double Fluid2D<Tparticle>::InterpolateSupport(const glm::dvec2& diff, double invCellSize)
{
	glm::dvec2 scaled = diff * invCellSize;

	return BSpline(scaled.x) * BSpline(scaled.y);
}

template <typename Tparticle>
double Fluid2D<Tparticle>::BSpline(double input)
{
	double output = 0.0;

	// If input is in [-1.5, 0.5)
	if (input >= -1.5 && input < -0.5)
	{
		double a = input + 1.5;
		output = 0.5 * a * a;
	}
	// if input is in [-0.5, 0.5)
	else if (input >= -0.5 && input < 0.5)
	{
		output = 0.75 - input * input;
	}
	// If input is in [0.5, 1.5)
	else if (input >= 0.5 && input < 1.5)
	{
		double a = 1.5 - input;
		output = 0.5 * a * a;
	}

	return output;
}