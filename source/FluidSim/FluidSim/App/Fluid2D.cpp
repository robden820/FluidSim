#include "Fluid2D.h"
#include <iostream>

#include "oneapi/tbb.h"
#include "GLFW/glfw3.h"

Fluid2D::Fluid2D(const ApplicationData& inData)
{
	float start = glfwGetTime();
	std::cout << "Initializing particles: ";
	// Initialize particles
	mParticles.reserve(inData.GetNumParticles());
	mParticlePositions.reserve(inData.GetNumParticles());

	for (int x = 0; x < 10; x++)
	{
		for (int y = 0; y < 10; y++)
		{
				glm::vec2 position((x - 5) * 0.25f, (y - 5) * 0.25f);

				Particle2D particle(position);

				mParticles.push_back(particle);
				mParticlePositions.push_back(position);
		}
	}

	std::cout << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	std::cout << "------------------------ \n";
	std::cout << "Initializing MAC Grid: \n";

	// Initialize Grid.
	MACGrid2D grid(inData);
	mMACGrid = grid;

	std::cout << "Total : " << glfwGetTime() - start << "\n";
	std::cout << "------------------------ \n";
}

void Fluid2D::Update(ApplicationData& inOutData)
{
	float deltaTime = inOutData.GetDeltaTime();
	// Transfer particle velocities to grid.
	InterpolateToGrid();

	// Update grid velocities.
	mMACGrid.Update(inOutData);

	// Interpolate velocities back to particles.
	InterpolateFromGrid();

	// Advect particles.
	for (int p = 0; p < GetNumParticles(); p++)
	{
		mParticles[p].StepParticle(deltaTime);
		mParticlePositions[p] = mParticles[p].GetPosition();
	}

	/* TO DO: this should happen automatically with correct handling of solid pressures.
	// Ensure particles stay inside the simulation domain.
	for (int p = 0; p < GetNumParticles(); p++)
	{
		if (!mDomain.IsPointInDomain(mParticles[p].GetPosition()))
		{
			ClampParticleToDomain(mParticles[p]);
		}
	}
	*/

	inOutData.Set2DParticlePositions(mParticlePositions);
}

/* TO DO: this should happen automatically with correct handling of solid pressures.
void Fluid2D::ClampParticleToDomain(Particle2D& particle)
{
	glm::vec2 particlePos = particle.GetPosition();
	glm::vec2 particleVel = particle.GetVelocity();

	if (particlePos.x < mDomain.GetLeft())
	{
		particlePos.x = mDomain.GetLeft();
		particleVel.x = 0.0f;
	}
	else if (particlePos.x > mDomain.GetRight())
	{
		particlePos.x = mDomain.GetRight();
		particleVel.x = 0.0f;
	}

	if (particlePos.y < mDomain.GetBottom())
	{
		particlePos.y = mDomain.GetBottom();
		particleVel.y = 0.0f;
	}
	else if (particlePos.y > mDomain.GetTop())
	{
		particlePos.y = mDomain.GetTop();
		particleVel.y = 0.0f;
	}

	particle.SetPosition(particlePos);
	particle.SetVelocity(particleVel);
}
*/

void Fluid2D::InterpolateToGrid()
{
	for (int c = 0; c < mMACGrid.GetNumCells(); c++)
	{
		if (mMACGrid.GetCellType(c) != CellType::eSOLID)
		{
			mMACGrid.SetCellType(c, CellType::eAIR);
		}
	}

	std::vector<float> contributedXWeights;
	std::vector<float> contributedYWeights;

	std::vector<float> interpXVelocities;
	std::vector<float> interpYVelocities;

	interpXVelocities.assign(mMACGrid.GetNumCells(), 0.f);
	interpYVelocities.assign(mMACGrid.GetNumCells(), 0.f);

	contributedXWeights.assign(mMACGrid.GetNumCells(), 0.f);
	contributedYWeights.assign(mMACGrid.GetNumCells(), 0.f);

	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(mParticles[p]);

		glm::vec2 particlePos = mParticles[p].GetPosition();
		glm::vec2 cellPos = mMACGrid.GetCellCenter(cellIndex);

		float xDiff = particlePos.x - cellPos.x;
		float yDiff = particlePos.y - cellPos.y;

		float xWeight = (xDiff / mMACGrid.GetCellSize()) + 0.5f;
		float yWeight = (yDiff / mMACGrid.GetCellSize()) + 0.5f;

		float velocityX = mMACGrid.GetCellXVelocity(cellIndex) * xWeight;
		float velocityY = mMACGrid.GetCellYVelocity(cellIndex) * yWeight;

		if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
		{
			mMACGrid.SetCellType(cellIndex, CellType::eFLUID);

			interpXVelocities[cellIndex] += mParticles[p].GetVelocity().x * xWeight;
			interpYVelocities[cellIndex] += mParticles[p].GetVelocity().y * yWeight;

			contributedXWeights[cellIndex] += xWeight;
			contributedYWeights[cellIndex] += yWeight;

			int x, y;
			std::tie(x, y) = mMACGrid.GetXYFromIndex(cellIndex);

			if (x < 0)
			{
				int neighbourLeft = mMACGrid.GetIndexFromXY(x - 1, y);

				if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
				{
					interpXVelocities[neighbourLeft] += mParticles[p].GetVelocity().x * (1 - xWeight);
					interpYVelocities[neighbourLeft] += mParticles[p].GetVelocity().y * (1 - yWeight);

					contributedXWeights[neighbourLeft] += (1 - xWeight);
					contributedYWeights[neighbourLeft] += (1 - yWeight);
				}
			}
			if (y < 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXY(x, y - 1);

				if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
				{
					interpXVelocities[neighbourBottom] += mParticles[p].GetVelocity().x * (1 - xWeight);
					interpYVelocities[neighbourBottom] += mParticles[p].GetVelocity().y * (1 - yWeight);

					contributedXWeights[neighbourBottom] += (1 - xWeight);
					contributedYWeights[neighbourBottom] += (1 - yWeight);
				}
			}
		}
	}

	for (int c = 0; c < mMACGrid.GetNumCells(); c++)
	{
		if (contributedXWeights[c] != 0.f)
		{
			interpXVelocities[c] *= 1.0f / contributedXWeights[c];
		}
		if (contributedYWeights[c] != 0.f)
		{
			interpYVelocities[c] *= 1.0f / contributedYWeights[c];
		}

		if (mMACGrid.GetCellType(c) != CellType::eSOLID)
		{
			mMACGrid.SetCellXVelocity(c, interpXVelocities[c]);
			mMACGrid.SetCellYVelocity(c, interpYVelocities[c]);
		}
	}
}

void Fluid2D::InterpolateFromGrid()
{
	for (int p = 0; p < GetNumParticles(); p++)
	{
		int cellIndex = ClosestCellToParticle(mParticles[p]);

		if (mMACGrid.GetCellType(cellIndex) != CellType::eSOLID)
		{
			int x, y;
			std::tie(x, y) = mMACGrid.GetXYFromIndex(cellIndex);

			glm::vec2 particlePos = mParticles[p].GetPosition();
			glm::vec2 cellPos = mMACGrid.GetCellCenter(cellIndex);

			float xDiff = particlePos.x - cellPos.x;
			float yDiff = particlePos.y - cellPos.y;

			float xWeight = (xDiff / mMACGrid.GetCellSize()) + 0.5f;
			float yWeight = (yDiff / mMACGrid.GetCellSize()) + 0.5f;

			float velocityX = mMACGrid.GetCellXVelocity(cellIndex) * xWeight;
			float velocityY = mMACGrid.GetCellYVelocity(cellIndex) * yWeight;

			if (x > 0)
			{
				int neighbourLeft = mMACGrid.GetIndexFromXY(x - 1, y);

				velocityX += mMACGrid.GetCellXVelocity(neighbourLeft) * (1 - xWeight);
			}
			if (y > 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXY(x, y - 1);

				velocityY += mMACGrid.GetCellYVelocity(neighbourBottom) * (1 - yWeight);
			}

			// Handle bottom solid boundary
			if (y > 0)
			{
				int neighbourBottom = mMACGrid.GetIndexFromXY(x, y - 1);

				if (mMACGrid.GetCellType(neighbourBottom) == CellType::eSOLID)
				{
					mParticles[p].ApplyForce(glm::vec2{ 0.f, 9.8f } * mParticles[p].GetMass());
				}
			}

			// Handle left and right solid boundaries
			if (x > 0)
			{
				int neighbourLeft = mMACGrid.GetIndexFromXY(x - 1, y);

				if (mMACGrid.GetCellType(neighbourLeft) == CellType::eSOLID)
				{
					velocityX = 0.f;
				}
			}
			if (x < mMACGrid.GetNumCellsWidth() - 1)
			{
				int neighbourRight = mMACGrid.GetIndexFromXY(x + 1, y);

				if (mMACGrid.GetCellType(neighbourRight) == CellType::eSOLID)
				{
					velocityX = 0.f;
				}
			}

			mParticles[p].SetVelocity(glm::vec2(velocityX, velocityY));
		}
		else
		{
			std::cout << "WARNING - particle may have penetrated solid boundary. \n";
			std::cout << "Fluid2D.pp - InterpolateFromGrid function. \n";
		}
	}
}

int Fluid2D::ClosestCellToParticle(const Particle2D& particle)
{
	glm::vec2 particlePos = particle.GetPosition();

	return mMACGrid.GetClosestCell(particlePos);
}