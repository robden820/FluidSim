#include "MACGrid3D.h"

#include <iostream>
#include "GLFW/glfw3.h"

#define TOLERANCE 0.0000000001f
#define PRECON_TUNER 0.97f
#define PRECON_SAFETY 0.25f

using namespace oneapi;

MACGrid3D::MACGrid3D(float inLeft, float inBottom, float inBack, float inWidth, float inHeight, float inLength, const std::vector<glm::vec3>& inParticlePositions, int inGridResolution)
{
	InitializeFromDimensions(inLeft, inBottom, inBack, inWidth, inHeight, inLength, inGridResolution);
	InitializeCellsFromParticles(inParticlePositions);
}

void MACGrid3D::InitializeFromDimensions(float inLeft, float inBottom, float inBack, float inWidth, float inHeight, float inLength, int inGridResolution)
{
	dLeft = inLeft;
	dBack = inBack;
	dBottom = inBottom;

	float dLength = inLength;
	float dWidth = inWidth;
	float dHeight = inHeight;

	mNumCellLength = inGridResolution;
	mNumCellWidth = inGridResolution;
	mNumCellHeight = inGridResolution;

	mNumCells = mNumCellLength * mNumCellWidth * mNumCellHeight;

	mCellSize = dLength / inGridResolution;
	mInvCellSize = 1 / mCellSize;

	float halfCell = mCellSize * 0.5f;

	mCellCenters.reserve(mNumCells);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		float centerX = dLeft + (x * mCellSize) + halfCell;
		float centerY = dBottom + (y * mCellSize) + halfCell;
		float centerZ = dBack + (z * mCellSize) + halfCell;

		mCellCenters.push_back(glm::vec3(centerX, centerY, centerZ));
	}

	mCellPressures.assign(mNumCells, 0.f);
	mCellXVelocities.assign(mNumCells, 0.f);
	mCellYVelocities.assign(mNumCells, 0.f);
	mCellZVelocities.assign(mNumCells, 0.f);
	mCellDivergence.assign(mNumCells, 0.f);
	mCellType.assign(mNumCells, CellType::eNONE);

	mIntXVelocities.assign(mNumCells, 0.f);
	mIntYVelocities.assign(mNumCells, 0.f);
	mIntZVelocities.assign(mNumCells, 0.f);

	mDensity = 1000.0f;
}

void MACGrid3D::InitializeCellsFromParticles(const std::vector<glm::vec3>& inParticlePositions)
{
	// Set edge cells to be solid.
	tbb::parallel_for(0, mNumCells, 1, [&](int cIndex)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(cIndex);

		if (x == 0 || x == mNumCellWidth - 1 || y == 0 || y == mNumCellHeight - 1 || z == 0 || z == mNumCellLength - 1)
		{
			mCellType[cIndex] = CellType::eSOLID;

		}
	});

	// For each particle, find the closest cell and mark it as a fluid.
	tbb::parallel_for(0, (int)inParticlePositions.size(), 1, [&](int pIndex)
	{
		int cellIndex = GetClosestCell(inParticlePositions[pIndex]);

		mCellType[cellIndex] = CellType::eFLUID;
		mCellPressures[cellIndex] = 1.0f;

	});
}

void MACGrid3D::Update(float deltaTime)
{
	//Advection
	float start = glfwGetTime();
	AdvectCellVelocity(deltaTime);

	std::cout << "MAC: advect: " << glfwGetTime() - start << "\n";

	// Projection step
	start = glfwGetTime();
	
	UpdateCellPressure(deltaTime, 50);

	std::cout << "MAC: Pressure: " << glfwGetTime() - start << "\n";

	// Velocity update
	start = glfwGetTime();
	UpdateCellVelocity(deltaTime);

	std::cout << "MAC: cell update: " << glfwGetTime() - start << "\n";
}

void MACGrid3D::AdvectCellVelocity(float deltaTime)
{
	tbb::parallel_for(0, mNumCells, 1, [&](int index)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		float xVelocity = mCellXVelocities[index];
		float yVelocity = mCellYVelocities[index];
		float zVelocity = mCellZVelocities[index];

		if (x > 0)
		{
			int neighbourLeft = GetIndexFromXYZ(x - 1, y, z);

			xVelocity += mCellXVelocities[neighbourLeft];
			xVelocity *= 0.5f;
		}

		if (y > 0)
		{
			int neighbourBottom = GetIndexFromXYZ(x, y - 1, z);

			yVelocity += mCellYVelocities[neighbourBottom];
			yVelocity *= 0.5f;
		}

		if (z > 0)
		{
			int neighbourBack = GetIndexFromXYZ(x, y, z - 1);

			zVelocity += mCellZVelocities[neighbourBack];
			zVelocity *= 0.5f;
		}

		glm::vec3 avgVelocity(xVelocity, yVelocity, zVelocity);

		// Want to find previous position, so go backwards in time.
		glm::vec3 prevPosition = mCellCenters[index] - avgVelocity * deltaTime;

		int prevCellIndex = GetClosestCell(prevPosition);

		mIntXVelocities[index] = mCellXVelocities[prevCellIndex];
		mIntYVelocities[index] = mCellYVelocities[prevCellIndex];
		mIntZVelocities[index] = mCellZVelocities[prevCellIndex];

		int xPrev = x, yPrev = y, zPrev = z;
		std::tie(xPrev, yPrev, zPrev) = GetXYZFromIndex(index);
		if (x > 0)
		{
			int prevNeighbourLeft = GetIndexFromXYZ(x - 1, y, z);
			mIntXVelocities[index] += mCellXVelocities[prevNeighbourLeft];
			mIntXVelocities[index] *= 0.5f;
		}
		if (y > 0)
		{
			int prevNeighbourBottom = GetIndexFromXYZ(x, y - 1, z);
			mIntYVelocities[index] += mCellYVelocities[prevNeighbourBottom];
			mIntYVelocities[index] *= 0.5f;
		}
		
		if (z > 0)
		{
			int prevNeighbourBack = GetIndexFromXYZ(x, y, z - 1);
			mIntZVelocities[index] += mCellZVelocities[prevNeighbourBack];
			mIntZVelocities[index] *= 0.5f;
		}
	});
}

int MACGrid3D::GetClosestCell(const glm::vec3& inPosition)
{
	int x = floor((inPosition.x - dLeft - (mCellSize * 0.5f)) * mInvCellSize);
	int y = floor((inPosition.y - dBottom - (mCellSize * 0.5f)) * mInvCellSize);
	int z = floor((inPosition.z - dBack - (mCellSize * 0.5f)) * mInvCellSize);

	int approxIndex = GetIndexFromXYZ(x, y, z);

	glm::vec3 vec = inPosition - mCellCenters[approxIndex];

	float distSqr = (vec.x * vec.x) + (vec.y * vec.y) + (vec.z * vec.z);

	float closestDistSqr = distSqr;
	int closestIndex = approxIndex;

	if (x > 0)
	{
		int neighbourLeft = GetIndexFromXYZ(x - 1, y, z);

		glm::vec3 cToP = inPosition - mCellCenters[neighbourLeft];

		float dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourLeft;
		}
	}
	if (x < mNumCellWidth - 1)
	{
		int neighbourRight = GetIndexFromXYZ(x + 1, y, z);

		glm::vec3 cToP = inPosition - mCellCenters[neighbourRight];

		float dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourRight;
		}
	}

	if (y > 0)
	{
		int neighbourBottom = GetIndexFromXYZ(x, y - 1, z);

		glm::vec3 cToP = inPosition - mCellCenters[neighbourBottom];

		float dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourBottom;
		}
	}
	if (y < mNumCellHeight - 1)
	{
		int neighbourTop = GetIndexFromXYZ(x, y + 1, z);

		glm::vec3 cToP = inPosition - mCellCenters[neighbourTop];

		float dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourTop;
		}
	}

	if (z > 0)
	{
		int neighbourBack = GetIndexFromXYZ(x, y, z - 1);

		glm::vec3 cToP = inPosition - mCellCenters[neighbourBack];

		float dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourBack;
		}
	}
	if (z < mNumCellWidth - 1)
	{
		int neighbourFront = GetIndexFromXYZ(x, y, z + 1);

		glm::vec3 cToP = inPosition - mCellCenters[neighbourFront];

		float dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourFront;
		}
	}

	return closestIndex;
}

void MACGrid3D::UpdateCellVelocity(float deltaTime)
{
	float scale = deltaTime * mInvCellSize * (1.0f / mDensity);

	float solidXVel = 0.0f;
	float solidYVel = 0.0f;
	float solidZVel = 0.0f;

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (x > 0)
		{
			int neighbourLeft = GetIndexFromXYZ(x - 1, y, z);

			if (mCellType[neighbourLeft] == CellType::eFLUID || mCellType[index] == CellType::eFLUID)
			{
				if (mCellType[neighbourLeft] == CellType::eSOLID || mCellType[index] == CellType::eSOLID)
				{
					mCellXVelocities[index] = solidXVel;
				}
				else
				{
					mCellXVelocities[index] -= scale * (mCellPressures[index] - mCellPressures[neighbourLeft]);
				}
			}
			//else mark as unknown
		}
		if (y > 0)
		{
			int neighbourBottom = GetIndexFromXYZ(x, y - 1, z);

			if (mCellType[neighbourBottom] == CellType::eFLUID || mCellType[index] == CellType::eFLUID)
			{
				if (mCellType[neighbourBottom] == CellType::eSOLID || mCellType[index] == CellType::eSOLID)
				{
					mCellYVelocities[index] = solidYVel;
				}
				else
				{
					mCellYVelocities[index] -= scale * (mCellPressures[index] - mCellPressures[neighbourBottom]);
				}
			}
			//else mark as unknown
		}
		if (z > 0)
		{
			int neighbourBack = GetIndexFromXYZ(x, y, z - 1);

			if (mCellType[neighbourBack] == CellType::eFLUID || mCellType[index] == CellType::eFLUID)
			{
				if (mCellType[neighbourBack] == CellType::eSOLID || mCellType[index] == CellType::eSOLID)
				{
					mCellZVelocities[index] = solidYVel;
				}
				else
				{
					mCellZVelocities[index] -= scale * (mCellPressures[index] - mCellPressures[neighbourBack]);
				}
			}
			//else mark as unknown
		}	
	}
}

void MACGrid3D::CalculateCellDivergence(float deltaTime)
{
	float scale = mInvCellSize;

	float solidXVel = 0.0f;
	float solidYVel = 0.0f;
	float solidZVel = 0.0f;

	mCellDivergence.assign(mNumCells, 0.f);

	// Calculate negative divergence of each cell.
	for (int index = 0; index < mNumCells; index++)
	{
		if (mCellType[index] == CellType::eFLUID)
		{
			int x, y, z;
			std::tie(x, y, z) = GetXYZFromIndex(index);

			float divergence = 0.f;

			int neighbourRight = GetIndexFromXYZ(x + 1, y, z);
			divergence += mIntXVelocities[neighbourRight] - mIntXVelocities[index];

			int neighbourTop = GetIndexFromXYZ(x, y + 1, z);
			divergence += mIntYVelocities[neighbourTop] - mIntYVelocities[index];

			int neighbourFront = GetIndexFromXYZ(x, y, z + 1);
			divergence += mIntZVelocities[neighbourFront] - mIntZVelocities[index];

			divergence *= -scale;

			mCellDivergence[index] = divergence;
		}
	}

	// Update divergence to account for solid cells.
	for (int index = 0; index < mNumCells; index++)
	{
		if (mCellType[index] == CellType::eFLUID)
		{
			int x, y, z;
			std::tie(x, y, z) = GetXYZFromIndex(index);

			int neighbourLeft = GetIndexFromXYZ(x - 1, y, z);

			if (mCellType[neighbourLeft] == CellType::eSOLID)
			{
				mCellDivergence[index] -= scale * (mIntXVelocities[index] - solidXVel);
			}

			int neighbourRight = GetIndexFromXYZ(x + 1, y, z);

			if(mCellType[neighbourRight] == CellType::eSOLID)
			{
				mCellDivergence[index] += scale * (mIntXVelocities[neighbourRight] - solidXVel);
			}

			int neighbourBottom = GetIndexFromXYZ(x, y - 1, z);

			if (mCellType[neighbourBottom] == CellType::eSOLID)
			{
				mCellDivergence[index] -= scale * (mIntYVelocities[index] - solidYVel);
			}

			int neighbourTop = GetIndexFromXYZ(x, y + 1, z);

			if (mCellType[neighbourTop] == CellType::eSOLID)
			{
				mCellDivergence[index] += scale * (mIntYVelocities[neighbourTop] - solidYVel);
			}

			int neighbourBack = GetIndexFromXYZ(x, y, z - 1);

			if (mCellType[neighbourBack] == CellType::eSOLID)
			{
				mCellDivergence[index] -= scale * (mIntZVelocities[index] - solidZVel);
			}

			int neighbourFront = GetIndexFromXYZ(x, y, z + 1);

			if (mCellType[neighbourFront] == CellType::eSOLID)
			{
				mCellDivergence[index] += scale * (mIntZVelocities[neighbourFront] - solidZVel);
			}
		}
	}
}

void MACGrid3D::UpdateCellPressure(float deltaTime, int maxIterations)
{
	std::vector<float> Adiagonal;
	std::vector<float> Ax;
	std::vector<float> Ay;
	std::vector<float> Az;

	Adiagonal.assign(mNumCells, 0.0f);
	Ax.assign(mNumCells, 0.0f);
	Ay.assign(mNumCells, 0.0f);
	Az.assign(mNumCells, 0.0f);

	std::cout << "--------------------------------\n";
	float start = glfwGetTime();
	std::cout << "Initialize linear system: ";

	InitializeLinearSystem(deltaTime, Adiagonal, Ax, Ay, Az);

	std::cout << glfwGetTime() - start << "\n";

	// Preconditioned Conjugate Gradient.
	std::vector<float> newPressure;
	newPressure.assign(mNumCells, 0.0f);

	start = glfwGetTime();
	std::cout << "Calculate cell divergence: ";

	CalculateCellDivergence(deltaTime);

	std::cout << glfwGetTime() - start << "\n";
	
	std::vector<float> residuals;
	residuals = mCellDivergence;

	float maxResidual = -1.0f;

	for (int index = 0; index < mNumCells; index++)
	{
		if (abs(residuals[index]) > maxResidual)
		{
			maxResidual = residuals[index];
		}
	}

	if (maxResidual < TOLERANCE)
	{
		// Reached convergence;
		return;
	}


	std::vector<float> z;
	z.assign(mNumCells, 0.0f);

	start = glfwGetTime();
	std::cout << "Calculate preconditioner: ";

	std::vector<float> precon;
	CalculatePreconditioner(precon, Adiagonal, Ax, Ay, Az);

	std::cout << glfwGetTime() - start << "\n";
	start = glfwGetTime();
	std::cout << "Apply preconditioner: ";

	ApplyPreconditioner(z, residuals, precon, Ax, Ay, Az);

	std::cout << glfwGetTime() - start << "\n";

	std::vector<float> search = z;

	float theta = 0.0f;
	for (int i = 0; i < mNumCells; i++)
	{
		theta += z[i] * residuals[i];
	}

	int iteration;

	for (iteration = 1; iteration <= maxIterations; ++iteration)
	{
		ApplyA(deltaTime, z, search, Adiagonal, Ax, Ay, Az);

		float phi = 0.0f;
		for (int i = 0; i < mNumCells; i++)
		{
			phi += z[i] * search[i];
		}

		if (abs(phi) < TOLERANCE)
		{
			phi = TOLERANCE;
		}

		float alpha = theta / phi;

		for (int i = 0; i < mNumCells; i++)
		{
			newPressure[i] += alpha * search[i];
			float diff = alpha * z[i];
			residuals[i] -= diff;
		}

		maxResidual = -1.0f;

		for (int index = 0; index < mNumCells; index++)
		{
			if (abs(residuals[index]) > maxResidual)
			{
				maxResidual = abs(residuals[index]);
			}
		}

		if (maxResidual < TOLERANCE)
		{
			// Reached convergence;
			break;
		}

		ApplyPreconditioner(z, residuals, precon, Ax, Ay, Az);

		float thetaNew = 0.0f;
		for (int i = 0; i < mNumCells; i++)
		{
			thetaNew += z[i] * residuals[i];
			if (thetaNew > 0.f)
			{
				std::cout << "|" << thetaNew;
			}
		}

		if (abs(theta) < TOLERANCE)
		{
			theta = TOLERANCE;
		}
		float beta = thetaNew / theta;

		for (int i = 0; i < mNumCells; i++)
		{
			search[i] = z[i] + beta * search[i];
		}

		theta = thetaNew;
		std::cout << ".";
	}
	std::cout << "\n";
	if (iteration >= maxIterations)
	{
		std::cout << "WARNING: MAX NUMBER OF ITERATIONS REACHED IN PRESSURE SOLVE" << "\n";
		std::cout << "Check pressure solver for potential issues, MACGrid3D.cpp" << "\n";
	}
	std::cout << "NUM PRESSURE SOLVE ITERATIONS: " << iteration << "\n";
	std::cout << "------------------------------\n";

	// Set pressure
	for (int i = 0; i < mNumCells; i++)
	{
		mCellPressures[i] = newPressure[i];
	}
}

void MACGrid3D::InitializeLinearSystem(float deltaTime, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ)
{
	float scale = -deltaTime * mInvCellSize * mInvCellSize * (1.0f / mDensity);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			if (x > 0)
			{
				int neighbourLeft = GetIndexFromXYZ(x - 1, y, z);

				if (mCellType[neighbourLeft] == CellType::eFLUID)
				{
					inDiag[index] += scale;
				}
			}

			if (x < mNumCellWidth - 1)
			{
				int neighbourRight = GetIndexFromXYZ(x + 1, y, z);

				if (mCellType[neighbourRight] == CellType::eFLUID)
				{
					inDiag[index] += scale;
					inX[index] -= scale;
				}
				else if (mCellType[neighbourRight] != CellType::eSOLID)
				{
					inDiag[index] += scale;
				}
			}

			if (y > 0)
			{
				int neighbourBottom = GetIndexFromXYZ(x, y - 1, z);

				if (mCellType[neighbourBottom] == CellType::eFLUID)
				{
					inDiag[index] += scale;
				}
			}

			if (y < mNumCellHeight - 1)
			{
				int neighbourTop = GetIndexFromXYZ(x, y + 1, z);

				if (mCellType[neighbourTop] == CellType::eFLUID)
				{
					inDiag[index] += scale;
					inY[index] -= scale;
				}
				else if (mCellType[neighbourTop] != CellType::eSOLID)
				{
					inDiag[index] += scale;
				}
			}

			if (z > 0)
			{
				int neighbourBack = GetIndexFromXYZ(x, y, z - 1);

				if (mCellType[neighbourBack] == CellType::eFLUID)
				{
					inDiag[index] += scale;
				}
			}

			if (z < mNumCellWidth - 1)
			{
				int neighbourFront = GetIndexFromXYZ(x, y, z + 1);

				if (mCellType[neighbourFront] == CellType::eFLUID)
				{
					inDiag[index] += scale;
					inZ[index] -= scale;
				}
				else if (mCellType[neighbourFront] != CellType::eSOLID)
				{
					inDiag[index] += scale;
				}
			}
		}
	}
}

void MACGrid3D::ApplyA(float deltaTime, std::vector<float>& outResult, const std::vector<float>& inVec, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ)
{
	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		float value = 0.f;

		if (x > 0)
		{
			int neighbourLeft = GetIndexFromXYZ(x - 1, y, z);
			value += inVec[neighbourLeft];
		}

		if (x < mNumCellWidth - 1)
		{
			int neighbourRight = GetIndexFromXYZ(x + 1, y, z);
			value += inVec[neighbourRight];
		}

		if (y > 0)
		{
			int neighbourBottom = GetIndexFromXYZ(x, y - 1, z);
			value += inVec[neighbourBottom];
		}

		if (y < mNumCellHeight - 1)
		{
			int neighbourTop = GetIndexFromXYZ(x, y + 1, z);
			value += inVec[neighbourTop];
		}

		if (z > 0)
		{
			int neighbourBack = GetIndexFromXYZ(x, y, z - 1);
			value += inVec[neighbourBack];
		}

		if (z < mNumCellLength - 1)
		{
			int neighbourFront = GetIndexFromXYZ(x, y, z + 1);
			value += inVec[neighbourFront];
		}

		value += inDiag[index] * inVec[index];
		outResult[index] = value;
	}
}

void MACGrid3D::CalculatePreconditioner(std::vector<float>& inOutPrecon, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ)
{
	inOutPrecon.assign(mNumCells, 0.f);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			float Axi = 0.0f; // Aplusi_iminus1
			float Axj = 0.0f; // Aplusi_jminus1
			float Axk = 0.0f; // Aplusi_kminus1

			float Ayi = 0.0f; // Aplusj_iminus1
			float Ayj = 0.0f; // Aplusj_jminus1
			float Ayk = 0.0f; // Aplusj_kminus1

			float Azi = 0.0f; // Aplusk_iminus1
			float Azj = 0.0f; // Aplusk_jminus1
			float Azk = 0.0f; // Aplusk_kminus1

			float preconi = 0.0f; // precon_iminus1
			float preconj = 0.0f; // precon_jminus1
			float preconk = 0.0f; // precon_kminus1

			if (x > 0)
			{
				int neighbourLeft = GetIndexFromXYZ(x - 1, y, z);;

				Axi = inX[neighbourLeft];
				Ayi = inY[neighbourLeft];
				Azi = inZ[neighbourLeft];

				preconi = inOutPrecon[neighbourLeft];
			}
					
			if (y > 0)
			{
				int neighbourBottom = GetIndexFromXYZ(x, y - 1, z);;

				Axj = inX[neighbourBottom];
				Ayj = inY[neighbourBottom];
				Azj = inZ[neighbourBottom];

				preconj = inOutPrecon[neighbourBottom];
			}

			if (z > 0)
			{
				int neighbourBack = GetIndexFromXYZ(x, y, z - 1);

				Axk = inX[neighbourBack];
				Ayk = inY[neighbourBack];
				Azk = inZ[neighbourBack];

				preconk = inOutPrecon[neighbourBack];
			}

			float a = Axi * preconi;
			float b = Ayj * preconj;
			float c = Azk * preconk;

			float termOne = a * a + b * b + c * c;

			float d = Axi * (Ayi + Azi) * preconi * preconi;
			float e = Ayj * (Axj + Azj) * preconj * preconj;
			float f = Azk * (Axk + Ayk) * preconk * preconk;

			float termTwo = d + e + f;

			float newPrecon = inDiag[index] - termOne - PRECON_TUNER * termTwo;

			if (newPrecon < PRECON_SAFETY * inDiag[index])
			{
				newPrecon = inDiag[index];
			}
					
			if (newPrecon > TOLERANCE)
			{
				inOutPrecon[index] = 1.0f / sqrt(newPrecon);
			}	
		}
	}
}

void MACGrid3D::ApplyPreconditioner(std::vector<float>& outResult, const std::vector<float>& inResidual, const std::vector<float>& inPrecon, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ)
{
	std::vector<float> intermediate;  // q
	intermediate.assign(mNumCells, 0.f);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{

			float Axi = 0.0f; // Aplusi_iminus1
			float Ayj = 0.0f; // Aplusj_jminus1
			float Azk = 0.0f; // Aplusk_kminus1

			float preconi = 0.0f; // precon_iminus1
			float preconj = 0.0f; // precon_jminus1
			float preconk = 0.0f; // precon_kminus1

			float intermediatei = 0.0f; // qminusi
			float intermediatej = 0.0f; // qminusj
			float intermediatek = 0.0f; // qminusk

			int neighbourLeft = GetIndexFromXYZ(x - 1, y, z);

			Axi = inX[neighbourLeft];
			preconi = inPrecon[neighbourLeft];
			intermediatei = intermediate[neighbourLeft];

			int neighbourBottom = GetIndexFromXYZ(x, y - 1, z);

			Ayj = inY[neighbourBottom];
			preconj = inPrecon[neighbourBottom];
			intermediatej = intermediate[neighbourBottom];

			int neighbourBack = GetIndexFromXYZ(x, y, z - 1);

			Azk = inZ[neighbourBack];
			preconk = inPrecon[neighbourBack];
			intermediatek = intermediate[neighbourBack];

			float t = inResidual[index] - Axi * preconi * intermediatei
										- Ayj * preconj * intermediatej
										- Azk * preconk * intermediatek;

			intermediate[index] = t * inPrecon[index];
		}
	}

	for (int index = mNumCells - 1; index >= 0; index--)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			float zi = 0.0f;
			float zj = 0.0f;
			float zk = 0.0f;

			int neighbourRight = GetIndexFromXYZ(x + 1, y, z);

			zi = outResult[neighbourRight];

			int neighbourTop = GetIndexFromXYZ(x, y + 1, z);

			zj = outResult[neighbourTop];

			int neighbourFront= GetIndexFromXYZ(x, y, z + 1);

			zk = outResult[neighbourFront];

			float t = intermediate[index] - inX[index] * inPrecon[index] * zi
										  - inY[index] * inPrecon[index] * zj
										  - inZ[index] * inPrecon[index] * zk;

			outResult[index] = t * inPrecon[index];
		}
	}
}

const CellType MACGrid3D::GetCellTypeFromPosition(const glm::vec3& inPosition)
{
	int index = GetClosestCell(inPosition);

	return GetCellType(index);
}

std::tuple<int, int, int> MACGrid3D::GetXYZFromIndex(int index)
{
	int x = 0;
	int y = 0;
	int z = 0;

	z = index % mNumCellLength;
	y = static_cast<int>(round((index - z) / mNumCellLength)) % mNumCellHeight;
	x = static_cast<int>(round((index - y * mNumCellHeight - z) / (mNumCellLength * mNumCellHeight))) % mNumCellWidth;

	return std::tuple<int, int, int>(x, y, z);
}

int MACGrid3D::GetIndexFromXYZ(int X, int Y, int Z)
{
	return  Z + Y * mNumCellLength + X * mNumCellHeight * mNumCellLength;
}