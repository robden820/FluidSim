#include "MACGrid3D.h"

#include <iostream>
#include "GLFW/glfw3.h"

#define TOLERANCE 0.0000000001f
#define PRECON_TUNER 0.97f
#define PRECON_SAFETY 0.25f

using namespace oneapi;

MACGrid3D::MACGrid3D(const ApplicationData& inData)
{
	InitializeGrid(inData);
	InitializeCellsFromParticles(inData.Get3DParticlePositions());
}

void MACGrid3D::InitializeGrid(const ApplicationData& inData)
{
	dLeft = inData.GetGridLeft();
	dBack = inData.GetGridBack();
	dBottom = inData.GetGridBottom();

	double dLength = inData.GetGridLength();
	double dWidth = inData.GetGridWidth();
	double dHeight = inData.GetGridHeight();

	mNumCellLength = inData.GetNumGridCellsLength();
	mNumCellWidth = inData.GetNumGridCellsWidth();
	mNumCellHeight = inData.GetNumGridCellsHeight();

	mNumCells = inData.GetNumGridCells();

	mCellSize = inData.GetGridCellSize();
	mInvCellSize = 1 / mCellSize;

	double halfCell = mCellSize * 0.5f;

	mCellCenters.reserve(mNumCells);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		double centerX = dLeft + (x * mCellSize) + halfCell;
		double centerY = dBottom + (y * mCellSize) + halfCell;
		double centerZ = dBack + (z * mCellSize) + halfCell;

		mCellCenters.push_back(glm::dvec3(centerX, centerY, centerZ));
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

	mDensity = inData.GetFluidDensity();
}

void MACGrid3D::InitializeCellsFromParticles(const std::vector<glm::dvec3>& inParticlePositions)
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
		mCellPressures[cellIndex] = 1.0;

	});
}

void MACGrid3D::Update(ApplicationData& inOutData)
{
	double deltaTime = inOutData.GetDeltaTime();

	//Advection
	double start = glfwGetTime();
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

void MACGrid3D::AdvectCellVelocity(double deltaTime)
{
	tbb::parallel_for(0, mNumCells, 1, [&](int index)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		double xVelocity = mCellXVelocities[index];
		double yVelocity = mCellYVelocities[index];
		double zVelocity = mCellZVelocities[index];

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

		glm::dvec3 avgVelocity(xVelocity, yVelocity, zVelocity);

		// Want to find previous position, so go backwards in time.
		glm::dvec3 prevPosition = mCellCenters[index] - avgVelocity * deltaTime;

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

int MACGrid3D::GetClosestCell(const glm::dvec3& inPosition)
{
	int x = static_cast<int>(floor((inPosition.x - dLeft - (mCellSize * 0.5f)) * mInvCellSize));
	int y = static_cast<int>(floor((inPosition.y - dBottom - (mCellSize * 0.5f)) * mInvCellSize));
	int z = static_cast<int>(floor((inPosition.z - dBack - (mCellSize * 0.5f)) * mInvCellSize));

	int approxIndex = GetIndexFromXYZ(x, y, z);

	glm::dvec3 vec = inPosition - mCellCenters[approxIndex];

	double distSqr = (vec.x * vec.x) + (vec.y * vec.y) + (vec.z * vec.z);

	double closestDistSqr = distSqr;
	int closestIndex = approxIndex;

	if (x > 0)
	{
		int neighbourLeft = GetIndexFromXYZ(x - 1, y, z);

		glm::dvec3 cToP = inPosition - mCellCenters[neighbourLeft];

		double dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourLeft;
		}
	}
	if (x < mNumCellWidth - 1)
	{
		int neighbourRight = GetIndexFromXYZ(x + 1, y, z);

		glm::dvec3 cToP = inPosition - mCellCenters[neighbourRight];

		double dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourRight;
		}
	}

	if (y > 0)
	{
		int neighbourBottom = GetIndexFromXYZ(x, y - 1, z);

		glm::dvec3 cToP = inPosition - mCellCenters[neighbourBottom];

		double dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourBottom;
		}
	}
	if (y < mNumCellHeight - 1)
	{
		int neighbourTop = GetIndexFromXYZ(x, y + 1, z);

		glm::dvec3 cToP = inPosition - mCellCenters[neighbourTop];

		double dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourTop;
		}
	}

	if (z > 0)
	{
		int neighbourBack = GetIndexFromXYZ(x, y, z - 1);

		glm::dvec3 cToP = inPosition - mCellCenters[neighbourBack];

		double dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourBack;
		}
	}
	if (z < mNumCellWidth - 1)
	{
		int neighbourFront = GetIndexFromXYZ(x, y, z + 1);

		glm::dvec3 cToP = inPosition - mCellCenters[neighbourFront];

		double dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = neighbourFront;
		}
	}

	return closestIndex;
}

void MACGrid3D::UpdateCellVelocity(double deltaTime)
{
	double scale = deltaTime * mInvCellSize * (1.0 / mDensity);

	double solidXVel = 0.0;
	double solidYVel = 0.0;
	double solidZVel = 0.0;

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

void MACGrid3D::CalculateCellDivergence()
{
	double scale = mInvCellSize;

	double solidXVel = 0.0;
	double solidYVel = 0.0;
	double solidZVel = 0.0;

	mCellDivergence.assign(mNumCells, 0.f);

	// Calculate negative divergence of each cell.
	for (int index = 0; index < mNumCells; index++)
	{
		if (mCellType[index] == CellType::eFLUID)
		{
			int x, y, z;
			std::tie(x, y, z) = GetXYZFromIndex(index);

			double divergence = 0.f;

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

void MACGrid3D::UpdateCellPressure(double deltaTime, int maxIterations)
{
	std::vector<double> Adiagonal;
	std::vector<double> Ax;
	std::vector<double> Ay;
	std::vector<double> Az;

	Adiagonal.assign(mNumCells, 0.0);
	Ax.assign(mNumCells, 0.0);
	Ay.assign(mNumCells, 0.0);
	Az.assign(mNumCells, 0.0);

	std::cout << "--------------------------------\n";
	double start = glfwGetTime();
	std::cout << "Initialize linear system: ";

	InitializeLinearSystem(deltaTime, Adiagonal, Ax, Ay, Az);

	std::cout << glfwGetTime() - start << "\n";

	// Preconditioned Conjugate Gradient.
	std::vector<double> newPressure;
	newPressure.assign(mNumCells, 0.0);

	start = glfwGetTime();
	std::cout << "Calculate cell divergence: ";

	CalculateCellDivergence();

	std::cout << glfwGetTime() - start << "\n";
	
	std::vector<double> residuals;
	residuals = mCellDivergence;

	double maxResidual = -1.0;

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


	std::vector<double> z;
	z.assign(mNumCells, 0.0);

	start = glfwGetTime();
	std::cout << "Calculate preconditioner: ";

	std::vector<double> precon;
	CalculatePreconditioner(precon, Adiagonal, Ax, Ay, Az);

	std::cout << glfwGetTime() - start << "\n";
	start = glfwGetTime();
	std::cout << "Apply preconditioner: ";

	ApplyPreconditioner(z, residuals, precon, Ax, Ay, Az);

	std::cout << glfwGetTime() - start << "\n";

	std::vector<double> search = z;

	double theta = 0.0;
	for (int i = 0; i < mNumCells; i++)
	{
		theta += z[i] * residuals[i];
	}

	int iteration;

	for (iteration = 1; iteration <= maxIterations; ++iteration)
	{
		ApplyA(deltaTime, z, search, Adiagonal, Ax, Ay, Az);

		double phi = 0.0;
		for (int i = 0; i < mNumCells; i++)
		{
			phi += z[i] * search[i];
		}

		if (abs(phi) < TOLERANCE)
		{
			phi = TOLERANCE;
		}

		double alpha = theta / phi;

		for (int i = 0; i < mNumCells; i++)
		{
			newPressure[i] += alpha * search[i];
			double diff = alpha * z[i];
			residuals[i] -= diff;
		}

		maxResidual = -1.0;

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

		double thetaNew = 0.0;
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
		double beta = thetaNew / theta;

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

void MACGrid3D::InitializeLinearSystem(double deltaTime, std::vector<double>& inDiag, std::vector<double>& inX, std::vector<double>& inY, std::vector<double>& inZ)
{
	double scale = -deltaTime * mInvCellSize * mInvCellSize * (1.0 / mDensity);

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

void MACGrid3D::ApplyA(double deltaTime, std::vector<double>& outResult, const std::vector<double>& inVec, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY, const std::vector<double>& inZ)
{
	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		double value = 0.f;

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

void MACGrid3D::CalculatePreconditioner(std::vector<double>& inOutPrecon, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY, const std::vector<double>& inZ)
{
	inOutPrecon.assign(mNumCells, 0.f);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			double Axi = 0.0; // Aplusi_iminus1
			double Axj = 0.0; // Aplusi_jminus1
			double Axk = 0.0; // Aplusi_kminus1

			double Ayi = 0.0; // Aplusj_iminus1
			double Ayj = 0.0; // Aplusj_jminus1
			double Ayk = 0.0; // Aplusj_kminus1

			double Azi = 0.0; // Aplusk_iminus1
			double Azj = 0.0; // Aplusk_jminus1
			double Azk = 0.0; // Aplusk_kminus1

			double preconi = 0.0; // precon_iminus1
			double preconj = 0.0; // precon_jminus1
			double preconk = 0.0; // precon_kminus1

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

			double a = Axi * preconi;
			double b = Ayj * preconj;
			double c = Azk * preconk;

			double termOne = a * a + b * b + c * c;

			double d = Axi * (Ayi + Azi) * preconi * preconi;
			double e = Ayj * (Axj + Azj) * preconj * preconj;
			double f = Azk * (Axk + Ayk) * preconk * preconk;

			double termTwo = d + e + f;

			double newPrecon = inDiag[index] - termOne - PRECON_TUNER * termTwo;

			if (newPrecon < PRECON_SAFETY * inDiag[index])
			{
				newPrecon = inDiag[index];
			}
					
			if (newPrecon > TOLERANCE)
			{
				inOutPrecon[index] = 1.0 / sqrt(newPrecon);
			}	
		}
	}
}

void MACGrid3D::ApplyPreconditioner(std::vector<double>& outResult, const std::vector<double>& inResidual, const std::vector<double>& inPrecon, const std::vector<double>& inX, const std::vector<double>& inY, const std::vector<double>& inZ)
{
	std::vector<double> intermediate;  // q
	intermediate.assign(mNumCells, 0.f);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{

			double Axi = 0.0; // Aplusi_iminus1
			double Ayj = 0.0; // Aplusj_jminus1
			double Azk = 0.0; // Aplusk_kminus1

			double preconi = 0.0; // precon_iminus1
			double preconj = 0.0; // precon_jminus1
			double preconk = 0.0; // precon_kminus1

			double intermediatei = 0.0; // qminusi
			double intermediatej = 0.0; // qminusj
			double intermediatek = 0.0; // qminusk

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

			double t = inResidual[index] - Axi * preconi * intermediatei
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
			double zi = 0.0;
			double zj = 0.0;
			double zk = 0.0;

			int neighbourRight = GetIndexFromXYZ(x + 1, y, z);

			zi = outResult[neighbourRight];

			int neighbourTop = GetIndexFromXYZ(x, y + 1, z);

			zj = outResult[neighbourTop];

			int neighbourFront= GetIndexFromXYZ(x, y, z + 1);

			zk = outResult[neighbourFront];

			double t = intermediate[index] - inX[index] * inPrecon[index] * zi
										  - inY[index] * inPrecon[index] * zj
										  - inZ[index] * inPrecon[index] * zk;

			outResult[index] = t * inPrecon[index];
		}
	}
}

const CellType MACGrid3D::GetCellTypeFromPosition(const glm::dvec3& inPosition)
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