#include "MACGrid2D.h"

#include <iostream>
#include "GLFW/glfw3.h"

#define TOLERANCE 0.000001
#define PRECON_TUNER 0.97
#define PRECON_SAFETY 0.25

MACGrid2D::MACGrid2D(const ApplicationData& inData)
{
	InitializeGrid(inData);
	InitializeCellsFromParticles(inData.Get2DParticlePositions());
}

void MACGrid2D::InitializeGrid(const ApplicationData& inData)
{
	dLeft = inData.GetGridLeft();
	dBottom = inData.GetGridBottom();

	mNumCellWidth = inData.GetNumGridCellsWidth();
	mNumCellHeight = inData.GetNumGridCellsHeight();

	mNumCells = inData.GetNumGridCells();

	mCellSize = inData.GetGridCellSize();
	mInvCellSize = 1 / mCellSize;

	float halfCell = mCellSize * 0.5f;

	mCellCenters.reserve(mNumCells);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		float centerX = dLeft + (x * mCellSize) + halfCell;
		float centerY = dBottom + (y * mCellSize) + halfCell;

		mCellCenters.push_back(glm::vec2(centerX, centerY));
	}

	mCellPressures.assign(mNumCells, 0.0);
	mCellXVelocities.assign(mNumCells, 0.0);
	mCellYVelocities.assign(mNumCells, 0.0);
	mCellDivergence.assign(mNumCells, 0.0);
	mCellType.assign(mNumCells, CellType::eNONE);

	mIntXVelocities.assign(mNumCells, 0.0);
	mIntYVelocities.assign(mNumCells, 0.0);

	mDensity = inData.GetFluidDensity();
}

void MACGrid2D::InitializeCellsFromParticles(const std::vector<glm::vec2>& inParticlePositions)
{
	// Set edge cells to be solid.
	tbb::parallel_for(0, mNumCells, 1, [&](int cIndex)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(cIndex);

		if (x == 0 || x == mNumCellWidth - 1 || y == 0 || y == mNumCellHeight - 1 )
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

void MACGrid2D::Update(ApplicationData& inOutData)
{
	float deltaTime = inOutData.GetDeltaTime();

	//Advection
	double start = glfwGetTime();
	AdvectCellVelocity(deltaTime);

	std::cout << "MAC: advect: " << glfwGetTime() - start << "\n";

	// Calculate cell divergence
	start = glfwGetTime();

	CalculateCellDivergence(deltaTime);

	std::cout << "Calculate cell divergence: " << glfwGetTime() - start << "\n";

	// Projection step
	start = glfwGetTime();
	
	UpdateCellPressure(deltaTime, 200);

	std::cout << "MAC: Pressure: " << glfwGetTime() - start << "\n";
	std::cout << "------------------------------\n";

	// Velocity update
	start = glfwGetTime();
	UpdateCellVelocity(deltaTime);

	std::cout << "MAC: cell update: " << glfwGetTime() - start << "\n";

	inOutData.SetCellTypes(mCellType);
}

void MACGrid2D::AdvectCellVelocity(float deltaTime)
{
	tbb::parallel_for(0, mNumCells, 1, [&](int index)
	{
			if (mCellType[index] == CellType::eFLUID)
			{
				int x, y;
				std::tie(x, y) = GetXYFromIndex(index);

				glm::vec2 cellVelocity = { mCellXVelocities[index], mCellYVelocities[index] };

				//Runge Kutta
				// First take a half step backwards over half the timestep to find the average velocity
				glm::vec2 halfPrevPosition = mCellCenters[index] - cellVelocity * 0.5f * deltaTime;

				int midStepIndex = GetClosestCell(halfPrevPosition);

				int midStepX, midStepY;
				std::tie(midStepX, midStepY) = GetXYFromIndex(midStepIndex);
				  
				int midStepNeighbourRight = GetIndexFromXY(midStepX + 1, midStepY);
				int midStepNeighbourTop = GetIndexFromXY(midStepX, midStepY + 1);

				// Compute the avgerage velocity by interpolating the velocity of the mid step position.
				glm::vec2 diff = mCellCenters[midStepIndex] - halfPrevPosition;

				glm::dvec2 weight = (diff * mInvCellSize) + 0.5f;

				double avgStepXVelocity = (mCellXVelocities[midStepIndex] * weight.x) + (mCellXVelocities[midStepNeighbourRight] * (1 - weight.x));
				double avgStepYVelocity = (mCellYVelocities[midStepIndex] * weight.y) + (mCellYVelocities[midStepNeighbourTop] * (1 - weight.y));

				glm::vec2 avgStepVelocity = { avgStepXVelocity, avgStepYVelocity };

				// Now that we have the average velocity of the step, we can work out the original position at the beginning of the step.
				glm::vec2 prevPosition = mCellCenters[index] - avgStepVelocity * deltaTime;

				int prevCellIndex = GetClosestCell(prevPosition);

				int prevX, prevY;
				std::tie(prevX, prevY) = GetXYFromIndex(prevCellIndex);

				int prevNeighbourRight = GetIndexFromXY(prevX + 1, prevY);
				int prevNeighbourTop = GetIndexFromXY(prevX, prevY + 1);


				if (prevCellIndex >= 0)
				{
					// Recalculate our weights for interpolation
					diff = mCellCenters[prevCellIndex] - prevPosition;

					weight = (diff * mInvCellSize) + 0.5f;

					double prevXVelocity = (mCellXVelocities[prevCellIndex] * weight.x) + (mCellXVelocities[prevNeighbourRight] * (1 - weight.x));
					double prevYVelocity = (mCellYVelocities[prevCellIndex] * weight.y) + (mCellYVelocities[prevNeighbourTop] * (1 - weight.y));

					mIntXVelocities[index] = prevXVelocity;
					mIntYVelocities[index] = prevYVelocity;
				}
			}
	});
}

int MACGrid2D::GetClosestCell(const glm::vec2& inPos)
{
	float tempA = (inPos.x - dLeft - (mCellSize * 0.5f)) * mInvCellSize;
	float tempB = (inPos.y - dBottom - (mCellSize * 0.5f)) * mInvCellSize;
	int x = static_cast<int>(round(tempA));
	int y = static_cast<int>(round(tempB));

	int approxIndex = GetIndexFromXY(x, y);

	if (approxIndex >= mNumCells || approxIndex < 0)
	{
		approxIndex = -1;
	}

	return approxIndex;
}

void MACGrid2D::UpdateCellVelocity(float deltaTime)
{
	float scale = deltaTime * mInvCellSize * (1.0f / mDensity);

	float solidXVel = 0.0f;
	float solidYVel = 0.0f;
	float solidZVel = 0.0f;

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (x > 0)
		{
			int neighbourLeft = GetIndexFromXY(x - 1, y);

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
		}
		if (y > 0)
		{
			int neighbourBottom = GetIndexFromXY(x, y - 1);

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
		}
	}
}

void MACGrid2D::CalculateCellDivergence(float deltaTime)
{
	float scale = mInvCellSize;

	float solidXVel = 0.0f;
	float solidYVel = 0.0f;

	// Calculate negative divergence of each cell.
	for (int index = 0; index < mNumCells; index++)
	{
		if (mCellType[index] == CellType::eFLUID)
		{
			int x, y;
			std::tie(x, y) = GetXYFromIndex(index);

			double divergence = 0.0;

			int neighbourRight = GetIndexFromXY(x + 1, y);
			divergence += mIntXVelocities[neighbourRight] - mIntXVelocities[index];

			int neighbourTop = GetIndexFromXY(x, y + 1);
			divergence += mIntYVelocities[neighbourTop] - mIntYVelocities[index];

			divergence *= -scale;

			mCellDivergence[index] = divergence;
		}
	}

	// Update divergence to account for solid cells.
	for (int index = 0; index < mNumCells; index++)
	{
		if (mCellType[index] == CellType::eFLUID)
		{
			int x, y;
			std::tie(x, y) = GetXYFromIndex(index);

			int neighbourLeft = GetIndexFromXY(x - 1, y);

			if (mCellType[neighbourLeft] == CellType::eSOLID)
			{
				mCellDivergence[index] -= scale * (mIntXVelocities[index] - solidXVel);
			}

			int neighbourRight = GetIndexFromXY(x + 1, y);

			if (mCellType[neighbourRight] == CellType::eSOLID)
			{
				mCellDivergence[index] += scale * (mIntXVelocities[neighbourRight] - solidXVel);
			}

			int neighbourBottom = GetIndexFromXY(x, y - 1);

			if (mCellType[neighbourBottom] == CellType::eSOLID)
			{
				mCellDivergence[index] -= scale * (mIntYVelocities[index] - solidYVel);
			}

			int neighbourTop = GetIndexFromXY(x, y + 1);

			if (mCellType[neighbourTop] == CellType::eSOLID)
			{
				mCellDivergence[index] += scale * (mIntYVelocities[neighbourTop] - solidYVel);
			}
		}
	}
}

void MACGrid2D::UpdateCellPressure(float deltaTime, int maxIterations)
{
	std::vector<double> Adiagonal;
	std::vector<double> Ax;
	std::vector<double> Ay;

	Adiagonal.assign(mNumCells, 0.0);
	Ax.assign(mNumCells, 0.0);
	Ay.assign(mNumCells, 0.0);

	std::cout << "--------------------------------\n";
	double start = glfwGetTime();
	std::cout << "Initialize linear system: ";

	InitializeLinearSystem(deltaTime, Adiagonal, Ax, Ay);

	std::cout << glfwGetTime() - start << "\n";

	// Preconditioned Conjugate Gradient.
	Eigen::VectorXd newPressure(mNumCells);
	newPressure.fill(0.0);

	Eigen::VectorXd residuals(mNumCells);
	// TO DO: make cell divergence an eigen vector.
	for (int i = 0; i < mNumCells; i++)
	{
		residuals[i] = mCellDivergence[i];
	}

	double maxResidual = -1.0;

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
		std::cout << "Zero cell divergence. Skipping pressure update. \n";
		return;
	}

	Eigen::VectorXd z(mNumCells);
	z.fill(0.0);

	start = glfwGetTime();
	std::cout << "Calculate preconditioner: ";

	std::vector<double> precon;
	precon.assign(mNumCells, 0.0);
	CalculatePreconditioner(precon, Adiagonal, Ax, Ay);

	std::cout << glfwGetTime() - start << "\n";
	start = glfwGetTime();
	std::cout << "Apply preconditioner: ";

	ApplyPreconditioner(z, residuals, precon, Ax, Ay);

	std::cout << glfwGetTime() - start << "\n";

	Eigen::VectorXd search = z;

	double sigma = z.dot(residuals);

	int iteration;

	for (iteration = 0; iteration < maxIterations; ++iteration)
	{
		ApplyA(deltaTime, z, search, Adiagonal, Ax, Ay);

		double phi = z.dot(search);
		double alpha = sigma / phi;

		// Update pressure and residual vectors.
		residuals -= z * alpha;
		newPressure += search * alpha;

		// Calculate maximum residual error.
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

		ApplyPreconditioner(z, residuals, precon, Ax, Ay);

		// SigmaNew = dot product of z and residuals.

		double sigmaNew = z.dot(residuals);
		double beta = sigmaNew / sigma;

		// Update search vector for next iteration.
		search *= beta;
		search += z;

		sigma = sigmaNew;
	}
	std::cout << "------------------------------\n";
	if (iteration >= maxIterations)
	{
		std::cout << "WARNING: MAX NUMBER OF ITERATIONS REACHED IN PRESSURE SOLVE" << "\n";
		std::cout << "Check pressure solver for potential issues, MACGrid2D.cpp" << "\n";
	}

	std::cout << "NUM PRESSURE SOLVE ITERATIONS: " << iteration << "\n";
	std::cout << "PRESSURE SOLVE ERROR: " << maxResidual << "\n";

	// Set pressure
	for (int i = 0; i < mNumCells; i++)
	{
		mCellPressures[i] = newPressure[i];
	}
}

void MACGrid2D::InitializeLinearSystem(float deltaTime, std::vector<double>& inDiag, std::vector<double>& inX, std::vector<double>& inY)
{
	double scale = deltaTime * mInvCellSize * mInvCellSize * (1.0f / mDensity);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y; 
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			// Handle left neighbour
			if (x > 0)
			{
				int neighbourLeft = GetIndexFromXY(x - 1, y);

				if (mCellType[neighbourLeft] == CellType::eFLUID)
				{
					inDiag[index] += scale;
				}
			}
			// Handle right neighbour
			if (x < mNumCellWidth - 1)
			{
				int neighbourRight = GetIndexFromXY(x + 1, y);

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
			// Handle bottom neighbour
			if (y > 0)
			{
				int neighbourBottom = GetIndexFromXY(x, y - 1);

				if (mCellType[neighbourBottom] == CellType::eFLUID)
				{
					inDiag[index] += scale;
				}
			}
			// Handle top neighbour
			if (y < mNumCellHeight - 1)
			{
				int neighbourTop = GetIndexFromXY(x, y + 1);

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
		}
	}
}

void MACGrid2D::ApplyA(float deltaTime, Eigen::VectorXd& outResult, const Eigen::VectorXd& inVec, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY)
{
	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		double value = 0.0;

		if (x > 0)
		{
			int neighbourLeft = GetIndexFromXY(x - 1, y);
			//value += inVec[neighbourLeft] *inX[neighbourLeft];
			value += inVec[neighbourLeft] * inX[index];
		}

		if (x < mNumCellWidth - 1)
		{
			int neighbourRight = GetIndexFromXY(x + 1, y);  
			//value += inVec[neighbourRight] *inX[index];
			value += inVec[neighbourRight] * inX[neighbourRight];

		}

		if (y > 0)
		{
			int neighbourBottom = GetIndexFromXY(x, y - 1);
			//value += inVec[neighbourBottom] *inY[neighbourBottom];
			value += inVec[neighbourBottom] * inY[index];
		}

		if (y < mNumCellHeight - 1)
		{
			int neighbourTop = GetIndexFromXY(x, y + 1);
			//value += inVec[neighbourTop] *inY[index];
			value += inVec[neighbourTop] * inY[neighbourTop];
		}

  		value += inDiag[index] * inVec[index];
		outResult[index] = value;
	}
}

void MACGrid2D::CalculatePreconditioner(std::vector<double>& inOutPrecon, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY)
{
	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			double Axi = 0.0f; // Aplusi_iminus1
			double Axj = 0.0f; // Aplusi_jminus1

			double Ayi = 0.0f; // Aplusj_iminus1
			double Ayj = 0.0f; // Aplusj_jminus1

			double preconi = 0.0f; // precon_iminus1
			double preconj = 0.0f; // precon_jminus1

			if (x > 0)
			{
				int neighbourLeft = GetIndexFromXY(x - 1, y);

				Axi = inX[neighbourLeft];
				Ayi = inY[neighbourLeft];

				preconi = inOutPrecon[neighbourLeft];
			}

			if (y > 0)
			{
				int neighbourBottom = GetIndexFromXY(x, y - 1);

				Axj = inX[neighbourBottom];
				Ayj = inY[neighbourBottom];

				preconj = inOutPrecon[neighbourBottom];
			}

			double a = Axi * preconi;
			double b = Ayj * preconj;

			double termOne = a * a + b * b;

			double d = Axi * Ayi * preconi * preconi;
			double e = Ayj * Axj * preconj * preconj;

			double termTwo = d + e;

			double newPrecon = inDiag[index] - termOne - PRECON_TUNER * termTwo;

			if (newPrecon < PRECON_SAFETY * inDiag[index])
			{
				newPrecon =  inDiag[index];
			}

			inOutPrecon[index] = 1.0 / sqrt(newPrecon);
		}
	}
}

void MACGrid2D::ApplyPreconditioner(Eigen::VectorXd& outResult, const Eigen::VectorXd& inResidual, const std::vector<double>& inPrecon, const std::vector<double>& inX, const std::vector<double>& inY)
{
	std::vector<double> intermediate;  // q
	intermediate.assign(mNumCells, 0.f);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{

			double Axi = 0.0f; // Aplusi_iminus1
			double Ayj = 0.0f; // Aplusj_jminus1

			double preconi = 0.0f; // precon_iminus1
			double preconj = 0.0f; // precon_jminus1

			double intermediatei = 0.0f; // qminusi
			double intermediatej = 0.0f; // qminusj

			int neighbourLeft = GetIndexFromXY(x - 1, y);

			Axi = inX[neighbourLeft];
			preconi = inPrecon[neighbourLeft];
			intermediatei = intermediate[neighbourLeft];

			int neighbourBottom = GetIndexFromXY(x, y - 1);

			Ayj = inY[neighbourBottom];
			preconj = inPrecon[neighbourBottom];
			intermediatej = intermediate[neighbourBottom];

			double t = inResidual[index] - Axi * preconi * intermediatei
										- Ayj * preconj * intermediatej;

			intermediate[index] = t * inPrecon[index];
		}
	}

	for (int index = mNumCells - 1; index >= 0; index--)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			double zi = 0.0f;
			double zj = 0.0f;

			int neighbourRight = GetIndexFromXY(x + 1, y);

			zi = outResult[neighbourRight];

			int neighbourTop = GetIndexFromXY(x, y + 1);

			zj = outResult[neighbourTop];

			double t = intermediate[index] - inX[index] * inPrecon[index] * zi
										  - inY[index] * inPrecon[index] * zj;

			outResult[index] = t * inPrecon[index];
		}
	}
}

const CellType MACGrid2D::GetCellTypeFromPosition(const glm::vec2& inPos)
{
	int index = GetClosestCell(inPos);

	return GetCellType(index);
}

std::tuple<int, int> MACGrid2D::GetXYFromIndex(int index)
{
	int x = 0;
	int y = 0;

	y = index % mNumCellHeight;
	x = static_cast<int>(round((index - y) / mNumCellHeight)) % mNumCellWidth;

	return std::tuple<int, int>(x, y);
}

int MACGrid2D::GetIndexFromXY(int X, int Y)
{
	return  Y + X * mNumCellHeight;
}