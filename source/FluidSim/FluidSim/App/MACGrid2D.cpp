#include "MACGrid2D.h"

#include <iostream>
#include "GLFW/glfw3.h"

#define TOLERANCE 0.000001f
#define PRECON_TUNER 0.97f
#define PRECON_SAFETY 0.25f

//using namespace oneapi;

MACGrid2D::MACGrid2D(const ApplicationData& inData)
{
	InitializeGrid(inData);
	InitializeCellsFromParticles(inData.Get2DParticlePositions());
}

void MACGrid2D::InitializeGrid(const ApplicationData& inData)
{
	dLeft = inData.GetGridLeft();
	dBottom = inData.GetGridBottom();

	float dWidth = inData.GetGridWidth();
	float dHeight = inData.GetGridHeight();

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

	mCellPressures.assign(mNumCells, 0.f);
	mCellXVelocities.assign(mNumCells, 0.f);
	mCellYVelocities.assign(mNumCells, 0.f);
	mCellDivergence.assign(mNumCells, 0.f);
	mCellType.assign(mNumCells, CellType::eNONE);

	mIntXVelocities.assign(mNumCells, 0.f);
	mIntYVelocities.assign(mNumCells, 0.f);

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
		mCellPressures[cellIndex] = 10.0f;
	});
}

void MACGrid2D::Update(ApplicationData& inOutData)
{
	float deltaTime = inOutData.GetDeltaTime();

	//Advection
	float start = glfwGetTime();
	AdvectCellVelocity(deltaTime);

	std::cout << "MAC: advect: " << glfwGetTime() - start << "\n";

	// Projection step
	start = glfwGetTime();

//	UpdateCellPressure(deltaTime, 100);

	std::cout << "MAC: Pressure: " << glfwGetTime() - start << "\n";

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
	int x = round(tempA);
	int y = round(tempB);

	int approxIndex = GetIndexFromXY(x, y);

	if (approxIndex >= mNumCells)
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

	mCellDivergence.assign(mNumCells, 0.f);

	// Calculate negative divergence of each cell.
	for (int index = 0; index < mNumCells; index++)
	{
		if (mCellType[index] == CellType::eFLUID)
		{
			int x, y;
			std::tie(x, y) = GetXYFromIndex(index);

			float divergence = 0.f;

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
	std::vector<float> Adiagonal;
	std::vector<float> Ax;
	std::vector<float> Ay;

	Adiagonal.assign(mNumCells, 0.0f);
	Ax.assign(mNumCells, 0.0f);
	Ay.assign(mNumCells, 0.0f);

	std::cout << "--------------------------------\n";
	float start = glfwGetTime();
	std::cout << "Initialize linear system: ";

	InitializeLinearSystem(deltaTime, Adiagonal, Ax, Ay);

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
			maxResidual = abs(residuals[index]);
		}
	}

	if (maxResidual < TOLERANCE)
	{
		// Reached convergence;
		std::cout << "Zero cell divergence. Skipping pressure update. \n";
		return;
	}


	std::vector<float> z;
	z.assign(mNumCells, 0.0f);

	start = glfwGetTime();
	std::cout << "Calculate preconditioner: ";

	std::vector<float> precon;
	CalculatePreconditioner(precon, Adiagonal, Ax, Ay);

	std::cout << glfwGetTime() - start << "\n";
	start = glfwGetTime();
	std::cout << "Apply preconditioner: ";

	ApplyPreconditioner(z, residuals, precon, Ax, Ay);

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
		ApplyA(deltaTime, z, search, Adiagonal, Ax, Ay);

		float phi = 0.0f;
		for (int i = 0; i < mNumCells; i++)
		{
			phi += z[i] * search[i];
		}

		float alpha = theta / phi;

		for (int i = 0; i < mNumCells; i++)
		{
			newPressure[i] += alpha * search[i];
			residuals[i] -= alpha * z[i];
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

		ApplyPreconditioner(z, residuals, precon, Ax, Ay);

		float thetaNew = 0.0f;
		for (int i = 0; i < mNumCells; i++)
		{
			thetaNew += z[i] * residuals[i];
		}

		float beta = thetaNew / theta;

		for (int i = 0; i < mNumCells; i++)
		{
			search[i] = z[i] + beta * search[i];
		}

		theta = thetaNew;
	}

	if (iteration >= maxIterations)
	{
		std::cout << "WARNING: MAX NUMBER OF ITERATIONS REACHED IN PRESSURE SOLVE" << "\n";
		std::cout << "Check pressure solver for potential issues, MACGrid2D.cpp" << "\n";
	}
	std::cout << "NUM PRESSURE SOLVE ITERATIONS: " << iteration << "\n";
	std::cout << "------------------------------\n";

	// Set pressure
	mCellPressures = newPressure;
}

void MACGrid2D::InitializeLinearSystem(float deltaTime, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY)
{
	float scale = -deltaTime * mInvCellSize * mInvCellSize * (1.0f / mDensity);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y; 
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			if (x > 0)
			{
				int neighbourLeft = GetIndexFromXY(x - 1, y);

				if (mCellType[neighbourLeft] == CellType::eFLUID)
				{
					inDiag[index] += scale;
				}
			}

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

			if (y > 0)
			{
				int neighbourBottom = GetIndexFromXY(x, y - 1);

				if (mCellType[neighbourBottom] == CellType::eFLUID)
				{
					inDiag[index] += scale;
				}
			}

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

void MACGrid2D::ApplyA(float deltaTime, std::vector<float>& outResult, const std::vector<float>& inVec, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY)
{
	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		float value = 0.f;

		if (x > 0)
		{
			int neighbourLeft = GetIndexFromXY(x - 1, y);
			value += inVec[neighbourLeft] * inX[neighbourLeft];
		}

		if (x < mNumCellWidth - 1)
		{
			int neighbourRight = GetIndexFromXY(x + 1, y);  
			value += inVec[neighbourRight] * inX[index];
		}

		if (y > 0)
		{
			int neighbourBottom = GetIndexFromXY(x, y - 1);
			value += inVec[neighbourBottom] * inY[neighbourBottom];
		}

		if (y < mNumCellHeight - 1)
		{
			int neighbourTop = GetIndexFromXY(x, y + 1);
			value += inVec[neighbourTop] * inY[index];
		}

  		value += inDiag[index] * inVec[index];
		outResult[index] = value;
	}
}

void MACGrid2D::CalculatePreconditioner(std::vector<float>& inOutPrecon, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY)
{
	inOutPrecon.assign(mNumCells, 0.f);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			float Axi = 0.0f; // Aplusi_iminus1
			float Axj = 0.0f; // Aplusi_jminus1

			float Ayi = 0.0f; // Aplusj_iminus1
			float Ayj = 0.0f; // Aplusj_jminus1

			float preconi = 0.0f; // precon_iminus1
			float preconj = 0.0f; // precon_jminus1

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

			float a = Axi * preconi;
			float b = Ayj * preconj;

			float termOne = a * a + b * b;

			float d = Axi * Ayi * preconi * preconi;
			float e = Ayj * Axj * preconj * preconj;

			float termTwo = d + e;

			float newPrecon = - inDiag[index] - termOne - PRECON_TUNER * termTwo;

			if (newPrecon < PRECON_SAFETY * inDiag[index])
			{
				newPrecon = - inDiag[index];
			}

			inOutPrecon[index] = 1.0f / sqrt(newPrecon);
		}
	}
}

void MACGrid2D::ApplyPreconditioner(std::vector<float>& outResult, const std::vector<float>& inResidual, const std::vector<float>& inPrecon, const std::vector<float>& inX, const std::vector<float>& inY)
{
	std::vector<float> intermediate;  // q
	intermediate.assign(mNumCells, 0.f);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{

			float Axi = 0.0f; // Aplusi_iminus1
			float Ayj = 0.0f; // Aplusj_jminus1

			float preconi = 0.0f; // precon_iminus1
			float preconj = 0.0f; // precon_jminus1

			float intermediatei = 0.0f; // qminusi
			float intermediatej = 0.0f; // qminusj

			int neighbourLeft = GetIndexFromXY(x - 1, y);

			Axi = inX[neighbourLeft];
			preconi = inPrecon[neighbourLeft];
			intermediatei = intermediate[neighbourLeft];

			int neighbourBottom = GetIndexFromXY(x, y - 1);

			Ayj = inY[neighbourBottom];
			preconj = inPrecon[neighbourBottom];
			intermediatej = intermediate[neighbourBottom];

			float t = inResidual[index] - Axi * preconi * intermediatei
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
			float zi = 0.0f;
			float zj = 0.0f;

			int neighbourRight = GetIndexFromXY(x + 1, y);

			zi = outResult[neighbourRight];

			int neighbourTop = GetIndexFromXY(x, y + 1);

			zj = outResult[neighbourTop];

			float t = intermediate[index] - inX[index] * inPrecon[index] * zi
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