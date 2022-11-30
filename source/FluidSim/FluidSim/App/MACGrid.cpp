#include "MACGrid.h"

#include <iostream>

#define TOLERANCE 0.000001f
#define PRECON_TUNER 0.97f
#define PRECON_SAFETY 0.25f

MACGrid::MACGrid(const Domain& inDomain, const std::vector<glm::vec3>& inParticlePositions, int inGridResolution)
{
	InitializeFromDomain(inDomain, inGridResolution);
	InitializeCellsFromParticles(inParticlePositions);
}

void MACGrid::InitializeFromDomain(const Domain& inDomain, int inGridResolution)
{
	glm::vec3 dCenter = inDomain.GetCenter();

	float dLeft = inDomain.GetLeft();
	float dBack = inDomain.GetBack();
	float dBottom = inDomain.GetBottom();

	float dLength = inDomain.GetLength();
	float dWidth = inDomain.GetWidth();
	float dHeight = inDomain.GetHeight();

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
	mFluidCell.assign(mNumCells, false);

	mIntXVelocities.assign(mNumCells, 0.f);
	mIntYVelocities.assign(mNumCells, 0.f);
	mIntZVelocities.assign(mNumCells, 0.f);
}

void MACGrid::InitializeCellsFromParticles(const std::vector<glm::vec3>& inParticlePositions)
{
	for (int cIndex = 0; cIndex < GetNumCells(); cIndex++)
	{
		for (int pIndex = 0; pIndex < inParticlePositions.size(); pIndex++)
		{
			glm::vec3 cToP = inParticlePositions[pIndex] - mCellCenters[cIndex];

			if (glm::length(cToP) < 1.0f)
			{
				mFluidCell[cIndex] = true;
				mCellPressures[cIndex] = 1.0f;

				break;
			}
		}
	}
}

void MACGrid::SetCellVelocity(int index, const glm::vec3& inVelocity)
{
	mCellXVelocities[index] = inVelocity.x;
	mCellYVelocities[index] = inVelocity.y;
	mCellZVelocities[index] = inVelocity.z;
}

void MACGrid::Update(float deltaTime)
{
	//Advection
	AdvectCellVelocity(deltaTime);

	// Projection step
	UpdateCellPressure(deltaTime, 200);

	UpdateCellVelocity(deltaTime);
}

void MACGrid::AdvectCellVelocity(float deltaTime)
{
	for (int index = 0; index < mNumCells; index++)
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
	}

}

int MACGrid::GetClosestCell(const glm::vec3& inPosition)
{
	float closestDistSqr = 10000.0f;
	int closestIndex = -1;

	for (int index = 0; index < mNumCells; index++)
	{
		glm::vec3 cToP = inPosition - mCellCenters[index];

		float dist = (cToP.x * cToP.x) + (cToP.y * cToP.y) + (cToP.z * cToP.z);

		if (dist < closestDistSqr)
		{
			closestDistSqr = dist;
			closestIndex = index;
		}
	}

	return closestIndex;
}

void MACGrid::UpdateCellVelocity(float deltaTime)
{
	float scale = deltaTime * mInvCellSize; // TO DO: water density

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (x < mNumCellWidth - 1)
		{
			int neighbourRight = GetIndexFromXYZ(x + 1, y, z);
					
			float velocity = mCellXVelocities[index];
					
			velocity -= scale * (mCellPressures[neighbourRight] - mCellPressures[index]);

			mCellXVelocities[index] = velocity;
		}

		if (y < mNumCellHeight - 1)
		{
			int neighbourTop = GetIndexFromXYZ(x, y + 1, z);

			float velocity = mCellYVelocities[index];

			velocity -= scale * (mCellPressures[neighbourTop] - mCellPressures[index]);

			mCellYVelocities[index] = velocity;
		}

		if (z < mNumCellLength - 1)
		{
			int neighbourFront = GetIndexFromXYZ(x, y, z + 1);

			float velocity = mCellZVelocities[index];

			velocity -= scale * (mCellPressures[neighbourFront] - mCellPressures[index]);

			mCellZVelocities[index] = velocity;
		}
	}
}

void MACGrid::CalculateCellDivergence(float deltaTime)
{
	float scale = mInvCellSize;

	mCellDivergence.assign(mNumCells, 0.f);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		float divergence = 0.f;

		if (x < mNumCellWidth - 1)
		{
			int neighbourRight = GetIndexFromXYZ(x + 1, y, z);

			//divergence += mCellXVelocities[neighbourRight] - mCellXVelocities[index];
			divergence += mIntXVelocities[neighbourRight] - mIntXVelocities[index];
		}

		if (y < mNumCellHeight - 1)
		{
			int neighbourTop = GetIndexFromXYZ(x, y + 1, z);

			//divergence += mCellYVelocities[neighbourTop] - mCellYVelocities[index];
			divergence += mIntYVelocities[neighbourTop] - mIntYVelocities[index];
		}

		if (z < mNumCellLength - 1)
		{
			int neighbourFront = GetIndexFromXYZ(x, y, z + 1);

			//divergence += mCellZVelocities[neighbourFront] - mCellZVelocities[index];
			divergence += mIntZVelocities[neighbourFront] - mIntZVelocities[index];
		}

		divergence *= -scale;

		mCellDivergence[index] = divergence;
	}
}

void MACGrid::UpdateCellPressure(float deltaTime, int maxIterations)
{
	int numCells = mNumCellHeight * mNumCellWidth * mNumCellLength;

	std::vector<float> Adiagonal;
	std::vector<float> Ax;
	std::vector<float> Ay;
	std::vector<float> Az;

	Adiagonal.assign(numCells, 0.0f);
	Ax.assign(numCells, 0.0f);
	Ay.assign(numCells, 0.0f);
	Az.assign(numCells, 0.0f);

	InitializeLinearSystem(deltaTime, Adiagonal, Ax, Ay, Az);

	// Preconditioned Conjugate Gradient.

	std::vector<float> newPressure;
	newPressure.assign(numCells, 0.0f);

	CalculateCellDivergence(deltaTime);
	std::vector<float> residuals;
	residuals = mCellDivergence;

	std::vector<float> z;
	z.assign(numCells, 0.0f);

	std::vector<float> precon;
	CalculatePreconditioner(precon, Adiagonal, Ax, Ay, Az);

	ApplyPreconditioner(z, residuals, precon, Ax, Ay, Az);

	std::vector<float> search = z;

	float theta = 0.0f;
	for (int i = 0; i < numCells; i++)
	{
		theta += z[i] * residuals[i];
	}

	int iteration;
	bool converged = false;

	for (iteration = 0; iteration < maxIterations; ++iteration)
	{
		ApplyA(deltaTime, z, search, Adiagonal, Ax, Ay, Az);

		float phi = 0.0f;
		for (int i = 0; i < numCells; i++)
		{
			phi += z[i] * search[i];
		}

		if (phi < TOLERANCE)
		{
			phi = TOLERANCE;
		}

		float alpha = theta / phi;

		for (int i = 0; i < numCells; i++)
		{
			newPressure[i] += alpha * search[i];
			residuals[i] -= alpha * z[i];
		}

		float maxResidual = -1.0f;

		for (int index = 0; index < numCells; index++)
		{
			if (abs(residuals[index]) > maxResidual)
			{
				maxResidual = residuals[index];
			}
		}

		if (maxResidual < TOLERANCE)
		{
			// Reached convergence;
			break;
		}

		ApplyPreconditioner(z, residuals, precon, Ax, Ay, Az);

		float thetaNew = 0.0f;
		for (int i = 0; i < numCells; i++)
		{
			thetaNew += z[i] * residuals[i];
		}

		float beta = thetaNew / theta;

		for (int i = 0; i < numCells; i++)
		{
			search[i] = z[i] + beta * search[i];
		}

		theta = thetaNew;

		
	}

	if (iteration == maxIterations)
	{
		std::cout << "WARNING: MAX NUMBER OF ITERATIONS REACHED IN PRESSURE SOLVE" << "\n";
		std::cout << "Check pressure solver for potential issues, MACGrid.cpp" << "\n";
	}
	std::cout << "NUM PRESSURE SOLVE ITERATIONS: " << iteration << "\n";

	// Set pressure

	for (int i = 0; i < numCells; i++)
	{
		//mGridCells[i].SetPressure(newPressure[i]);
		if (mFluidCell[i])
		{
			mCellPressures[i] = newPressure[i];
		}
		else
		{
			mCellPressures[i] = 0.1f;
		}
	}
}

void MACGrid::InitializeLinearSystem(float deltaTime, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ)
{
	float scale = deltaTime * mInvCellSize * mInvCellSize; // TO DO include fluid density.

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (mFluidCell[index])
		{
			if (x > 0)
			{
				int neighbourLeft = GetIndexFromXYZ(x - 1, y, z);

				if (mFluidCell[neighbourLeft])
				{
					inDiag[index] += scale;
				}
			}

			if (x < mNumCellWidth - 1)
			{
				int neighbourRight = GetIndexFromXYZ(x + 1, y, z);

				if (mFluidCell[neighbourRight])
				{
					inDiag[index] += scale;
					inX[index] -= scale;
				}
				else
				{
					inDiag[index] += scale;
				}
			}

			if (y > 0)
			{
				int neighbourBottom = GetIndexFromXYZ(x, y - 1, z);

				if (mFluidCell[neighbourBottom])
				{
					inDiag[index] += scale;
				}
			}

			if (y < mNumCellHeight - 1)
			{
				int neighbourTop = GetIndexFromXYZ(x, y + 1, z);

				if (mFluidCell[neighbourTop])
				{
					inDiag[index] += scale;
					inY[index] -= scale;
				}
				else
				{
					inDiag[index] += scale;
				}
			}

			if (z > 0)
			{
				int neighbourBack = GetIndexFromXYZ(x, y, z - 1);

				if (mFluidCell[neighbourBack])
				{
					inDiag[index] += scale;
				}
			}

			if (z < mNumCellWidth - 1)
			{
				int neighbourFront = GetIndexFromXYZ(x, y, z + 1);

				if (mFluidCell[neighbourFront])
				{
					inDiag[index] += scale;
					inZ[index] -= scale;
				}
				else
				{
					inDiag[index] += scale;
				}
			}
		}
	}
}

void MACGrid::ApplyA(float deltaTime, std::vector<float>& outResult, const std::vector<float>& inVec, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ)
{
	int numCells = mNumCellHeight * mNumCellWidth * mNumCellLength;
	float scale = deltaTime * mInvCellSize * mInvCellSize; // TO DO;

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

		value *= -scale;
		value += inDiag[index] * scale * inVec[index];

		outResult[index] = value;
	}
}

void MACGrid::CalculatePreconditioner(std::vector<float>& inOutPrecon, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ)
{
	inOutPrecon.assign(mNumCells, 0.f);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (mFluidCell[index])
		{
			float Axi = 0.0f;
			float Axj = 0.0f;
			float Axk = 0.0f;

			float Ayi = 0.0f;
			float Ayj = 0.0f;
			float Ayk = 0.0f;

			float Azi = 0.0f;
			float Azj = 0.0f;
			float Azk = 0.0f;

			float preconi = 0.0f;
			float preconj = 0.0f;
			float preconk = 0.0f;

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
				inOutPrecon[index] = 1 / sqrt(newPrecon);
			}
					
		}
	}
}

void MACGrid::ApplyPreconditioner(std::vector<float>& outResult, const std::vector<float>& inResidual, const std::vector<float>& inPrecon, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ)
{
	std::vector<float> intermediate;
	intermediate.assign(mNumCells, 0.f);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (mFluidCell[index])
		{

			float Axi = 0.0f;
			float Ayj = 0.0f;
			float Azk = 0.0f;

			float preconi = 0.0f;
			float preconj = 0.0f;
			float preconk = 0.0f;

			float intermediatei = 0.0f;
			float intermediatej = 0.0f;
			float intermediatek = 0.0f;

			if (x > 0)
			{
				int neighbourLeft = GetIndexFromXYZ(x - 1, y, z);

				Axi = inX[neighbourLeft];
				preconi = inPrecon[neighbourLeft];
				intermediatei = intermediate[neighbourLeft];
			}

			if (y > 0)
			{
				int neighbourBottom = GetIndexFromXYZ(x, y - 1, z);

				Ayj = inY[neighbourBottom];
				preconj = inPrecon[neighbourBottom];
				intermediatej = intermediate[neighbourBottom];
			}

			if (z > 0)
			{
				int neighbourBack = GetIndexFromXYZ(x, y, z - 1);

				Azk = inZ[neighbourBack];
				preconk = inPrecon[neighbourBack];
				intermediatek = intermediate[neighbourBack];
			}

			float t = inResidual[index] - Axi * preconi * intermediatei
										- Ayj * preconj * intermediatej
										- Azk * preconk * intermediatek;

			intermediate[index] = t * inPrecon[index];
		}
	}

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y, z;
		std::tie(x, y, z) = GetXYZFromIndex(index);

		if (mFluidCell[index])
		{
			float zi = 0.0f;
			float zj = 0.0f;
			float zk = 0.0f;

			if (x < mNumCellWidth - 1)
			{
				int neighbourRight = GetIndexFromXYZ(x + 1, y, z);

				zi = outResult[neighbourRight];
			}

			if (y < mNumCellHeight - 1)
			{
				int neighbourTop = GetIndexFromXYZ(x, y + 1, z);

				zj = outResult[neighbourTop];
			}

			if (z < mNumCellLength - 1)
			{
				int neighbourFront= GetIndexFromXYZ(x, y, z + 1);

				zk = outResult[neighbourFront];
			}

			float t = intermediate[index] - inX[index] * inPrecon[index] * zi
											- inY[index] * inPrecon[index] * zj
											- inZ[index] * inPrecon[index] * zk;

			outResult[index] = t * inPrecon[index];
		}
	}
}

std::tuple<int, int, int> MACGrid::GetXYZFromIndex(int index)
{
	int x = 0;
	int y = 0;
	int z = 0;

	z = index % mNumCellLength;
	y = static_cast<int>(round((index - z) / mNumCellLength)) % mNumCellHeight;
	x = static_cast<int>(round(index - y * mNumCellHeight - z) / (mNumCellLength * mNumCellHeight)) % mNumCellWidth;

	return std::tuple<int, int, int>(x, y, z);
}

int MACGrid::GetIndexFromXYZ(int X, int Y, int Z)
{
	return  Z + Y * mNumCellLength + X * mNumCellHeight * mNumCellLength;
}

int MACGrid::GetNumNeighbourOfFluidCells(int index)
{
	int x, y, z;
	std::tie(x, y, z) = GetXYZFromIndex(index);

	int back = z > 0 ? GetIndexFromXYZ(x, y, z - 1) : -1;
	int front = z < mNumCellLength - 1 ? GetIndexFromXYZ(x, y, z + 1) : -1;

	int numFluidNeighbours = 0;

	if (x > 0)
	{
		int left = GetIndexFromXYZ(x - 1, y, z);

		if (mFluidCell[left])
		{
			++numFluidNeighbours;
		}
		
	}
	if (x < mNumCellWidth - 1)
	{
		int right = GetIndexFromXYZ(x + 1, y, z);

		if (mFluidCell[right])
		{
			++numFluidNeighbours;
		}
	}

	if (y > 0)
	{
		int bottom = GetIndexFromXYZ(x, y - 1, z);

		if (mFluidCell[bottom])
		{
			++numFluidNeighbours;
		}
	}
	if (y < mNumCellHeight - 1)
	{
		int top = GetIndexFromXYZ(x, y + 1, z);

		if (mFluidCell[top])
		{
			++numFluidNeighbours;
		}
	}
	
	if (z > 0)
	{
		int back = GetIndexFromXYZ(x, y, z - 1);

		if (mFluidCell[back])
		{
			++numFluidNeighbours;
		}
	}
	if (z < mNumCellLength - 1)
	{
		int front = GetIndexFromXYZ(x, y, z + 1);

		if (mFluidCell[front])
		{
			++numFluidNeighbours;
		}
	}

	return numFluidNeighbours;

}