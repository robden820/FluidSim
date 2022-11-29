#include "MACGrid.h"

#include <iostream>

#define TOLERANCE 0.000001f
#define PRECON_TUNER 0.97f
#define PRECON_SAFETY 0.25f

MACGrid::MACGrid(const Domain& inDomain, int inGridResolution)
{
	InitializeFromDomain(inDomain, inGridResolution);
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

	int numCells = mNumCellLength * mNumCellWidth * mNumCellHeight;

	mCellSize = dLength / inGridResolution;
	mInvCellSize = 1 / mCellSize;

	float halfCell = mCellSize * 0.5f;

	mGridCells.reserve(numCells);
	mCellCenters.reserve(numCells);

	for (int x = 0; x < mNumCellWidth; x++)
	{
		float centerX = dLeft + (x * mCellSize) + halfCell;

		for (int y = 0; y < mNumCellHeight; y++)
		{
			float centerY = dBottom + (y * mCellSize) + halfCell;

			for (int z = 0; z < mNumCellLength; z++)
			{
				float centerZ = dBack + (z * mCellSize) + halfCell;

				MACGridCell gridCell;

				mGridCells.push_back(gridCell);
				mCellCenters.push_back(glm::vec3(centerX, centerY, centerZ));
			}
		}
	}
}

void MACGrid::Update(float deltaTime)
{
	UpdateCellVelocity(deltaTime);

	UpdateCellPressure(deltaTime, 200);
}

void MACGrid::UpdateCellVelocity(float deltaTime)
{
	int index = 0;
	
	float scale = deltaTime * mInvCellSize; // TO DO: water density

	for (int x = 0; x < mNumCellWidth; x++)
	{
		for (int y = 0; y < mNumCellHeight; y++)
		{
			for (int z = 0; z < mNumCellLength; z++)
			{
				if (x < mNumCellWidth - 1)
				{
					int neighbourRight = z + y * mNumCellLength + (x + 1) * mNumCellHeight * mNumCellLength;
					
					float velocity = mGridCells[index].GetFaceVelocity(MACGridCell::eRIGHT);
					
					velocity -= scale * (mGridCells[neighbourRight].GetPressure() - mGridCells[index].GetPressure());

					mGridCells[index].SetFaceVelocity(MACGridCell::eRIGHT, velocity);
				}

				if (y < mNumCellHeight - 1)
				{
					int neighbourTop = z + (y + 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					float velocity = mGridCells[index].GetFaceVelocity(MACGridCell::eTOP);

					velocity -= scale * (mGridCells[neighbourTop].GetPressure() - mGridCells[index].GetPressure());

					mGridCells[index].SetFaceVelocity(MACGridCell::eTOP, velocity);
				}

				if (z < mNumCellLength - 1)
				{
					int neighbourFront = (z + 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					float velocity = mGridCells[index].GetFaceVelocity(MACGridCell::eFRONT);

					velocity -= scale * (mGridCells[neighbourFront].GetPressure() - mGridCells[index].GetPressure());

					mGridCells[index].SetFaceVelocity(MACGridCell::eFRONT, velocity);
				}

				++index;
			}
		}
	}

}

void MACGrid::CalculateCellDivergence(float deltaTime)
{
	int index = 0;

	float scale = mInvCellSize;

	mCellDivergence.assign(mNumCellWidth * mNumCellHeight * mNumCellLength, 0.f);

	for (int x = 0; x < mNumCellWidth; x++)
	{
		for (int y = 0; y < mNumCellHeight; y++)
		{
			for (int z = 0; z < mNumCellLength; z++)
			{
				float divergence = 0.f;

				if (x < mNumCellWidth - 1)
				{
					int neighbourRight = z + y * mNumCellLength + (x + 1) * mNumCellHeight * mNumCellLength;

					divergence += mGridCells[neighbourRight].GetFaceVelocity(MACGridCell::eRIGHT) - mGridCells[index].GetFaceVelocity(MACGridCell::eRIGHT);
				}

				if (y < mNumCellHeight - 1)
				{
					int neighbourTop = z + (y + 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					divergence += mGridCells[neighbourTop].GetFaceVelocity(MACGridCell::eTOP) - mGridCells[index].GetFaceVelocity(MACGridCell::eTOP);
				}

				if (z < mNumCellLength - 1)
				{
					int neighbourFront = (z + 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					divergence += mGridCells[neighbourFront].GetFaceVelocity(MACGridCell::eFRONT) - mGridCells[index].GetFaceVelocity(MACGridCell::eFRONT);
				}

				divergence *= -scale;

				mCellDivergence[index] = divergence;

				++index;
			}
		}
	}
}

void MACGrid::UpdateCellPressure(float deltaTime, int maxIterations)
{
	int numCells = mNumCellHeight * mNumCellWidth * mNumCellLength;

	std::vector<float> Adiagonal;
	std::vector<float> Ax(numCells);
	std::vector<float> Ay(numCells);
	std::vector<float> Az(numCells);

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

	for (iteration = 0; iteration < maxIterations; iteration++)
	{
		ApplyA(z, search, Adiagonal, Ax, Ay, Az);

		float phi = 0.0f;
		for (int i = 0; i < numCells; i++)
		{
			phi += z[i] * search[i];
		}

		float alpha = theta / phi;

		for (int i = 0; i < numCells; i++)
		{
			newPressure[i] += alpha * search[i];
			residuals[i] -= alpha * z[i];
		}

		float maxResidual = 0.0f;

		for (int index = 0; index < numCells; index++)
		{
			if (residuals[index] > maxResidual)
			{
				maxResidual = residuals[index];
			}
		}

		if (maxResidual < TOLERANCE)
		{
			converged = true;
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

		if (converged)
		{
			break;
		}
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
		mGridCells[i].SetPressure(newPressure[i]);

		/*
		if (mGridCells[i].GetCellType() == MACGridCell::eFLUID)
		{
			mGridCells[i].SetPressure(newPressure[i]);
		}
		else
		{
			mGridCells[i].SetPressure(0.1f);
		}
		*/
	}
}

void MACGrid::InitializeLinearSystem(float deltaTime, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ)
{
	float scale = deltaTime * mInvCellSize * mInvCellSize; // TO DO include fluid density.

	int index = 0;

	for (int x = 0; x < mNumCellWidth; x++)
	{
		for (int y = 0; y < mNumCellHeight; y++)
		{
			for (int z = 0; z < mNumCellLength; z++)
			{
				if (x > 0)
				{
					int neighbourLeft = z + y * mNumCellLength + (x - 1) * mNumCellHeight * mNumCellLength;

					if (mGridCells[neighbourLeft].GetCellType() == MACGridCell::eFLUID)
					{
						inDiag[index] += scale;
					}
				}

				if (x < mNumCellWidth - 1)
				{
					int neighbourRight = z + y * mNumCellLength + (x + 1) * mNumCellHeight * mNumCellLength;

					if (mGridCells[neighbourRight].GetCellType() == MACGridCell::eFLUID)
					{
						inDiag[index] += scale;
						inX[index] -= scale;
					}
					else if (mGridCells[neighbourRight].GetCellType() == MACGridCell::eEMPTY)
					{
						inDiag[index] += scale;
					}
				}

				if (y > 0)
				{
					int neighbourBottom = z + (y - 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					if (mGridCells[neighbourBottom].GetCellType() == MACGridCell::eFLUID)
					{
						inDiag[index] += scale;
					}
				}

				if (y < mNumCellHeight - 1)
				{
					int neighbourTop = z + (y + 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					if (mGridCells[neighbourTop].GetCellType() == MACGridCell::eFLUID)
					{
						inDiag[index] += scale;
						inY[index] -= scale;
					}
					else if (mGridCells[neighbourTop].GetCellType() == MACGridCell::eEMPTY)
					{
						inDiag[index] += scale;
					}
				}

				if (z > 0)
				{
					int neighbourBack = (z - 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					if (mGridCells[neighbourBack].GetCellType() == MACGridCell::eFLUID)
					{
						inDiag[index] += scale;
					}
				}

				if (z < mNumCellWidth - 1)
				{
					int neighbourFront = (z + 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					if (mGridCells[neighbourFront].GetCellType() == MACGridCell::eFLUID)
					{
						inDiag[index] += scale;
						inZ[index] -= scale;
					}
					else if (mGridCells[neighbourFront].GetCellType() == MACGridCell::eEMPTY)
					{
						inDiag[index] += scale;
					}
				}

				++index;
			}
		}
	}
}

void MACGrid::ApplyA(std::vector<float>& outResult, std::vector<float>& inVec, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ)
{
	int numCells = mNumCellHeight * mNumCellWidth * mNumCellLength;
	float scale = 0.002f * mInvCellSize * mInvCellSize; // TO DO;

	int index = 0;

	for (int x = 0; x < mNumCellWidth; x++)
	{
		for (int y = 0; y < mNumCellHeight; y++)
		{
			for (int z = 0; z < mNumCellLength; z++)
			{

				float value = 0.f;

				if (x > 0)
				{
					int neighbourLeft = z + y * mNumCellLength + (x - 1) * mNumCellHeight * mNumCellLength;
					value += inVec[neighbourLeft];
				}

				if (x < mNumCellWidth - 1)
				{
					int neighbourRight = z + y * mNumCellLength + (x + 1) * mNumCellHeight * mNumCellLength;
					value += inVec[neighbourRight];
				}

				if (y > 0)
				{
					int neighbourBottom = z + (y - 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;
					value += inVec[neighbourBottom];
				}

				if (y < mNumCellHeight - 1)
				{
					int neighbourTop = z + (y + 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;
					value += inVec[neighbourTop];
				}

				if (z > 0)
				{
					int neighbourBack = (z - 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;
					value += inVec[neighbourBack];
				}

				if (z < mNumCellLength - 1)
				{
					int neighbourFront = (z + 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;
					value += inVec[neighbourFront];
				}

				value += inDiag[index] * scale * inVec[index];

				outResult[index] = value;

				++index;
			}
		}
	}
}

void MACGrid::CalculatePreconditioner(std::vector<float>& inOutPrecon, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ)
{
	inOutPrecon.assign(mNumCellHeight * mNumCellLength * mNumCellWidth, 0.f);

	int index = 0;

	for (int x = 0; x < mNumCellWidth; x++)
	{
		for (int y = 0; y < mNumCellHeight; y++)
		{
			for (int z = 0; z < mNumCellLength; z++)
			{
				if (mGridCells[index].GetCellType() == MACGridCell::eFLUID)
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
						int neighbourLeft = z + y * mNumCellLength + (x - 1) * mNumCellHeight * mNumCellLength;

						Axi = inX[neighbourLeft];
						Ayi = inY[neighbourLeft];
						Azi = inZ[neighbourLeft];

						preconi = inOutPrecon[neighbourLeft];
					}
					
					if (y > 0)
					{
						int neighbourBottom = z + (y - 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;

						Axj = inX[neighbourBottom];
						Ayj = inY[neighbourBottom];
						Azj = inZ[neighbourBottom];

						preconj = inOutPrecon[neighbourBottom];
					}

					if (z > 0)
					{
						int neighbourBack = (z - 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;

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
					

					inOutPrecon[index] = 1 / sqrt(newPrecon);
				}

				++index;
			}
		}
	}
}

void MACGrid::ApplyPreconditioner(std::vector<float>& outResult, const std::vector<float>& inResidual, const std::vector<float>& inPrecon, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ)
{

	std::vector<float> intermediate;
	intermediate.assign(mNumCellHeight * mNumCellLength * mNumCellWidth, 0.f);

	int index = 0;

	for (int x = 0; x < mNumCellWidth; x++)
	{
		for (int y = 0; y < mNumCellHeight; y++)
		{
			for (int z = 0; z < mNumCellLength; z++)
			{
				if (mGridCells[index].GetCellType() == MACGridCell::eFLUID)
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
						int neighbourLeft = z + y * mNumCellLength + (x - 1) * mNumCellHeight * mNumCellLength;

						Axi = inX[neighbourLeft];
						preconi = inPrecon[neighbourLeft];
						intermediatei = intermediate[neighbourLeft];

					}

					if (y > 0)
					{
						int neighbourBottom = z + (y - 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;

						Ayj = inY[neighbourBottom];
						preconj = inPrecon[neighbourBottom];
						intermediatej = intermediate[neighbourBottom];
					}

					if (z > 0)
					{
						int neighbourBack = (z - 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;

						Azk = inZ[neighbourBack];
						preconk = inPrecon[neighbourBack];
						intermediatek = intermediate[neighbourBack];
					}

					float t = inResidual[index] - Axi * preconi * intermediatei
												- Ayj * preconj * intermediatej
												- Azk * preconk * intermediatek;

					intermediate[index] = t * inPrecon[index];
				}
				
				++index;
			}
		}
	}

	index = 0;

	for (int x = 0; x < mNumCellWidth; x++)
	{
		for (int y = 0; y < mNumCellHeight; y++)
		{
			for (int z = 0; z < mNumCellLength; z++)
			{
				if (mGridCells[index].GetCellType() == MACGridCell::eFLUID)
				{
					float zi = 0.0f;
					float zj = 0.0f;
					float zk = 0.0f;

					if (x < mNumCellWidth - 1)
					{
						int neighbourRight = z + y * mNumCellLength + (x + 1) * mNumCellHeight * mNumCellLength;

						zi = outResult[neighbourRight];
					}

					if (y < mNumCellHeight - 1)
					{
						int neighbourTop = z + (y + 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;

						zj = outResult[neighbourTop];
					}

					if (z < mNumCellLength - 1)
					{
						int neighbourFront= (z + 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;

						zk = outResult[neighbourFront];
					}

					float t = intermediate[index] - inX[index] * inPrecon[index] * zi
													- inY[index] * inPrecon[index] * zj
													- inZ[index] * inPrecon[index] * zk;

					outResult[index] = t * inPrecon[index];
				}

				++index;
			}
		}
	}
}