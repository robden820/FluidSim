# include "MACGrid.h"

MACGrid::MACGrid(const Fluid& inFluid, float inGridResolution)
{
	InitializeFromDomain(inFluid.GetDomain(), inGridResolution);
}

void MACGrid::InitializeFromDomain(const Domain& inDomain, float inGridResolution)
{
	glm::vec3 dCenter = inDomain.GetCenter();

	float dLeft = inDomain.GetLeft();
	float dBack = inDomain.GetBack();
	float dBottom = inDomain.GetBottom();

	float dLength = inDomain.GetLength();
	float dWidth = inDomain.GetWidth();
	float dHeight = inDomain.GetHeight();

	mNumCellLength = floor(dLength / inGridResolution);
	mNumCellWidth = floor(dWidth / inGridResolution);
	mNumCellHeight = float(dHeight / inGridResolution);

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

	UpdateCellPressure(deltaTime);
}

void MACGrid::UpdateCellVelocity(float deltaTime)
{
	int index = 0;
	
	float scale = deltaTime * mInvCellSize; // TO DO: multiply by 1 / fluid density.

	for (int x = 0; x < mNumCellWidth; x++)
	{
		for (int y = 0; y < mNumCellHeight; y++)
		{
			for (int z = 0; z < mNumCellLength; z++)
			{
				if (x < mNumCellWidth - 1)
				{
					int neighbourRight = z + y * mNumCellLength + (x + 1) * mNumCellHeight * mNumCellLength;
					
					glm::vec3 velocity = mGridCells[index].GetVelocity(MACGridCell::eRIGHT);
					
					velocity -= scale * (mGridCells[neighbourRight].GetPressure() - mGridCells[index].GetPressure());

					mGridCells[index].SetVelocity(MACGridCell::eRIGHT, velocity);
				}

				if (y < mNumCellHeight - 1)
				{
					int neighbourTop = z + (y + 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					glm::vec3 velocity = mGridCells[index].GetVelocity(MACGridCell::eTOP);

					velocity -= scale * (mGridCells[neighbourTop].GetPressure() - mGridCells[index].GetPressure());

					mGridCells[index].SetVelocity(MACGridCell::eTOP, velocity);
				}

				if (z < mNumCellLength - 1)
				{
					int neighbourFront = (z + 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					glm::vec3 velocity = mGridCells[index].GetVelocity(MACGridCell::eFRONT);

					velocity -= scale * (mGridCells[neighbourFront].GetPressure() - mGridCells[index].GetPressure());

					mGridCells[index].SetVelocity(MACGridCell::eFRONT, velocity);
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

	for (int x = 0; x < mNumCellWidth; x++)
	{
		for (int y = 0; y < mNumCellHeight; y++)
		{
			for (int z = 0; z < mNumCellLength; z++)
			{
				glm::vec3 divergence(0.f);

				if (x < mNumCellWidth - 1)
				{
					int neighbourRight = z + y * mNumCellLength + (x + 1) * mNumCellHeight * mNumCellLength;

					divergence += mGridCells[neighbourRight].GetVelocity(MACGridCell::eRIGHT) - mGridCells[index].GetVelocity(MACGridCell::eRIGHT);
				}

				if (y < mNumCellHeight - 1)
				{
					int neighbourTop = z + (y + 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					divergence += mGridCells[neighbourTop].GetVelocity(MACGridCell::eTOP) - mGridCells[index].GetVelocity(MACGridCell::eTOP);
				}

				if (z < mNumCellLength - 1)
				{
					int neighbourFront = (z + 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					divergence += mGridCells[neighbourFront].GetVelocity(MACGridCell::eFRONT) - mGridCells[index].GetVelocity(MACGridCell::eFRONT);
				}

				divergence *= -scale;

				mCellDivergence[index] = divergence;

				++index;
			}
		}
	}
}

void MACGrid::UpdateCellPressure(float deltaTime)
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
					else if (mGridCells[neighbourRight].GetCellType() == MACGridCell::eNONE)
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
					else if (mGridCells[neighbourTop].GetCellType() == MACGridCell::eNONE)
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
					else if (mGridCells[neighbourFront].GetCellType() == MACGridCell::eNONE)
					{
						inDiag[index] += scale;
					}
				}

				++index;
			}
		}
	}
}