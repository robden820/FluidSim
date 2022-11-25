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

	float cellSize = dLength / inGridResolution;
	float halfCell = cellSize * 0.5f;

	mGridCells.reserve(numCells);
	mCellCenters.reserve(numCells);

	for (int x = 0; x < mNumCellWidth; x++)
	{
		float centerX = dLeft + (x * cellSize) + halfCell;

		for (int y = 0; y < mNumCellHeight; y++)
		{
			float centerY = dBottom + (y * cellSize) + halfCell;

			for (int z = 0; z < mNumCellLength; z++)
			{
				float centerZ = dBack + (z * cellSize) + halfCell;

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
	for (int x = 0; x < mNumCellWidth; x++)
	{
		for (int y = 0; y < mNumCellHeight; y++)
		{
			for (int z = 0; z < mNumCellLength; z++)
			{
				int index = z + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;

				if (x < mNumCellWidth - 1)
				{
					int neighbour = z + y * mNumCellLength + (x + 1) * mNumCellHeight * mNumCellLength;
					
					glm::vec3 velocity = mGridCells[index].GetVelocity(MACGridCell::eRIGHT);
					
					velocity -= deltaTime * (mGridCells[neighbour].GetPressure() - mGridCells[index].GetPressure());

					mGridCells[index].SetVelocity(MACGridCell::eRIGHT, velocity);
				}

				if (y < mNumCellHeight - 1)
				{
					int neighbour = z + (y + 1) * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					glm::vec3 velocity = mGridCells[index].GetVelocity(MACGridCell::eTOP);

					velocity -= deltaTime * (mGridCells[neighbour].GetPressure() - mGridCells[index].GetPressure());

					mGridCells[index].SetVelocity(MACGridCell::eTOP, velocity);
				}

				if (z < mNumCellLength - 1)
				{
					int neighbour = (z + 1) + y * mNumCellLength + x * mNumCellHeight * mNumCellLength;

					glm::vec3 velocity = mGridCells[index].GetVelocity(MACGridCell::eFRONT);

					velocity -= deltaTime * (mGridCells[neighbour].GetPressure() - mGridCells[index].GetPressure());

					mGridCells[index].SetVelocity(MACGridCell::eFRONT, velocity);
				}
			}
		}
	}

}

void MACGrid::UpdateCellPressure(float deltaTime)
{

}