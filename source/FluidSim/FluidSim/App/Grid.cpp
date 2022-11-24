#include "Grid.h"

Grid::Grid(float inLength, float inWidth, float inHeight, float inCellSize)
{
	mLength = inLength;
	mWidth = inWidth;
	mHeight = inHeight;

	mCellSize = inCellSize;

	InitializeGridNodes();
}


Grid::Grid(const Domain& inDomain, float inCellSize)
{
	mLength = inDomain.GetLength();
	mWidth = inDomain.GetWidth();
	mHeight = inDomain.GetHeight();

	mCellSize = inCellSize;

	InitializeGridNodes();
}

void Grid::InitializeGridNodes()
{
	int numGridNodes = floor((mLength * mWidth * mHeight) / mCellSize);

	mGridNodes.reserve(numGridNodes);
	mGridNodesTemp.reserve(numGridNodes);

	for (int x = 0; x < mWidth; x++)
	{
		for (int y = 0; y < mHeight; y++)
		{
			for (int z = 0; z < mLength; z++)
			{
				GridNode node(glm::vec3(x * mCellSize, y * mCellSize, z * mCellSize));
				mGridNodes.push_back(node);
			}
		}
	}

	mGridNodesTemp = mGridNodes;
}

void Grid::StepGrid(float deltaTime)
{
	mGridNodesTemp = mGridNodes;

	for (int x = 0; x < mWidth; x++)
	{
		for (int y = 0; y < mHeight; y++)
		{
			for (int z = 0; z < mLength; z++)
			{
				int index = z + y * mLength + x * mHeight;

				int numNeighbours = 0;

				glm::vec3 velocity(0.0f, 0.0f, 0.0f);

				if (x != 0)
				{
					int neighbour = z + y * mLength + (x - 1) * mHeight;
					velocity += mGridNodes[neighbour].GetVelocity();
					numNeighbours++;
				}
				if (x != mWidth - 1)
				{
					int neighbour = z + y * mLength + (x + 1) * mHeight;
					velocity += mGridNodes[neighbour].GetVelocity();
					numNeighbours++;
				}

				if (y != 0)
				{
					int neighbour = z + (y - 1) * mLength + x * mHeight;
					velocity += mGridNodes[neighbour].GetVelocity();
					numNeighbours++;
				}
				if (y != mHeight - 1)
				{
					int neighbour = z + (y + 1) * mLength + x * mHeight;
					velocity += mGridNodes[neighbour].GetVelocity();
					numNeighbours++;
				}

				if (z != 0)
				{
					int neighbour = (z - 1) + y * mLength + x * mHeight;
					velocity += mGridNodes[neighbour].GetVelocity();
					numNeighbours++;
				}
				if (z != mLength - 1)
				{
					int neighbour = (z + 1) + y * mLength + x * mHeight;
					velocity += mGridNodes[neighbour].GetVelocity();
					numNeighbours++;
				}

				velocity *= 1 / numNeighbours;
				//velocity += glm::vec3(0.f, -9.8f, 0.f);

				mGridNodesTemp[index].SetVelocity(velocity);
			}
		}
	}

	mGridNodes = mGridNodesTemp;
}