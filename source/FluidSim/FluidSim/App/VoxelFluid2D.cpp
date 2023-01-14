#include "VoxelFluid2D.h"

#include <iostream>

VoxelFluid2D::VoxelFluid2D(const ApplicationData& inData)
{
	mVoxelSize = inData.GetGridCellSize();

	Initialize(inData);
}

void VoxelFluid2D::Initialize(const ApplicationData& inData)
{
	float dLeft = inData.GetGridLeft();
	float dBottom = inData.GetGridBottom();

	mNumVoxelsWidth = inData.GetNumGridCellsWidth();
	mNumVoxelsHeight = inData.GetNumGridCellsHeight();

	mNumVoxels = mNumVoxelsWidth * mNumVoxelsHeight;

	float halfVoxel = mVoxelSize * 0.5f;

	mVoxelCenters.reserve(mNumVoxels);
	mVoxelStates.reserve(mNumVoxels);

	for (int x = 0; x < mNumVoxelsWidth; x++)
	{
		float centerX = dLeft + (x * mVoxelSize) + halfVoxel;

		for (int y = 0; y < mNumVoxelsHeight; y++)
		{
			float centerY = dBottom + (y * mVoxelSize) + halfVoxel;

			mVoxelCenters.push_back(glm::vec2(centerX, centerY));
			mVoxelStates.push_back(VoxelState::eNONE);
		}
	}
}