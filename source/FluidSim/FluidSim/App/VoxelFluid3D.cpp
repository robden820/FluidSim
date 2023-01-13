#include "VoxelFluid3D.h"

#include <iostream>

VoxelFluid3D::VoxelFluid3D(const ApplicationData& inData)
{
	mVoxelSize = inData.GetGridCellSize();

	Initialize(inData);
}

void VoxelFluid3D::Initialize(const ApplicationData& inData)
{
	float dLeft = inData.GetGridLeft();
	float dBack = inData.GetGridBack();
	float dBottom = inData.GetGridBottom();

	float dLength = inData.GetGridLength();
	float dWidth = inData.GetGridWidth();
	float dHeight = inData.GetGridHeight();

	mNumVoxelsLength = inData.GetNumGridCellsLength();
	mNumVoxelsWidth = inData.GetNumGridCellsWidth();
	mNumVoxelsHeight = inData.GetNumGridCellsHeight();

	mNumVoxels = mNumVoxelsLength * mNumVoxelsWidth * mNumVoxelsHeight;

	float halfVoxel = mVoxelSize * 0.5f;

	mVoxelCenters.reserve(mNumVoxels);
	mVoxelStates.reserve(mNumVoxels);

	for (int x = 0; x < mNumVoxelsWidth; x++)
	{
		float centerX = dLeft + (x * mVoxelSize) + halfVoxel;

		for (int y = 0; y < mNumVoxelsHeight; y++)
		{
			float centerY = dBottom + (y * mVoxelSize) + halfVoxel;

			for (int z = 0; z < mNumVoxelsLength; z++)
			{
				float centerZ = dBack + (z * mVoxelSize) + halfVoxel;

				mVoxelCenters.push_back(glm::vec3(centerX, centerY, centerZ));
				mVoxelStates.push_back(VoxelState::eNONE);
			}
		}
	}
}