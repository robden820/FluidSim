#include "VoxelFluid.h"

#include <iostream>

VoxelFluid::VoxelFluid(const Fluid& inFluid, float inVoxelSize)
{
	mVoxelSize = inVoxelSize;

	InitializeFromDomain(inFluid.GetDomain());
}

void VoxelFluid::InitializeFromDomain(const Domain& inDomain)
{
	glm::vec3 dCenter = inDomain.GetCenter();

	float dLeft = inDomain.GetLeft();
	float dBack = inDomain.GetBack();
	float dBottom = inDomain.GetBottom();

	float dLength = inDomain.GetLength();
	float dWidth = inDomain.GetWidth();
	float dHeight = inDomain.GetHeight();

	mNumVoxelsLength = floor(dLength / mVoxelSize);
	mNumVoxelsWidth = floor(dWidth / mVoxelSize);
	mNumVoxelsHeight = floor(dHeight / mVoxelSize);

	int numVoxels = mNumVoxelsLength * mNumVoxelsWidth * mNumVoxelsHeight;

	float halfVoxel = mVoxelSize * 0.5f;

	mVoxelCenters.reserve(numVoxels);
	mFluidVoxel.reserve(numVoxels);

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
				mFluidVoxel.push_back(false);
			}
		}
	}
}

void VoxelFluid::UpdateVoxelStates(const Fluid& inFluid)
{
	int index = 0;

	for (int x = 0; x < mNumVoxelsWidth; x++)
	{
		for (int y = 0; y < mNumVoxelsHeight; y++)
		{
			for (int z = 0; z < mNumVoxelsLength; z++)
			{
				mFluidVoxel[index] = false;

				if (inFluid.GetMACGrid().GetCellType(index))
				{
					mFluidVoxel[index] = true;
				}

				++index;
			}
		}
	}
}