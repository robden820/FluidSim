#include "VoxelFluid.h"

#include <iostream>

VoxelFluid::VoxelFluid(const Fluid& inFluid)
{
	mVoxelSize = inFluid.GetDomain().GetLength() / inFluid.GetMACGridResolution();

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

void VoxelFluid::UpdateVoxelStates(const Fluid& inFluid)
{
	for (int index = 0; index < mNumVoxels; index++)
	{
		mVoxelStates[index] = VoxelState::eNONE;

		if (inFluid.GetMACGrid().GetCellType(index) == MACGrid::CellType::eFLUID)
		{
			mVoxelStates[index] = VoxelState::eFLUID;
		}
		else if (inFluid.GetMACGrid().GetCellType(index) == MACGrid::CellType::eSOLID)
		{
			mVoxelStates[index] = VoxelState::eSOLID;
		}
	}
}