#include "VoxelFluid2D.h"

#include <iostream>

VoxelFluid2D::VoxelFluid2D(const Fluid2D& inFluid)
{
	mVoxelSize = inFluid.GetDomain().GetWidth() / inFluid.GetMACGridResolution();

	InitializeFromDomain(inFluid.GetDomain());
}

void VoxelFluid2D::InitializeFromDomain(const Domain2D& inDomain)
{
	glm::vec2 dCenter = inDomain.GetCenter();

	float dLeft = inDomain.GetLeft();
	float dBottom = inDomain.GetBottom();

	float dWidth = inDomain.GetWidth();
	float dHeight = inDomain.GetHeight();

	mNumVoxelsWidth = floor(dWidth / mVoxelSize);
	mNumVoxelsHeight = floor(dHeight / mVoxelSize);

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

void VoxelFluid2D::UpdateVoxelStates(const Fluid2D& inFluid)
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