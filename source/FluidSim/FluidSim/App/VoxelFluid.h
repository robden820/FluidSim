#pragma once

#include <vector>

#include "ApplicationData.h"

class VoxelFluid
{
	public:
		enum VoxelState
		{
			eFLUID = 0,
			eSOLID = 1,
			eNONE = 2
		};

		virtual ~VoxelFluid() = default;

	protected:
		virtual void Initialize(const ApplicationData& inData) = 0;

		float mVoxelSize;
		int mNumVoxels;

		std::vector<VoxelState> mVoxelStates;
};