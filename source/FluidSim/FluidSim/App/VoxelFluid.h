#pragma once

#include <glm/vec3.hpp>
#include <vector>

#include "Fluid3D.h"
#include "Domain3D.h"

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

		float mVoxelSize;
		int mNumVoxels;

		std::vector<VoxelState> mVoxelStates;
};