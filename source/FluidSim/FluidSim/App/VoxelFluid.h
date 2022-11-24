#pragma once

#include <glm/vec3.hpp>
#include <vector>

#include "Fluid.h"
#include "Domain.h"

class VoxelFluid
{
	public:
		VoxelFluid() = default;
		~VoxelFluid() = default;

		VoxelFluid(const Fluid& inFluid, float inVoxelSize);

		void UpdateVoxelStates(const Fluid& inFluid);

		const std::vector<glm::vec3>& GetVoxelCenters() const { return mVoxelCenters; }
		const glm::vec3& GetVoxelCenter(int index) const { return mVoxelCenters[index]; }
		
		const std::vector<bool>& GetVoxelStates() const { return mFluidVoxel; }
		const bool GetVoxelState(int index) const { return mFluidVoxel[index]; }
		void SetVoxelState(int index, bool inState) { mFluidVoxel[index] = inState; }

	private:

		void InitializeFromDomain(const Domain& inDomain);
		
		float mVoxelSize;

		int mNumVoxelsLength;
		int mNumVoxelsHeight;
		int mNumVoxelsWidth;

		std::vector<glm::vec3> mVoxelCenters;
		std::vector<bool> mFluidVoxel;

};