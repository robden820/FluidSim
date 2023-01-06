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

		VoxelFluid() = default;
		~VoxelFluid() = default;

		VoxelFluid(const Fluid3D& inFluid);

		void UpdateVoxelStates(const Fluid3D& inFluid);

		const std::vector<glm::vec3>& GetVoxelCenters() const { return mVoxelCenters; }
		const glm::vec3& GetVoxelCenter(int index) const { return mVoxelCenters[index]; }
		
		const std::vector<VoxelState>& GetVoxelStates() const { return mVoxelStates; }
		const VoxelState GetVoxelState(int index) const { return mVoxelStates[index]; }
		void SetVoxelState(int index, VoxelState inState) { mVoxelStates[index] = inState; }

	private:

		void InitializeFromDomain(const Domain3D& inDomain);

		float mVoxelSize;

		int mNumVoxelsLength;
		int mNumVoxelsHeight;
		int mNumVoxelsWidth;
		int mNumVoxels;

		std::vector<glm::vec3> mVoxelCenters;
		std::vector<VoxelState> mVoxelStates;

};