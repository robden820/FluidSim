#pragma once

#include <glm/vec3.hpp>
#include <vector>

#include "VoxelFluid.h"

#include "Fluid3D.h"

class VoxelFluid3D : public VoxelFluid
{
public:

	VoxelFluid3D() = default;
	~VoxelFluid3D() = default;

	VoxelFluid3D(const ApplicationData& inData);

	const std::vector<glm::vec3>& GetVoxelCenters() const { return mVoxelCenters; }
	const glm::vec3& GetVoxelCenter(int index) const { return mVoxelCenters[index]; }

	const std::vector<VoxelState>& GetVoxelStates() const { return mVoxelStates; }
	const VoxelState GetVoxelState(int index) const { return mVoxelStates[index]; }
	void SetVoxelState(int index, VoxelState inState) { mVoxelStates[index] = inState; }

private:

	void Initialize(const ApplicationData& inData) override;

	int mNumVoxelsLength;
	int mNumVoxelsHeight;
	int mNumVoxelsWidth;

	std::vector<glm::vec3> mVoxelCenters;
};
