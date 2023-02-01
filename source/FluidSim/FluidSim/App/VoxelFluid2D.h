#pragma once

#include <glm/vec2.hpp>
#include <vector>

#include "VoxelFluid.h"

class VoxelFluid2D : public VoxelFluid
{
public:

	VoxelFluid2D();
	~VoxelFluid2D() = default;

	VoxelFluid2D(const ApplicationData& inData);


	const std::vector<glm::vec2>& GetVoxelCenters() const { return mVoxelCenters; }
	const glm::vec2& GetVoxelCenter(int index) const { return mVoxelCenters[index]; }

	const std::vector<VoxelState>& GetVoxelStates() const { return mVoxelStates; }
	const VoxelState GetVoxelState(int index) const { return mVoxelStates[index]; }
	void SetVoxelState(int index, VoxelState inState) { mVoxelStates[index] = inState; }

private:
	void Initialize(const ApplicationData& inData) override;

	int mNumVoxelsHeight;
	int mNumVoxelsWidth;

	std::vector<glm::vec2> mVoxelCenters;

};
