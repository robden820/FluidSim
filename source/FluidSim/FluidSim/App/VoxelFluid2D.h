#pragma once

#include <glm/vec2.hpp>
#include <vector>

#include "VoxelFluid.h"

#include "Fluid2D.h"
#include "Domain2D.h"

class VoxelFluid2D : public VoxelFluid
{
public:

	VoxelFluid2D() = default;
	~VoxelFluid2D() = default;

	VoxelFluid2D(const Fluid2D& inFluid);

	void UpdateVoxelStates(const Fluid2D& inFluid);

	const std::vector<glm::vec2>& GetVoxelCenters() const { return mVoxelCenters; }
	const glm::vec2& GetVoxelCenter(int index) const { return mVoxelCenters[index]; }

	const std::vector<VoxelState>& GetVoxelStates() const { return mVoxelStates; }
	const VoxelState GetVoxelState(int index) const { return mVoxelStates[index]; }
	void SetVoxelState(int index, VoxelState inState) { mVoxelStates[index] = inState; }

private:

	void InitializeFromDomain(const Domain2D& inDomain);

	int mNumVoxelsLength;
	int mNumVoxelsHeight;
	int mNumVoxelsWidth;

	std::vector<glm::vec2> mVoxelCenters;

};
