#pragma once

#include <vector>

#include "MACGridCell.h"
#include "Fluid.h"
#include "Domain.h"

#include "glm/vec3.hpp"

class MACGrid
{
public:
	MACGrid() = default;
	~MACGrid() = default;

	MACGrid(const Fluid& inFluid, float inGridResolution);

	void Update(float deltaTime);

private:

	void InitializeFromDomain(const Domain& inDomain, float inGridResolution);

	void UpdateCellDivergence(float deltaTime);
	void UpdateCellPressure(float deltaTime);
	void UpdateCellVelocity(float deltaTime);

	float mNumCellWidth;
	float mNumCellLength;
	float mNumCellHeight;

	float mCellSize;    // deltaX
	float mInvCellSize; // 1 / deltaX

	std::vector<MACGridCell> mGridCells;
	std::vector<glm::vec3> mCellCenters;
	std::vector<glm::vec3> mCellDivergence;
};