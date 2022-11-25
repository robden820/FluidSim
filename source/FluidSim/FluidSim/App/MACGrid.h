#pragma once

#include <vector>

#include "MACGridCell.h"
#include "Domain.h"

class MACGrid
{
public:
	MACGrid() = default;
	~MACGrid() = default;

	MACGrid(const Domain& inDomain, float inCellSize);

	void Update(float deltaTime);

private:

	void UpdateCellPressure(float deltaTime);
	void UpdateCellVelocity(float deltaTime);

	std::vector<MACGridCell> mGridCells;
};