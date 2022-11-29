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

	void CalculateCellDivergence(float deltaTime);
	void UpdateCellPressure(float deltaTime, int maxIterations);
	void UpdateCellVelocity(float deltaTime);

	void InitializeLinearSystem(float deltaTime, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ);

	void ApplyA(std::vector<float>& inResult, std::vector<float>& inVec, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ);
	void ApplyPreconditioner();

	float mNumCellWidth;
	float mNumCellLength;
	float mNumCellHeight;

	float mCellSize;    // deltaX
	float mInvCellSize; // 1 / deltaX

	std::vector<MACGridCell> mGridCells;
	std::vector<glm::vec3> mCellCenters;
	std::vector<glm::vec3> mCellDivergence;
};