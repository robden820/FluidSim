#pragma once

#include <vector>

#include "MACGridCell.h"
#include "Domain.h"

#include "glm/vec3.hpp"

class MACGrid
{
public:
	MACGrid() = default;
	~MACGrid() = default;

	MACGrid(const Domain& inDomain, int inGridResolution);

	void Update(float deltaTime);

	int GetNumCells() const { return mNumCellHeight * mNumCellWidth * mNumCellLength; }

	const MACGridCell& GetGridCell(int index) const { return mGridCells[index]; }
	const glm::vec3& GetCellCenter(int index) const { return mCellCenters[index]; }

	void SetGridCellVelocity(int index, const glm::vec3& inVelocity) { mGridCells[index].SetCellVelocity(inVelocity); }
	void SetGridCellType(int index, MACGridCell::CellType inType) { mGridCells[index].SetCellType(inType); }

private:

	void InitializeFromDomain(const Domain& inDomain, int inGridResolution);

	void CalculateCellDivergence(float deltaTime);
	void UpdateCellPressure(float deltaTime, int maxIterations);
	void UpdateCellVelocity(float deltaTime);

	void InitializeLinearSystem(float deltaTime, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ);

	void CalculatePreconditioner(std::vector<float>& inOutPrecon, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ);

	void ApplyA(std::vector<float>& outResult, std::vector<float>& inVec, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ);
	void ApplyPreconditioner(std::vector<float>& outResult, const std::vector<float>& inResidual, const std::vector<float>& inPrecon, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ);

	int mNumCellWidth;
	int mNumCellLength;
	int mNumCellHeight;

	float mCellSize;    // deltaX
	float mInvCellSize; // 1 / deltaX

	std::vector<MACGridCell> mGridCells;
	std::vector<glm::vec3> mCellCenters;
	std::vector<float> mCellDivergence;
};