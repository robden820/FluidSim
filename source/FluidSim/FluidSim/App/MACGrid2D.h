#pragma once

#include <vector>

#include "MACGrid.h"

#include "glm/glm.hpp"
#include "oneapi/tbb.h"

class MACGrid2D : public MACGrid
{
public:
	MACGrid2D() = default;
	~MACGrid2D() = default;

	MACGrid2D(const ApplicationData& inData);

	void Update(ApplicationData& inOutData);

	int GetNumCells() const { return mNumCellHeight * mNumCellWidth; }
	int GetNumCellsWidth() const { return mNumCellWidth; }
	int GetNumCellsHeight() const { return mNumCellHeight; }

	const glm::vec2& GetCellCenter(int index) const { return mCellCenters[index]; }
	int GetClosestCell(const glm::vec2& inPos);

	const float GetCellXVelocity(int index) const { return mCellXVelocities[index]; }
	const float GetCellYVelocity(int index) const { return mCellYVelocities[index]; }

	void SetCellXVelocity(int index, float inVelocity) { mCellXVelocities[index] = inVelocity; }
	void SetCellYVelocity(int index, float inVelocity) { mCellYVelocities[index] = inVelocity; }

	const float GetCellPressure(int index) const { return mCellPressures[index]; }
	void SetCellPressure(int index, float inPressure) { mCellPressures[index] = inPressure; }

	const CellType GetCellType(int index) const { return mCellType[index]; }
	const CellType GetCellTypeFromPosition(const glm::vec2& inPos);
	void SetCellType(int index, CellType inCellType) { mCellType[index] = inCellType; }

	std::tuple<int, int> GetXYFromIndex(int index);
	int GetIndexFromXY(int X, int Y);

	float GetCellSize() const { return mCellSize; }
	float GetInverseCellSize() const { return mInvCellSize; }

private:

	void InitializeGrid(const ApplicationData& inData) override;
	void InitializeCellsFromParticles(const std::vector<glm::vec2>& inParticlePositions);

	void CalculateCellDivergence(float deltaTime);

	void AdvectCellVelocity(float deltaTime);
	void UpdateCellPressure(float deltaTime, int maxIterations);
	void UpdateCellVelocity(float deltaTime);

	void InitializeLinearSystem(float deltaTime, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY);

	void CalculatePreconditioner(std::vector<float>& inOutPrecon, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY);

	void ApplyA(float deltaTime, std::vector<float>& outResult, const std::vector<float>& inVec, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY);
	void ApplyPreconditioner(std::vector<float>& outResult, const std::vector<float>& inResidual, const std::vector<float>& inPrecon, const std::vector<float>& inX, const std::vector<float>& inY);

	int mNumCellWidth;
	int mNumCellHeight;

	float dLeft;
	float dBottom;

	std::vector<glm::vec2> mCellCenters;

	std::vector<float> mCellXVelocities;
	std::vector<float> mCellYVelocities;

	// Intermediate cell velocities
	std::vector<float> mIntXVelocities;
	std::vector<float> mIntYVelocities;

	std::vector<float> mCellPressures;
};