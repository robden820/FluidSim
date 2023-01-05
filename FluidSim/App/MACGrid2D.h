#pragma once

#include <vector>

#include "Domain2D.h"

#include "glm/glm.hpp"
//#include "oneapi/tbb.h"

class MACGrid2D
{
public:

	enum CellType
	{
		eFLUID = 0,
		eSOLID = 1,
		eAIR = 2,
		eNONE = 3
	};

	MACGrid2D() = default;
	~MACGrid2D() = default;

	MACGrid2D(const Domain2D& inDomain, const std::vector<std::tuple<float, float>>& inParticlePositions, int inGridResolution);

	void Update(float deltaTime);

	int GetNumCells() const { return mNumCellHeight * mNumCellWidth; }

	const std::tuple<float, float> GetCellCenter(int index) const { return mCellCenters[index]; }
	int GetClosestCell(const float inXPos, const float inYPos);

	const float GetCellXVelocity(int index) const { return mCellXVelocities[index]; }
	const float GetCellYVelocity(int index) const { return mCellYVelocities[index]; }

	void SetCellXVelocity(int index, float inVelocity) { mCellXVelocities[index] = inVelocity; }
	void SetCellYVelocity(int index, float inVelocity) { mCellYVelocities[index] = inVelocity; }

	const float GetCellPressure(int index) const { return mCellPressures[index]; }
	void SetCellPressure(int index, float inPressure) { mCellPressures[index] = inPressure; }

	const CellType GetCellType(int index) const { return mCellType[index]; }
	const CellType GetCellTypeFromPosition(const float inXPos, const float inYPos);
	void SetCellType(int index, CellType inCellType) { mCellType[index] = inCellType; }

	std::tuple<int, int> GetXYFromIndex(int index);
	int GetIndexFromXY(int X, int Y);

	float GetCellSize() { return mCellSize; }

private:

	void InitializeFromDomain(const Domain2D& inDomain, int inGridResolution);
	void InitializeCellsFromParticles(const std::vector<std::tuple<float, float>>& inParticlePositions);

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

	int mNumCells;

	float mCellSize;    // deltaX
	float mInvCellSize; // 1 / deltaX

	float mDensity;

	std::vector<std::tuple<float, float>> mCellCenters;
	std::vector<float> mCellDivergence;

	std::vector<float> mCellPressures;
	std::vector<float> mCellXVelocities;
	std::vector<float> mCellYVelocities;

	// Intermediate cell velocities
	std::vector<float> mIntXVelocities;
	std::vector<float> mIntYVelocities;

	std::vector<CellType> mCellType;
};