#pragma once

#include <vector>

#include "Domain.h"

#include "glm/glm.hpp"

class MACGrid
{
public:

	enum CellType
	{
		eFLUID = 0,
		eSOLID = 1,
		eAIR = 2,
		eNONE = 3
	};

	MACGrid() = default;
	~MACGrid() = default;

	MACGrid(const Domain& inDomain, const std::vector<glm::vec3>& inParticlePositions, int inGridResolution);

	void Update(float deltaTime);

	int GetNumCells() const { return mNumCellHeight * mNumCellWidth * mNumCellLength; }

	const glm::vec3& GetCellCenter(int index) const { return mCellCenters[index]; }
	int GetClosestCell(const glm::vec3& inPosition);

	const float GetCellXVelocity(int index) const { return mCellXVelocities[index]; }
	const float GetCellYVelocity(int index) const { return mCellYVelocities[index]; }
	const float GetCellZVelocity(int index) const { return mCellZVelocities[index]; }

	void SetCellXVelocity(int index, float inVelocity) { mCellXVelocities[index] = inVelocity; }
	void SetCellYVelocity(int index, float inVelocity) { mCellYVelocities[index] = inVelocity; }
	void SetCellZVelocity(int index, float inVelocity) { mCellZVelocities[index] = inVelocity; }

	void SetCellVelocity(int index, const glm::vec3& inVelocity);

	const float GetCellPressure(int index) const { return mCellPressures[index]; }
	void SetCellPressure(int index, float inPressure) { mCellPressures[index]= inPressure; }

	const CellType GetCellType(int index) const { return mCellType[index]; }
	void SetCellType(int index, CellType inCellType) { mCellType[index] = inCellType; }

	std::tuple<int, int, int> GetXYZFromIndex(int index);
	int GetIndexFromXYZ(int X, int Y, int Z);

private:

	void InitializeFromDomain(const Domain& inDomain, int inGridResolution);
	void InitializeCellsFromParticles(const std::vector<glm::vec3>& inParticlePositions);

	void CalculateCellDivergence(float deltaTime);

	void AdvectCellVelocity(float deltaTime);
	void UpdateCellPressure(float deltaTime, int maxIterations);
	void UpdateCellVelocity(float deltaTime);

	void InitializeLinearSystem(float deltaTime, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ);

	void CalculatePreconditioner(std::vector<float>& inOutPrecon, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ);

	void ApplyA(float deltaTime, std::vector<float>& outResult, const std::vector<float>& inVec, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ);
	void ApplyPreconditioner(std::vector<float>& outResult, const std::vector<float>& inResidual, const std::vector<float>& inPrecon, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ);

	int mNumCellWidth;
	int mNumCellLength;
	int mNumCellHeight;
	int mNumCells;

	float mCellSize;    // deltaX
	float mInvCellSize; // 1 / deltaX

	std::vector<glm::vec3> mCellCenters;
	std::vector<float> mCellDivergence;

	std::vector<float> mCellPressures;
	std::vector<float> mCellXVelocities;
	std::vector<float> mCellYVelocities;
	std::vector<float> mCellZVelocities;

	std::vector<float> mIntXVelocities;
	std::vector<float> mIntYVelocities;
	std::vector<float> mIntZVelocities;

	std::vector<CellType> mCellType;
};