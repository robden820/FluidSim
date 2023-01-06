#pragma once

#include <vector>

#include "Domain3D.h"
#include "MACGrid.h"

#include "glm/glm.hpp"
#include "onetbb/oneapi/tbb.h"

class MACGrid3D : public MACGrid
{
public:

	MACGrid3D() = default;
	~MACGrid3D() = default;

	MACGrid3D(const Domain3D& inDomain, const std::vector<glm::vec3>& inParticlePositions, int inGridResolution);

	void Update(float deltaTime) override;

	int GetNumCells() const override { return mNumCellWidth * mNumCellHeight * mNumCellLength; }

	const glm::vec3& GetCellCenter(int index) const { return mCellCenters[index]; }
	int GetClosestCell(const glm::vec3& inPosition);

	const float GetCellXVelocity(int index) const { return mCellXVelocities[index]; }
	const float GetCellYVelocity(int index) const { return mCellYVelocities[index]; }
	const float GetCellZVelocity(int index) const { return mCellZVelocities[index]; }

	void SetCellXVelocity(int index, float inVelocity) { mCellXVelocities[index] = inVelocity; }
	void SetCellYVelocity(int index, float inVelocity) { mCellYVelocities[index] = inVelocity; }
	void SetCellZVelocity(int index, float inVelocity) { mCellZVelocities[index] = inVelocity; }

	const CellType GetCellTypeFromPosition(const glm::vec3& inPosition);

	std::tuple<int, int, int> GetXYZFromIndex(int index);
	int GetIndexFromXYZ(int X, int Y, int Z);

private:

	void InitializeFromDomain(const Domain3D& inDomain, int inGridResolution);
	void InitializeCellsFromParticles(const std::vector<glm::vec3>& inParticlePositions);

	void CalculateCellDivergence(float deltaTime) override;

	void AdvectCellVelocity(float deltaTime) override;
	void UpdateCellPressure(float deltaTime, int maxIterations) override;
	void UpdateCellVelocity(float deltaTime) override;

	void InitializeLinearSystem(float deltaTime, std::vector<float>& inDiag, std::vector<float>& inX, std::vector<float>& inY, std::vector<float>& inZ);

	void CalculatePreconditioner(std::vector<float>& inOutPrecon, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ);

	void ApplyA(float deltaTime, std::vector<float>& outResult, const std::vector<float>& inVec, const std::vector<float>& inDiag, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ);
	void ApplyPreconditioner(std::vector<float>& outResult, const std::vector<float>& inResidual, const std::vector<float>& inPrecon, const std::vector<float>& inX, const std::vector<float>& inY, const std::vector<float>& inZ);

	int mNumCellWidth;
	int mNumCellLength;
	int mNumCellHeight;

	float dLeft;
	float dBottom;
	float dBack;

	std::vector<glm::vec3> mCellCenters;

	std::vector<float> mCellXVelocities;
	std::vector<float> mCellYVelocities;
	std::vector<float> mCellZVelocities;

	// Intermediate cell velocities for calculations.
	std::vector<float> mIntXVelocities;
	std::vector<float> mIntYVelocities;
	std::vector<float> mIntZVelocities;

};
