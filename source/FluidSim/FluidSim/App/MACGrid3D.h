#pragma once

#include <vector>

#include "MACGrid.h"

#include "glm/glm.hpp"
#include "onetbb/oneapi/tbb.h"

class MACGrid3D : public MACGrid
{
public:

	MACGrid3D() = default;
	~MACGrid3D() = default;

	MACGrid3D(const ApplicationData& inData);

	void Update(ApplicationData& inOutData) override;

	int GetNumCells() const override { return mNumCellWidth * mNumCellHeight * mNumCellLength; }

	const glm::dvec3& GetCellCenter(int index) const { return mCellCenters[index]; }
	int GetClosestCell(const glm::dvec3& inPosition);

	const double GetCellXVelocity(int index) const { return mCellXVelocities[index]; }
	const double GetCellYVelocity(int index) const { return mCellYVelocities[index]; }
	const double GetCellZVelocity(int index) const { return mCellZVelocities[index]; }

	void SetCellXVelocity(int index, double inVelocity) { mCellXVelocities[index] = inVelocity; }
	void SetCellYVelocity(int index, double inVelocity) { mCellYVelocities[index] = inVelocity; }
	void SetCellZVelocity(int index, double inVelocity) { mCellZVelocities[index] = inVelocity; }

	const CellType GetCellTypeFromPosition(const glm::dvec3& inPosition);

	std::tuple<int, int, int> GetXYZFromIndex(int index);
	int GetIndexFromXYZ(int X, int Y, int Z);

private:

	void InitializeGrid(const ApplicationData& inData) override;
	void InitializeCellsFromParticles(const std::vector<glm::dvec3>& inParticlePositions);

	void CalculateCellDivergence() override;

	void AdvectCellVelocity(double deltaTime) override;
	void UpdateCellPressure(double deltaTime, int maxIterations) override;
	void UpdateCellVelocity(double deltaTime) override;

	void InitializeLinearSystem(double deltaTime, std::vector<double>& inDiag, std::vector<double>& inX, std::vector<double>& inY, std::vector<double>& inZ);

	void CalculatePreconditioner(std::vector<double>& inOutPrecon, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY, const std::vector<double>& inZ);

	void ApplyA(double deltaTime, std::vector<double>& outResult, const std::vector<double>& inVec, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY, const std::vector<double>& inZ);
	void ApplyPreconditioner(std::vector<double>& outResult, const std::vector<double>& inResidual, const std::vector<double>& inPrecon, const std::vector<double>& inX, const std::vector<double>& inY, const std::vector<double>& inZ);

	int mNumCellWidth;
	int mNumCellLength;
	int mNumCellHeight;

	double dLeft;
	double dBottom;
	double dBack;

	std::vector<glm::dvec3> mCellCenters;

	std::vector<double> mCellXVelocities;
	std::vector<double> mCellYVelocities;
	std::vector<double> mCellZVelocities;

	// Intermediate cell velocities for calculations.
	std::vector<double> mIntXVelocities;
	std::vector<double> mIntYVelocities;
	std::vector<double> mIntZVelocities;

};
