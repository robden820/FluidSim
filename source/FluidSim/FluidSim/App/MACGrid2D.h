#pragma once

#include <vector>

#include "MACGrid.h"

#include "glm/glm.hpp"
#include "oneapi/tbb.h"

#include <Eigen/Sparse>

class MACGrid2D : public MACGrid
{
public:
	MACGrid2D();
	MACGrid2D(const ApplicationData& inData);
	~MACGrid2D() = default;

	void Update(ApplicationData& inOutData);

	void Advect(ApplicationData& inOutData);
	void ApplyForces(float deltaTime);
	void Project(ApplicationData& inOutData);

	void ExtrapolateVelocityField(bool extrapolateIntVelocities);

	int GetNumCells() const { return mNumCellHeight * mNumCellWidth; }
	int GetNumCellsWidth() const { return mNumCellWidth; }
	int GetNumCellsHeight() const { return mNumCellHeight; }

	const glm::vec2& GetCellCenter(int index) const { return mCellCenters[index]; }
	int GetClosestCell(const glm::vec2& inPos) const;

	const double GetCellXVelocity(int index) const { return mCellXVelocities[index]; }
	const double GetCellYVelocity(int index) const { return mCellYVelocities[index]; }

	const double GetPrevCellXVelocity(int index) const { return mCellXVelocitiesPrev[index]; }
	const double GetPrevCellYVelocity(int index) const { return mCellXVelocitiesPrev[index]; }

	void SetCellXVelocity(int index, double inVelocity) { mCellXVelocities[index] = inVelocity; }
	void SetCellYVelocity(int index, double inVelocity) { mCellYVelocities[index] = inVelocity; }

	void SetIntXVelocity(int index, double inVelocity) { mIntXVelocities[index] = inVelocity; }
	void SetIntYVelocity(int index, double inVelocity) { mIntYVelocities[index] = inVelocity; }

	const double GetCellPressure(int index) const { return mCellPressures[index]; }
	void SetCellPressure(int index, float inPressure) { mCellPressures[index] = inPressure; }

	const CellType GetCellType(int index) const { return mCellType[index]; }
	const std::vector<CellType>& GetCellTypes() const { return mCellType; }
	const CellType GetCellTypeFromPosition (const glm::vec2& inPos) const;
	void SetCellType(int index, CellType inCellType) { mCellType[index] = inCellType; }

	std::tuple<int, int> GetXYFromIndex(int index) const;
	int GetIndexFromXY(int X, int Y) const;

	float GetCellSize() const { return mCellSize; }
	float GetInverseCellSize() const { return mInvCellSize; }

	void UpdateCellTypesFromParticles(const std::vector<glm::vec2>& inParticlePositions);

	void UpdateApplicationData(ApplicationData& inOutData);

private:

	void InitializeGrid(const ApplicationData& inData) override;
	void InitializeGridPressure();

	void CalculateCellDivergence();

	void UpdateCellPressure(float deltaTime, int maxIterations);
	void UpdateCellPressureSpare(float deltaTime, int maxIterations);

	void AdvectCellVelocity(float deltaTime);
	void UpdateCellVelocity(float deltaTime);

	void InitializeLinearSystem(float deltaTime, std::vector<double>& inDiag, std::vector<double>& inX, std::vector<double>& inY);
	void InitializeLinearSystemSparse(float deltaTime, Eigen::SparseMatrix<double>& A);

	void CalculatePreconditioner(std::vector<double>& inOutPrecon, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY);

	void ApplyA(float deltaTime, Eigen::VectorXd& outResult, const Eigen::VectorXd &inVec, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY);
	void ApplyPreconditioner(Eigen::VectorXd& outResult, const Eigen::VectorXd& inResidual, const std::vector<double>& inPrecon, const std::vector<double>& inX, const std::vector<double>& inY);

	int mNumCellWidth;
	int mNumCellHeight;

	float dLeft;
	float dBottom;

	std::vector<glm::vec2> mCellCenters;

	std::vector<double> mCellXVelocities; // X velocity of the cell, staggered to the cells left edge.
	std::vector<double> mCellYVelocities; // Y velocity of the cell, staggered to the cells bottom edge.

	std::vector<double> mCellXVelocitiesPrev;
	std::vector<double> mCellYVelocitiesPrev;

	// Intermediate cell velocities
	std::vector<double> mIntXVelocities;
	std::vector<double> mIntYVelocities;

	std::vector<double> mCellPressures;
};