#pragma once

#include <vector>

#include "MACGrid.h"

#include "glm/glm.hpp"
#include "oneapi/tbb.h"

#include <Eigen/Sparse>

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

	const double GetCellXVelocity(int index) const { return mCellXVelocities[index]; }
	const double GetCellYVelocity(int index) const { return mCellYVelocities[index]; }

	void SetCellXVelocity(int index, double inVelocity) { mCellXVelocities[index] = inVelocity; }
	void SetCellYVelocity(int index, double inVelocity) { mCellYVelocities[index] = inVelocity; }

	const double GetCellPressure(int index) const { return mCellPressures[index]; }
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

	void UpdateCellPressureSparse(float deltaTime, int maxIterations);

	void CalculateCellDivergence(float deltaTime);

	void AdvectCellVelocity(float deltaTime);
	void UpdateCellPressure(float deltaTime, int maxIterations);
	void UpdateCellVelocity(float deltaTime);

	void InitializeLinearSystem(float deltaTime, std::vector<double>& inDiag, std::vector<double>& inX, std::vector<double>& inY);
	void InitializeLinearSystemSparse(float deltaTime, Eigen::SparseMatrix<float>& inOutA);

	void CalculatePreconditioner(std::vector<double>& inOutPrecon, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY);
	void CalculatePreconditionerSparse(Eigen::VectorXf& inOutPrecon, const Eigen::SparseMatrix<float>& inA);

	void ApplyA(float deltaTime, std::vector<double>& outResult, const std::vector<double>& inVec, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY);
	void ApplyPreconditioner(std::vector<double>& outResult, const std::vector<double>& inResidual, const std::vector<double>& inPrecon, const std::vector<double>& inX, const std::vector<double>& inY);

	int mNumCellWidth;
	int mNumCellHeight;

	float dLeft;
	float dBottom;

	std::vector<glm::vec2> mCellCenters;

	std::vector<double> mCellXVelocities; // X velocity of the cell, staggered to the cells left edge.
	std::vector<double> mCellYVelocities; // Y velocity of the cell, staggered to the cells bottom edge.

	// Intermediate cell velocities
	std::vector<double> mIntXVelocities;
	std::vector<double> mIntYVelocities;

	std::vector<double> mCellPressures;
};