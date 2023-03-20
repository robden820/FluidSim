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
	void ApplyForces(double deltaTime);
	void Project(ApplicationData& inOutData);

	void ExtrapolateVelocityField(bool extrapolateIntVelocities);

	int GetNumCells() const { return mNumCellHeight * mNumCellWidth; }
	int GetNumCellsWidth() const { return mNumCellWidth; }
	int GetNumCellsHeight() const { return mNumCellHeight; }

	const glm::dvec2& GetCellCenter(int index) const { return mCellCenters[index]; }
	int GetClosestCell(const glm::dvec2& inPos) const;

	const double GetCellXVelocity(int index) const { return mCellXVelocities[index]; }
	const double GetCellYVelocity(int index) const { return mCellYVelocities[index]; }

	const double GetCellXVelocityDiff(int index) const { return mCellXVelocitiesDiff[index]; }
	const double GetCellYVelocityDiff(int index) const { return mCellYVelocitiesDiff[index]; }

	const double GetCellXVelocityPrev(int index) const { return mCellXVelocitiesPrev[index]; }
	const double GetCellYVelocityPrev(int index) const { return mCellYVelocitiesPrev[index]; }

	const glm::dvec2& GetCellVelocity(int index) const;

	void SetCellXVelocity(int index, double inVelocity) { mCellXVelocities[index] = inVelocity; }
	void SetCellYVelocity(int index, double inVelocity) { mCellYVelocities[index] = inVelocity; }

	void SetIntXVelocity(int index, double inVelocity) { mIntXVelocities[index] = inVelocity; }
	void SetIntYVelocity(int index, double inVelocity) { mIntYVelocities[index] = inVelocity; }

	const double GetCellPressure(int index) const { return mCellPressures[index]; }
	void SetCellPressure(int index, double inPressure) { mCellPressures[index] = inPressure; }

	const double GetCellMassX(int index) const { return mCellMassX[index]; }
	void SetCellMassX(int index, double inMass) { mCellMassX[index] = inMass; }

	const double GetCellMassY(int index) const { return mCellMassY[index]; }
	void SetCellMassY(int index, double inMass) { mCellMassY[index] = inMass; }

	const CellType GetCellType(int index) const { return mCellType[index]; }
	const std::vector<CellType>& GetCellTypes() const { return mCellType; }
	const CellType GetCellTypeFromPosition (const glm::dvec2& inPos) const;
	void SetCellType(int index, CellType inCellType) { mCellType[index] = inCellType; }

	std::tuple<int, int> GetXYFromIndex(int index) const;
	int GetIndexFromXY(int X, int Y) const;

	double GetCellSize() const { return mCellSize; }
	double GetInverseCellSize() const { return mInvCellSize; }

	void UpdateCellTypesFromParticles(const std::vector<glm::dvec2>& inParticlePositions);

	void UpdateApplicationData(ApplicationData& inOutData);

	void SavePreviousVelocities();
	void CalculateVelocityChange();

private:

	void InitializeGrid(const ApplicationData& inData) override;
	void InitializeGridPressure();

	void CalculateCellDivergence();

	void UpdateCellPressure(double deltaTime, int maxIterations);
	void UpdateCellPressureSparse(double deltaTime, int maxIterations);

	void AdvectCellVelocity(double deltaTime);
	void UpdateCellVelocity(double deltaTime);

	void InitializeLinearSystem(double deltaTime, std::vector<double>& inDiag, std::vector<double>& inX, std::vector<double>& inY);
	void InitializeLinearSystemSparse(double deltaTime, Eigen::SparseMatrix<double>& A);

	void CalculatePreconditioner(std::vector<double>& inOutPrecon, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY);

	void ApplyA(double deltaTime, Eigen::VectorXd& outResult, const Eigen::VectorXd &inVec, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY);
	void ApplyPreconditioner(Eigen::VectorXd& outResult, const Eigen::VectorXd& inResidual, const std::vector<double>& inPrecon, const std::vector<double>& inX, const std::vector<double>& inY);

	int mNumCellWidth;
	int mNumCellHeight;

	double dLeft;
	double dBottom;

	std::vector<glm::dvec2> mCellCenters;

	std::vector<double> mCellXVelocities; // X velocity of the cell, staggered to the cells left edge.
	std::vector<double> mCellYVelocities; // Y velocity of the cell, staggered to the cells bottom edge.

	std::vector<double> mCellXVelocitiesDiff;
	std::vector<double> mCellYVelocitiesDiff;

	std::vector<double> mCellXVelocitiesPrev;
	std::vector<double> mCellYVelocitiesPrev;

	// Intermediate cell velocities
	std::vector<double> mIntXVelocities;
	std::vector<double> mIntYVelocities;

	std::vector<double> mCellPressures;

	// Mass contributed by particles to each face of the cell. Same face as the velocity.
	std::vector<double> mCellMassX;
	std::vector<double> mCellMassY;
};