#pragma once

#include <vector>

#include "glm/glm.hpp"
#include "onetbb/oneapi/tbb.h"

#include "ApplicationData.h"

class MACGrid
{
	public:

		virtual ~MACGrid() = default;

		virtual void Update(ApplicationData& inOutData) = 0;

		virtual int GetNumCells() const = 0;

		double GetCellPressure(int index) const { return mCellPressures[index]; }
		void SetCellPressure(int index, double inPressure) { mCellPressures[index]= inPressure; }

		const CellType GetCellType(int index) const { return mCellType[index]; }
		void SetCellType(int index, CellType inCellType) { mCellType[index] = inCellType; }

		double GetCellSize() { return mCellSize; }

	protected:

		virtual void InitializeGrid(const ApplicationData& inData) = 0;

		virtual void CalculateCellDivergence() = 0;

		virtual void AdvectCellVelocity(double deltaTime) = 0;
		virtual void UpdateCellPressure(double deltaTime, int maxIterations) = 0;
		virtual void UpdateCellVelocity(double deltaTime) = 0;

		int mNumCells;

		double mCellSize;    // deltaX
		double mInvCellSize; // 1 / deltaX

		double mDensity;
		double mInvDensity; // 1 / density;

		std::vector<double> mCellDivergence;

		std::vector<double> mCellPressures;

		std::vector<CellType> mCellType;
};