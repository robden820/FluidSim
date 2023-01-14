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

		float GetCellPressure(int index) const { return mCellPressures[index]; }
		void SetCellPressure(int index, float inPressure) { mCellPressures[index]= inPressure; }

		const CellType GetCellType(int index) const { return mCellType[index]; }
		void SetCellType(int index, CellType inCellType) { mCellType[index] = inCellType; }

		float GetCellSize() { return mCellSize; }

	protected:

		virtual void InitializeGrid(const ApplicationData& inData) = 0;

		virtual void CalculateCellDivergence(float deltaTime) = 0;

		virtual void AdvectCellVelocity(float deltaTime) = 0;
		virtual void UpdateCellPressure(float deltaTime, int maxIterations) = 0;
		virtual void UpdateCellVelocity(float deltaTime) = 0;

		int mNumCells;

		float mCellSize;    // deltaX
		float mInvCellSize; // 1 / deltaX

		float mDensity;

		std::vector<float> mCellDivergence;

		std::vector<float> mCellPressures;

		std::vector<CellType> mCellType;
};