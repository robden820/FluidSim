#pragma once

#include <vector>
#include <string>
#include "glm/glm.hpp"

// Class that holds all data relevant to different parts of the application.
// All application components read/write data here, so reduces dependencies on each other.

enum CellType
{
	eFLUID = 0,
	eSOLID = 1,
	eAIR = 2,
	eNONE = 3
};

class ApplicationData
{
	public:
		ApplicationData() = default;
		~ApplicationData() = default;

		// -------------------------------------------------- //
		// Simulation Data
	public:
		float GetDeltaTime() { return mDeltaTime; }
		void SetDeltaTime(float inDeltaTime) { mDeltaTime = inDeltaTime; }

	private:
		float mDeltaTime;

		// -------------------------------------------------- //
		// Fluid Data
	public:
		int GetNumParticles() const { return mNumParticles; }
		void SetNumParticles(int inNumParticles);

		const std::vector<glm::vec2>& Get2DParticlePositions() const { return mParticlePositions2D; }
		const glm::vec2& Get2DParticlePosition(int index) const;
		void Set2DParticlePositions(const std::vector<glm::vec2>& inPositions);

		const std::vector<glm::vec3>& Get3DParticlePositions() const { return mParticlePositions3D; }
		const glm::vec3& Get3DParticlePosition(int index) const;
		void Set3DParticlePositions(const std::vector<glm::vec3>& inPositions);

	private:
		int mNumParticles;

		std::vector<glm::vec2> mParticlePositions2D;
		std::vector<glm::vec3> mParticlePositions3D;

	// -------------------------------------------------- //
	// Grid Data
	public:
		int GetNumGridCells() const { return mNumGridCellsWidth * mNumGridCellsHeight * mNumGridCellsLength; }
		int GetNumGridCellsWidth() const { return mNumGridCellsWidth; }
		int GetNumGridCellsHeight() const { return mNumGridCellsHeight; }
		int GetNumGridCellsLength() const { return mNumGridCellsLength; }

		void SetNumGridCellsWidth(int inNumCellsWidth) { mNumGridCellsWidth = inNumCellsWidth; }
		void SetNumGridCellsHeight(int inNumCellsHeight) { mNumGridCellsHeight = inNumCellsHeight; }
		void SetNumGridCellsLength(int inNumCellsLength) { mNumGridCellsLength = inNumCellsLength; }

		float GetGridLeft() const { return mGridLeft; }
		float GetGridBottom() const { return mGridBottom; }
		float GetGridBack() const { return mGridBack; }

		void SetGridLeft(float inGridLeft) { mGridLeft = inGridLeft; }
		void SetGridBottom(float inGridBottom) { mGridBottom = inGridBottom; }
		void SetGridBack(float inGridBack) { mGridBack = inGridBack; }

		float GetGridCellSize() const { return mGridCellSize; }
		void SetGridCellSize(float inGridCellSize) { mGridCellSize = inGridCellSize; }

		float GetGridWidth() const { return mNumGridCellsWidth * mGridCellSize; }
		float GetGridHeight() const { return mNumGridCellsHeight * mGridCellSize; }
		float GetGridLength() const { return mNumGridCellsLength * mGridCellSize; }

		const std::vector<CellType>& GetCellTypes() const { return mCellTypes; }
		CellType GetCellType(int index) const;
		void SetCellTypes(const std::vector<CellType>& inCellTypes);

	private:
		int mNumGridCellsWidth;
		int mNumGridCellsHeight;
		int mNumGridCellsLength;

		float mGridLeft;
		float mGridBottom;
		float mGridBack;

		float mGridCellSize; // Only support grid with regular sized cells;

		std::vector<CellType> mCellTypes;

	// -------------------------------------------------- //
	// Rendering Data
	public:

	private:

	// -------------------------------------------------- //
	// Misc
	private:
		void LogError(const std::string& inErrorMsg) const;
};