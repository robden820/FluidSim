#pragma once

#include <algorithm>
#include <vector>
#include <string>
#include "glm/glm.hpp"

// Class that holds all data relevant to different parts of the application.
// All application components read/write data here, so reduces dependencies on each other.

enum class CellType : int
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
		double GetDeltaTime() { return mDeltaTime; }
		void SetDeltaTime(double inDeltaTime) { mDeltaTime = inDeltaTime; }

	private:
		double mDeltaTime;

		// -------------------------------------------------- //
		// Fluid Data
	public:
		int GetNumParticles() const { return mNumParticles; }
		void SetNumParticles(int inNumParticles);

		double GetFluidDensity() const { return mFluidDensity; }
		void SetFluidDensity(double inDensity) { mFluidDensity = inDensity; }

		double GetFLIPBlend() const { return mFLIPBlend; }
		void SetFLIPBlend(double inBlend) { mFLIPBlend = std::clamp(inBlend, 0.0, 1.0); }

		const std::vector<double>& GetParticleSpeeds() const { return mParticleSpeeds; }
		const double GetParticleSpeed(int index) const;
		void SetParticleSpeeds(const std::vector<double>& inSpeeds);

		const std::vector<glm::dvec2>& Get2DParticlePositions() const { return mParticlePositions2D; }
		const glm::dvec2& Get2DParticlePosition(int index) const;
		void Set2DParticlePositions(const std::vector<glm::dvec2>& inPositions);

		const std::vector<glm::dvec3>& Get3DParticlePositions() const { return mParticlePositions3D; }
		const glm::dvec3& Get3DParticlePosition(int index) const;
		void Set3DParticlePositions(const std::vector<glm::dvec3>& inPositions);

	private:
		int mNumParticles;

		double mFluidDensity;

		double mFLIPBlend;

		std::vector<double> mParticleSpeeds;

		std::vector<glm::dvec2> mParticlePositions2D;
		std::vector<glm::dvec3> mParticlePositions3D;

	// -------------------------------------------------- //
	// Grid Data
	public:
		int GetNumGridCells() const { return mNumGridCells; }
		int GetNumGridCellsWidth() const { return mNumGridCellsWidth; }
		int GetNumGridCellsHeight() const { return mNumGridCellsHeight; }
		int GetNumGridCellsLength() const { return mNumGridCellsLength; }

		void SetNumGridCellsWidth(int inNumCellsWidth) { mNumGridCellsWidth = inNumCellsWidth; }
		void SetNumGridCellsHeight(int inNumCellsHeight) { mNumGridCellsHeight = inNumCellsHeight; }
		void SetNumGridCellsLength(int inNumCellsLength) { mNumGridCellsLength = inNumCellsLength; }
		void UpdateNumGridCells(); // Call after setting number of grid cells in each dimension to reserve vector memory.

		double GetGridLeft() const { return mGridLeft; }
		double GetGridBottom() const { return mGridBottom; }
		double GetGridBack() const { return mGridBack; }

		void SetGridLeft(double inGridLeft) { mGridLeft = inGridLeft; }
		void SetGridBottom(double inGridBottom) { mGridBottom = inGridBottom; }
		void SetGridBack(double inGridBack) { mGridBack = inGridBack; }

		double GetGridCellSize() const { return mGridCellSize; }
		void SetGridCellSize(double inGridCellSize) { mGridCellSize = inGridCellSize; }

		double GetGridWidth() const { return mNumGridCellsWidth * mGridCellSize; }
		double GetGridHeight() const { return mNumGridCellsHeight * mGridCellSize; }
		double GetGridLength() const { return mNumGridCellsLength * mGridCellSize; }

		const std::vector<CellType>& GetCellTypes() const { return mCellTypes; }
		CellType GetCellType(int index) const;
		void SetCellTypes(const std::vector<CellType>& inCellTypes);

		const std::vector<glm::dvec2>& GetCellCenters2D() const { return mCellCenters2D; }
		const glm::dvec2& GetCellCenter2D(int index) const;
		void SetCellCenters2D(const std::vector<glm::dvec2>& inCellCenters);

		const std::vector<glm::dvec3>& GetCellCenters3D() const { return mCellCenters3D; }
		const glm::dvec3& GetCellCenter3D(int index) const;
		void SetCellCenters3D(const std::vector<glm::dvec3>& inCellCenters);

	private:
		int mNumGridCellsWidth = 1;
		int mNumGridCellsHeight = 1;
		int mNumGridCellsLength = 1;

		int mNumGridCells;

		double mGridLeft;
		double mGridBottom;
		double mGridBack;

		double mGridCellSize; // Only support grid with regular sized cells;

		std::vector<CellType> mCellTypes;
		std::vector<glm::dvec2> mCellCenters2D;
		std::vector<glm::dvec3> mCellCenters3D;

	// -------------------------------------------------- //
	// Rendering Data
	public:

	private:

	// -------------------------------------------------- //
	// Misc
	private:
		void LogError(const std::string& inErrorMsg) const;
};