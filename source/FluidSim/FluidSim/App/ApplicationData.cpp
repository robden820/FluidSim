#include "ApplicationData.h"

#include <iostream>

void ApplicationData::SetNumParticles(int inNumParticles)
{
	mNumParticles = inNumParticles;

	mParticlePositions2D.reserve(mNumParticles);
	mParticlePositions3D.reserve(mNumParticles);
}

const glm::vec2& ApplicationData::Get2DParticlePosition(int index) const
{
	if (index < 0 || index > mParticlePositions2D.size() - 1)
	{
		LogError("Attempting to access invalid 2D particle position index.");
		return glm::vec2{ 0.f, 0.f };
	}

	return mParticlePositions2D[index];
}

void ApplicationData::Set2DParticlePositions(const std::vector<glm::vec2>& inPositions)
{
	mParticlePositions2D.clear();
	mParticlePositions2D.reserve(inPositions.size());

	mParticlePositions2D = inPositions;
}

const glm::vec3& ApplicationData::Get3DParticlePosition(int index) const
{
	if (index < 0 || index > mParticlePositions3D.size() - 1)
	{
		LogError("Attempting to access invalid 3D particle position index.");
		return glm::vec3{ 0.f, 0.f, 0.f};
	}

	return mParticlePositions3D[index];
}

void ApplicationData::Set3DParticlePositions(const std::vector<glm::vec3>& inPositions)
{
	mParticlePositions3D.clear();
	mParticlePositions3D.reserve(inPositions.size());

	mParticlePositions3D = inPositions;
}

void ApplicationData::UpdateNumGridCells()
{
	mNumGridCells = mNumGridCellsWidth * mNumGridCellsHeight * mNumGridCellsLength;

	mCellTypes.clear();
	mCellTypes.reserve(mNumGridCells);
}

CellType ApplicationData::GetCellType(int index) const
{
	if (index < 0 || index > mCellTypes.size() - 1)
	{
		LogError("Attempting to access invalid cell type index.");
		return eNONE;
	}

	return mCellTypes[index];
}

void ApplicationData::SetCellTypes(const std::vector<CellType>& inCellTypes)
{
	mCellTypes.clear();
	mCellTypes.reserve(inCellTypes.size());

	mCellTypes = inCellTypes;
}

const glm::vec2& ApplicationData::GetCellCenter2D(int index) const
{
	if (index < 0 || index > mCellCenters2D.size() - 1)
	{
		LogError("Attempting to access invalid 2D cell center index.");
		return glm::vec2(0.0f, 0.0f);
	}

	return mCellCenters2D[index];
}

void ApplicationData::SetCellCenters2D(const std::vector<glm::vec2>& inCellCenters)
{
	mCellCenters2D.clear();
	mCellCenters2D.reserve(inCellCenters.size());

	mCellCenters2D = inCellCenters;
}

const glm::vec3& ApplicationData::GetCellCenter3D(int index) const
{
	if (index < 0 || index > mCellCenters3D.size() - 1)
	{
		LogError("Attempting to access invalid 2D cell center index.");
		return glm::vec3(0.0f, 0.0f, 0.0f);
	}

	return mCellCenters3D[index];
}

void ApplicationData::SetCellCenters3D(const std::vector<glm::vec3>& inCellCenters)
{
	mCellCenters3D.clear();
	mCellCenters3D.reserve(inCellCenters.size());

	mCellCenters3D = inCellCenters;
}

void ApplicationData::LogError(const std::string& inErrorMsg) const
{
	std::cout << "ERROR: " << inErrorMsg << "\n";
}