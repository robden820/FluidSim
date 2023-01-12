#include "ApplicationData.h"

#include <iostream>

const glm::vec2& ApplicationData::Get2DParticlePosition(int index) const
{
	if (index < 0 || index > mParticlePositions2D.size() - 1)
	{
		LogError("Attempting to access invalid 2D particle position.");
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
		LogError("Attempting to access invalid 3D particle position.");
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

void ApplicationData::LogError(const std::string& inErrorMsg) const
{
	std::cout << "ERROR: " << inErrorMsg << "\n";
}