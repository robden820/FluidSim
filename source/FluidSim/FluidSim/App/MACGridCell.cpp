#include "MACGridCell.h"

#define ERROR 0.0001f

MACGridCell::MACGridCell()
{
	// Only need to store velocities on 3 positive axis faces, as neighbours will store velocities of other faces.
	mNumFaces = 3;

	mPressure = 1.f;

	mFaceVelocities.reserve(mNumFaces);
	mFaceValid.reserve(mNumFaces);

	for (int face = 0; face < mNumFaces; face++)
	{
		mFaceVelocities.push_back(glm::vec3(0.f));

		// Initialise all faces as invalid, but must be updated during MAC grid creation.
		mFaceValid.push_back(false);
	}
}

void MACGridCell::GetCellVelocity(glm::vec3& outVelocity) const
{
	outVelocity = glm::vec3(mFaceVelocities[eRIGHT].x, mFaceVelocities[eTOP].y, mFaceVelocities[eFRONT].z);
}

void MACGridCell::SetCellVelocity(const glm::vec3& inVelocity)
{
	SetVelocity(eRIGHT, glm::vec3(inVelocity.x, 0.f, 0.f));
	SetVelocity(eTOP, glm::vec3(0.f, inVelocity.y, 0.f));
	SetVelocity(eFRONT, glm::vec3(0.f, 0.f, inVelocity.z));
}

void MACGridCell::EnforceZeroDivergence()
{
	glm::vec3 divergence(0.f);

	for (int face = 0; face < mNumFaces; face++)
	{
		divergence += mFaceVelocities[face];
	}

	if (glm::length(divergence) < ERROR)
	{
		return;
	}
}