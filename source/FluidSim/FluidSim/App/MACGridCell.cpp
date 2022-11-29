#include "MACGridCell.h"

#include <iostream>

MACGridCell::MACGridCell()
{
	// Only need to store velocities on 3 positive axis faces, as neighbours will store velocities of other faces.
	mNumFaces = 3;

	mPressure = 1.f;

	mFaceValid.reserve(mNumFaces);

	mFaceVelocities = { 0.f, 0.f, 0.f };

	for (int face = 0; face < mNumFaces; face++)
	{
		// Initialise all faces as invalid, but must be updated during MAC grid creation.
		mFaceValid.push_back(false);
	}

	mCellType = eEMPTY;
}