#pragma once

#include <vector>

#include "glm/glm.hpp"

class MACGridCell
{
public:
	enum CellType
	{
		eFLUID = 0,
		eSOLID = 1,
		eEMPTY = 2
	};

	enum Face
	{
		eRIGHT = 0, // x axis
		eTOP = 1,   // y axis
		eFRONT = 2, // z axis
		eNONE = 3
	};

		MACGridCell();
		~MACGridCell() = default;

		void SetFaceValidity(int inFace, bool isValid) { mFaceValid[inFace] = isValid; }

		const float GetPressure() const { return mPressure; }
		void SetPressure(float inPressure) { mPressure = inPressure; }

		float GetFaceVelocity(Face inFace) const { return mFaceVelocities[inFace]; }
		void SetFaceVelocity(Face inFace, float inVelocity) { mFaceVelocities[inFace] = inVelocity; }

		const glm::vec3& GetCellVelocity() const { return mFaceVelocities; }
		void SetCellVelocity(const glm::vec3& inVelocity) { mFaceVelocities = inVelocity; }

		const CellType& GetCellType() const { return mCellType; }
		void SetCellType(CellType inCellType) { mCellType = inCellType; }

	private:
		
		int mNumFaces;

		float mPressure;

		glm::vec3 mFaceVelocities;
		std::vector<bool> mFaceValid; // A face is invalid if there is no neighboring grid cell in that direction.

		CellType mCellType;
};