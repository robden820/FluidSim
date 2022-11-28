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
		eNONE = 2
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

		const std::vector<glm::vec3>& GetVelocities() const { return mFaceVelocities; }
		const glm::vec3& GetVelocity(Face inFace) const { return mFaceVelocities[inFace]; }
		void SetVelocity(Face inFace, const glm::vec3& inVelocity) { mFaceVelocities[inFace] = inVelocity; }

		const CellType& GetCellType() const { return mCellType; }
		void SetCellType(CellType inCellType) { mCellType = inCellType; }

		void EnforceZeroDivergence();

	private:
		
		int mNumFaces;

		float mPressure;

		std::vector<glm::vec3> mFaceVelocities;
		std::vector<bool> mFaceValid; // A face is invalid if there is no neighboring grid cell in that direction.

		CellType mCellType;
};