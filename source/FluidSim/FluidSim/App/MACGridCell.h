#pragma once

#include <vector>

#include "glm/glm.hpp"

class MACGridCell
{
	enum Face
	{
		eLEFT = 0,
		eRIGHT = 1,
		eTOP = 2,
		eBOTTOM = 3,
		eBACK = 4,
		eFRONT = 5,
		eNONE = 0
	};

	public:
		MACGridCell();
		~MACGridCell() = default;

		void SetFaceValidity(Face face, bool isValid) { mFaceValid[face] = isValid; }

		const float GetPressure() const { return mPressure; }
		void SetPressure(float inPressure) { mPressure = inPressure; }

		const std::vector<glm::vec3>& GetVelocities() const { return mFaceVelocities; }
		const glm::vec3& GetVelocity(Face face) const { return mFaceVelocities[face]; }

		void EnforceZeroDivergence();

	private:
		
		int mNumFaces;

		float mPressure;

		std::vector<glm::vec3> mFaceVelocities;
		std::vector<bool> mFaceValid; // A face is invalid if there is no neighboring grid cell in that direction.
};