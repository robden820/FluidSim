#pragma once

#include "glm/vec3.hpp"

class Domain
{
	public:
		Domain() = default;
		Domain(const glm::vec3& inCenter, const glm::vec3& inExtent);
		Domain(const glm::vec3& inCenter, float inLength, float inWidth, float inHeight);

		~Domain() = default;

		bool IsPointInDomain(const glm::vec3& point);

		float GetLength() const { return mLength; }
		float GetWidth() const { return mWidth; }
		float GetHeight() const { return mHeight; }

		float GetBottom() { return mBottom; }
		float GetTop() { return mTop; }
		float GetLeft() { return mLeft; }
		float GetRight() { return mRight; }
		float GetFront() { return mFront; }
		float GetBack() { return mBack; }

	private:
		void CalculateSides();

		glm::vec3 mCenter; // Position of domain center

		glm::vec3 mExtent; // Position of domain corner.

		float mLength; // z axis
		float mWidth;  // x axis
		float mHeight; // y axis

		float mBottom;
		float mTop;
		float mLeft;
		float mRight;
		float mFront; // Side closest to camera.
		float mBack;  // Side furthest from camera.
};