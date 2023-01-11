#pragma once

#include "Domain.h"

#include "glm/vec3.hpp"

class Domain3D : public Domain
{
public:
	Domain3D() = default;
	Domain3D(const glm::vec3& inCenter, const glm::vec3& inExtent);
	Domain3D(const glm::vec3& inCenter, float inLength, float inWidth, float inHeight);

	~Domain3D() = default;

	bool IsPointInDomain(const glm::vec3& point);

	const glm::vec3& GetCenter() const { return mCenter; }
	const glm::vec3& GetExtent() const { return mExtent; }

	float GetLength() const { return mLength; }
	float GetWidth() const { return mWidth; }
	float GetHeight() const { return mHeight; }

	float GetBottom() const { return mBottom; }
	float GetTop() const { return mTop; }
	float GetLeft() const { return mLeft; }
	float GetRight() const { return mRight; }
	float GetFront() const { return mFront; }
	float GetBack() const { return mBack; }

private:
	void CalculateSides() override;

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
