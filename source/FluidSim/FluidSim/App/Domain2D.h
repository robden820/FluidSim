#pragma once

#include "Domain.h"

#include "glm\glm.hpp"

class Domain2D : public Domain
{
public:
	Domain2D() = default;
	Domain2D(const glm::vec2& inCenter, const glm::vec2& inExtent);
	Domain2D(const glm::vec2& inCenter, float inWidth, float inHeight);

	~Domain2D() = default;

	bool IsPointInDomain(const glm::vec2& inPoint);

	const glm::vec2& GetCenter() const { return mCenter; }
	const glm::vec2& GetExtent() const { return mExtent; }

	float GetWidth() const { return mWidth; }
	float GetHeight() const { return mHeight; }

	float GetBottom() const { return mBottom; }
	float GetTop() const { return mTop; }
	float GetLeft() const { return mLeft; }
	float GetRight() const { return mRight; }

private:
	void CalculateSides() override;

	glm::vec2 mCenter;
	glm::vec2 mExtent;

	float mWidth;  // x axis
	float mHeight; // y axis

	float mBottom;
	float mTop;
	float mLeft;
	float mRight;
};