#pragma once

#include <tuple>

class Domain2D
{
public:
	Domain2D() = default;
	Domain2D(const float inCenterX, const float inCenterY, float inWidth, float inHeight);

	~Domain2D() = default;

	bool IsPointInDomain2D(const float pointX, const float pointY);

	const std::tuple<float, float> GetCenter() const { return std::tuple<float, float>(mCenterX, mCenterY); }
	const std::tuple<float, float> GetExtent() const { return std::tuple<float, float>(mExtentX, mExtentY); }

	float GetWidth() const { return mWidth; }
	float GetHeight() const { return mHeight; }

	float GetBottom() const { return mBottom; }
	float GetTop() const { return mTop; }
	float GetLeft() const { return mLeft; }
	float GetRight() const { return mRight; }

private:
	void CalculateSides();

	float mCenterX;
	float mCenterY;

	float mExtentX;
	float mExtentY;

	float mWidth;  // x axis
	float mHeight; // y axis

	float mBottom;
	float mTop;
	float mLeft;
	float mRight;
};