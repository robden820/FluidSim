#include "Domain2D.h"

Domain2D::Domain2D(const float inCenterX, const float inCenterY, float inWidth, float inHeight)
{
	mCenterX = inCenterX;
	mCenterY = inCenterY;

	mWidth = inWidth;
	mHeight = inHeight;

	CalculateSides();

	mExtentX = mRight;
	mExtentX = mTop;
}

void Domain2D::CalculateSides()
{
	// X axis
	mLeft = mCenterX - mWidth * 0.5f;
	mRight = mCenterX + mWidth * 0.5f;
	// Y axis
	mBottom = mCenterY - mHeight * 0.5f;
	mTop = mCenterY + mHeight * 0.5f;
}

bool Domain2D::IsPointInDomain2D(const float pointX, const float pointY)
{
	bool inside = true;

	if (pointX <= mLeft)
	{
		inside = false;
	}
	else if (pointX >= mRight)
	{
		inside = false;
	}
	else if (pointY <= mBottom)
	{
		inside = false;
	}
	else if (pointY >= mTop)
	{
		inside = false;
	}

	return inside;
}