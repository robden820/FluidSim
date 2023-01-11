#include "Domain2D.h"

Domain2D::Domain2D(const glm::vec2& inCenter, const glm::vec2& inExtent)
{
	mCenter = inCenter;
	mExtent = inExtent;

	mWidth = fabsf(mExtent.x * 2.0f);
	mHeight = fabsf(mExtent.y * 2.0f);

	CalculateSides();
}

Domain2D::Domain2D(const glm::vec2& inCenter, float inWidth, float inHeight)
{
	mCenter = inCenter;

	mWidth = inWidth;
	mHeight = inHeight;

	CalculateSides();

	mExtent.x = mRight;
	mExtent.y = mTop;
}

void Domain2D::CalculateSides()
{
	// X axis
	mLeft = mCenter.x - mWidth * 0.5f;
	mRight = mLeft + mWidth;
	// Y axis
	mBottom = mCenter.y - mHeight * 0.5f;
	mTop = mBottom + mHeight;
}

bool Domain2D::IsPointInDomain(const glm::vec2& inPoint)
{
	bool inside = true;

	if (inPoint.x <= mLeft)
	{
		inside = false;
	}
	else if (inPoint.x >= mRight)
	{
		inside = false;
	}
	else if (inPoint.y <= mBottom)
	{
		inside = false;
	}
	else if (inPoint.y >= mTop)
	{
		inside = false;
	}

	return inside;
}