#include "Domain.h"

Domain::Domain(const glm::vec3& inCenter, const glm::vec3& inExtent)
{
	mCenter = inCenter;
	mExtent = inExtent;

	mLength = fabsf(mExtent.z * 2.0f);
	mWidth = fabsf(mExtent.x * 2.0f);
	mHeight = fabsf(mExtent.y * 2.0f);

	CalculateSides();
}

Domain::Domain(const glm::vec3& inCenter, float inLength, float inWidth, float inHeight)
{
	mCenter = inCenter;
	mLength = inLength;
	mWidth = inWidth;
	mHeight = inHeight;

	CalculateSides();

	mExtent = { mRight, mTop, mFront };
}

void Domain::CalculateSides()
{
	// X axis
	mLeft = mCenter.x - mWidth * 0.5f;
	mRight = mCenter.x + mWidth * 0.5f;
	// Y axis
	mBottom = mCenter.y - mHeight * 0.5f;
	mTop = mCenter.y + mHeight * 0.5f;
	// Z axis
	mBack = mCenter.z - mLength * 0.5f;
	mFront = mCenter.z + mLength * 0.5f;
}

bool Domain::IsPointInDomain(const glm::vec3& point)
{
	bool inside = true;

	if(point.x <= mLeft)
	{
		inside = false;
	}
	else if (point.x >= mRight)
	{
		inside = false;
	}
	else if (point.y <= mBottom)
	{
		inside = false;
	}
	else if (point.y >= mTop)
	{
		inside = false;
	}
	else if (point.z <= mBack)
	{
		inside = false;
	}
	else if (point.z >= mFront)
	{
		inside = false;
	}

	return inside;
}