#include "ContactResolver.h"

ContactResolver::ContactResolver(int iterations)
{
	mIterations = iterations;
}

void ContactResolver::ResolveContacts(Contact* contactArray, int numContacts, float deltaTime)
{
	int iterationsUsed = 0;

	while (iterationsUsed < mIterations)
	{
		// Want to deal with the most 'severe' contacts first.

		// Find contact with the largest closing velocity (most negative separating velocity.
		float max = 0;
		int maxIndex = numContacts - 1;

		for (int i = 0; i < numContacts; i++)
		{
			float separatingVel = contactArray[i].CalculateSeparationVelocity();

			if (separatingVel < max)
			{
				max = separatingVel;
				maxIndex = i;
			}
		}

		// Resolve the contact found.
		contactArray[maxIndex].Resolve(deltaTime);
		iterationsUsed++;
	}
}