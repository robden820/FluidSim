#pragma once

#include "Contact.h"
#include <memory>

/*
* Contact resolution for particle contacts.
* One resolver instance is required for the whole physics simulation.
*/
class ContactResolver
{
public:

	ContactResolver(int iterations = 0);

	~ContactResolver() = default;

	void SetIterations(int iterations) { mIterations = iterations; }

	void ResolveContacts(Contact* contactArray, int numContacts, float deltaTime);

protected:

	int mIterations; // Number of allowed iterations for contact resolution.
};

