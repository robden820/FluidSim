#pragma once

#include "Particle.h"

class Contact
{
public:

	Contact(Particle& particle0, Particle& particle1, float restitution, float penetrationDepth, glm::vec3 contactNormal);

	~Contact() = default;

	Particle* mParticle[2]; // The particles involved in the contact. The second can be null.

	float mRestitution; // Coeffecient of restitution for the contact.

	glm::vec3 mContactNormal; // Direction of contact normal from the perspective of particle 0.

	float mPenetrationDepth; // Depth of interpenetration of contact.

	void Resolve(float deltaTime); // Resolves contact for velocity and interpenetration.

	float CalculateSeparationVelocity() const; // Calculates separation velocity of the contact.

private:

	void ResolveVelocity(float deltaTime); // Resolves the impulse calculations.

	void ResolveInterpenetration(float deltaTime); // Resolves interpenetration of the contact.
};