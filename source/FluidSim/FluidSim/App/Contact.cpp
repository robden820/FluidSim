
#include "Contact.h"
#include "glm/glm.hpp"

Contact::Contact(Particle& particle0, Particle& particle1, float restitution, float penetrationDepth, glm::vec3 contactNormal)
{
	mParticle[0] = &particle0;
	mParticle[1] = &particle1;

	mRestitution = restitution;
	mPenetrationDepth = penetrationDepth;
	mContactNormal = contactNormal;

}

void Contact::Resolve(float deltaTime)
{
	ResolveVelocity(deltaTime);
	ResolveInterpenetration(deltaTime);
}

float Contact::CalculateSeparationVelocity() const
{
	glm::vec3 relativeVelocity = mParticle[0]->GetVelocity();

	if (mParticle[1])
	{
		relativeVelocity -= mParticle[1]->GetVelocity();
	}

	return glm::dot(relativeVelocity, mContactNormal);
}

void Contact::ResolveVelocity(float deltaTime)
{
	float separationVelocity = CalculateSeparationVelocity();

	// Check if anything needs to be resolved.
	if (separationVelocity > 0)
	{
		// No impulse required as we are either separating or stationary.
		return;
	}

	// Calculate separation velocity after the collision.
	float newSeparationVelocity = -1 * separationVelocity * mRestitution;

	// Check the velocity change due to acceleration only.
	glm::vec3 accelerationCausedVelocity = mParticle[0]->GetAcceleration();

	if (mParticle[1])
	{
		accelerationCausedVelocity += mParticle[1]->GetAcceleration();
	}

	float accelerationCausedSepVel = glm::dot(accelerationCausedVelocity,(mContactNormal)) * deltaTime;

	// If we have velocity due to acceleration build up, remove it from the new seperation velocity.
	// This is to deal with resting contacts.
	if (accelerationCausedSepVel < 0)
	{
		newSeparationVelocity += mRestitution * accelerationCausedSepVel;

		// Don't reverse the direction of velocity.
		if (newSeparationVelocity < 0)
		{
			newSeparationVelocity = 0;
		}
	}

	float deltaVelocity = newSeparationVelocity - separationVelocity;

	// Apply the change in velocity proportional to the mass of each particle.
	float p1InvMass = 1 / mParticle[0]->GetMass();
	float p2InvMass = 1 / mParticle[1]->GetMass();

	float totalInverseMass = p1InvMass;

	if (mParticle[1])
	{
		totalInverseMass += p2InvMass;
	}

	if (totalInverseMass <= 0)
	{
		// If both particles have infinite mass then the impulses will have no effect.
		return;
	}

	// Calculate the impulse to be applied.
	float impulseStrength = deltaVelocity / totalInverseMass;

	glm::vec3 impulse = mContactNormal * impulseStrength;

	//Apply impulse to each particle remembering contact normal is from perspective of particle 0.
	glm::vec3 particle0NewVelocity = mParticle[0]->GetVelocity() + impulse * p1InvMass;
	mParticle[0]->SetVelocity(particle0NewVelocity);

	if (mParticle[1])
	{
		glm::vec3 particle1NewVelocity = mParticle[1]->GetVelocity() + impulse * p2InvMass * -1.0f;
		mParticle[1]->SetVelocity(particle1NewVelocity);
	}
}

void Contact::ResolveInterpenetration(float deltaTime)
{
	// Check if objects are penetrating.
	if (mPenetrationDepth <= 0.0f)
	{
		return;
	}

	// Apply the change in velocity proportional to the mass of each particle.
	float p1InvMass = 1 / mParticle[0]->GetMass();
	float p2InvMass = 1 / mParticle[1]->GetMass();

	float totalInverseMass = p1InvMass;

	if (mParticle[1])
	{
		totalInverseMass += p2InvMass;
	}

	if (totalInverseMass <= 0)
	{
		// If both particles have infinite mass then the impulses will have no effect.
		return;
	}

	glm::vec3 deltaPositionPerInvMass = mContactNormal * (-1 * mPenetrationDepth / totalInverseMass);

	// Apply delta position to each particle remembering contact normal is from perspective of particle 0.
	glm::vec3 particle0NewPosition = mParticle[0]->GetPosition() + deltaPositionPerInvMass * p1InvMass;
	mParticle[0]->SetPosition(particle0NewPosition);

	if (mParticle[1])
	{
		glm::vec3 particle1NewPosition = mParticle[1]->GetPosition() + deltaPositionPerInvMass * p2InvMass;
		mParticle[1]->SetPosition(particle1NewPosition);
	}
}