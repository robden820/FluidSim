#include "GridNode.h"

GridNode::GridNode(const glm::vec3& inPosition)
{
	mPosition = inPosition;
	mVelocity = glm::vec3(0.f, 0.f, 0.f);
}