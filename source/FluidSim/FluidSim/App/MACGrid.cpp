# include "MACGrid.h"

MACGrid::MACGrid(const Domain& inDomain, float inCellSize)
{

}

void MACGrid::Update(float deltaTime)
{
	UpdateCellVelocity(deltaTime);

	UpdateCellPressure(deltaTime);
}

void MACGrid::UpdateCellVelocity(float deltaTime)
{

}

void MACGrid::UpdateCellPressure(float deltaTime)
{

}