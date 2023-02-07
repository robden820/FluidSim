#include "Simulation2D.h"

#include <iostream>

Simulation2D::Simulation2D(ApplicationData& inOutData)
{
	// Initialise new staggered Mac grid.
	MACGrid2D grid(inOutData);
	mMACGrid = grid;

	// Initialise new fluid.
	Fluid2D fluid(inOutData);
	mFluid = fluid;

	inOutData.SetCellTypes(mMACGrid.GetCellTypes());
}

void Simulation2D::StepSimulation(ApplicationData& inOutData)
{
	float deltaTime = inOutData.GetDeltaTime();

	mFluid.InterpolateToGrid(mMACGrid);

//	mMACGrid.Advect(inOutData);

	mMACGrid.ApplyForces(deltaTime);

	mMACGrid.Project(inOutData);

	mFluid.InterpolateFromGrid(mMACGrid);

	//mFluid.StepParticles(deltaTime);
	mFluid.StepParticlesRK3(deltaTime, mMACGrid);

	mMACGrid.UpdateCellTypesFromParticles(inOutData.Get2DParticlePositions());

	mFluid.UpdateApplicationData(inOutData);
	mMACGrid.UpdateApplicationData(inOutData);
}