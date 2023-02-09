#include "Simulation2D.h"

#include <iostream>

Simulation2D::Simulation2D(ApplicationData& inOutData)
{
	// Initialise new staggered Mac grid.
	MACGrid2D grid(inOutData);
	mMACGrid = grid;

	mMACGrid.UpdateApplicationData(inOutData);

	// Initialise new fluid.
	Fluid2D fluid(inOutData);
	mFluid = fluid;

	mFluid.UpdateApplicationData(inOutData);
}

void Simulation2D::StepSimulation(ApplicationData& inOutData)
{
	float deltaTime = inOutData.GetDeltaTime();

//	mFluid.InterpolateToGrid(mMACGrid);
	mFluid.InterpolateToGridBSpline(mMACGrid);

//	mMACGrid.Advect(inOutData);
	mMACGrid.ExtrapolateVelocityField(true);
	mMACGrid.ApplyForces(deltaTime);
	
	mMACGrid.Project(inOutData);
	mMACGrid.ExtrapolateVelocityField(false);

	mFluid.InterpolateFromGrid(mMACGrid);

	//mFluid.StepParticles(deltaTime);
	mFluid.StepParticlesRK3(deltaTime, mMACGrid);

	mFluid.UpdateApplicationData(inOutData);
	mMACGrid.UpdateCellTypesFromParticles(inOutData.Get2DParticlePositions());

	
	mMACGrid.UpdateApplicationData(inOutData);
}