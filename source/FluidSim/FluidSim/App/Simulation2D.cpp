#include "Simulation2D.h"

#include <iostream>
#include "GLFW/glfw3.h"

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
	double start = glfwGetTime();

	double deltaTime = inOutData.GetDeltaTime();

	mFluid.InterpolateToGrid(mMACGrid);

	std::cout << "Interpolate to grid: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	mMACGrid.ExtrapolateVelocityField(true);

	std::cout << "Extrapolate velocity field: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	mMACGrid.ApplyForces(deltaTime);

	std::cout << "Apply forces: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();
	
	mMACGrid.Project(inOutData);

	std::cout << "Project fluid: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	mMACGrid.ExtrapolateVelocityField(false);

	std::cout << "Extrapolate velocity field: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	mMACGrid.CalculateVelocityChange();

	mFluid.InterpolateFromGridBSpline(mMACGrid);

	std::cout << "Interpolate from grid: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	//mFluid.StepParticles(deltaTime, mMACGrid);
	mFluid.StepParticlesEuler(deltaTime, mMACGrid);

	std::cout << "Step particles (RK3): " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	mFluid.DeleteBadParticles(mMACGrid);

	std::cout << "Remove bad particles: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	mFluid.UpdateApplicationData(inOutData);
	mMACGrid.UpdateCellTypesFromParticles(inOutData.Get2DParticlePositions());
	mMACGrid.UpdateApplicationData(inOutData);

	std::cout << "Update data: " << glfwGetTime() - start << "\n";
}