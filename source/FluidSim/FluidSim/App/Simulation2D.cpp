#include "Simulation2D.h"

#include <iostream>
#include "GLFW/glfw3.h"

#include "FLIPFluid2D.h"
#include "APICFluid2D.h"
#include "RPICFluid2D.h"

Simulation2D::Simulation2D(ApplicationData& inOutData)
{
	// Initialise new staggered Mac grid.
	MACGrid2D grid(inOutData);
	mMACGrid = grid;

	mMACGrid.UpdateApplicationData(inOutData);

	const int fluidType = 2;

	// Initialise new fluid.
	if (fluidType == 0)
	{
		FLIPFluid2D fluid(inOutData);
		mFluid = std::make_unique<FLIPFluid2D>(fluid);
	}
	else if (fluidType == 1)
	{
		APICFluid2D fluid(inOutData);
		mFluid = std::make_unique<APICFluid2D>(fluid);
	}
	else if (fluidType == 2)
	{
		RPICFluid2D fluid(inOutData);
		mFluid = std::make_unique<RPICFluid2D>(fluid);
	}

	mFluid->UpdateApplicationData(inOutData);
}

void Simulation2D::StepSimulation(ApplicationData& inOutData)
{
	double start = glfwGetTime();

	double deltaTime = inOutData.GetDeltaTime();

	// Interpolate the velocities of the particles to the grid.
	mFluid->InterpolateToGrid(mMACGrid);

	std::cout << "Interpolate to grid: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	// Extrapolate the velocities of the fluid cells to the entire grid (excluding solid cells).
	// - Set to true as we are extrapolating the intermediate cell velocities.
	mMACGrid.ExtrapolateVelocityField(true);

	// Save the velocity field for use later when calculating velocity change over the timestep.
	mMACGrid.SavePreviousVelocities();

	std::cout << "Extrapolate velocity field: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	// Apply forces to the grid (currently just gravitational force).
	mMACGrid.ApplyForces(deltaTime);

	std::cout << "Apply forces: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();
	
	// Calculate cell divergence, update cell pressures to enforce zero divergence, update cell velocities based on new pressures.
	mMACGrid.Project(inOutData);

	std::cout << "Project fluid: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	// Extrapolate the velocities of the fluid cells to the entire grid (excluding solid cells).
	// - Set to false as we are extrapolating the final cell velocities.
	mMACGrid.ExtrapolateVelocityField(false);

	std::cout << "Extrapolate velocity field: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	// Calculate the difference in final velocity vs initial velocity for this time step.
	mMACGrid.CalculateVelocityChange();

	// Interpolate the velocities of the grid to the particles.
	mFluid->InterpolateFromGridBSpline(mMACGrid);

	std::cout << "Interpolate from grid: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	// Advance the particles through the velocity field.
	mFluid->StepParticles(deltaTime, mMACGrid);
	//mFluid->StepParticlesEuler(deltaTime, mMACGrid);

	std::cout << "Step particles (Euler): " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	// Delete any particles that have left the grid.
	mFluid->DeleteBadParticles(mMACGrid);

	std::cout << "Remove bad particles: " << glfwGetTime() - start << "\n";
	start = glfwGetTime();

	// Update all the data for the simulation.
	mFluid->UpdateApplicationData(inOutData);
	mMACGrid.UpdateCellTypesFromParticles(inOutData.Get2DParticlePositions());
	mMACGrid.UpdateApplicationData(inOutData);

	std::cout << "Update data: " << glfwGetTime() - start << "\n";

	mFluid->CalculateSystemEnergy();
}