#pragma once

#include <memory>

#include "ApplicationData.h"

#include "Fluid2D.h"
#include "MACGrid2D.h"

class Simulation2D
{
	enum class FluidType : bool
	{
		eFLIP = false,
		eAPIC = true,
	};

public:
	Simulation2D(ApplicationData& inOutData);
	Simulation2D() = default;

	~Simulation2D() = default;

	void StepSimulation(ApplicationData& inOutData);

private:

	std::unique_ptr<IFluid2D> mFluid;
	MACGrid2D mMACGrid;
};