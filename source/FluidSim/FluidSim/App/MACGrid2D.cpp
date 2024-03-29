#include "MACGrid2D.h"

#include <iostream>

#define TOLERANCE 0.0000001
#define PRECON_TUNER 0.97
#define PRECON_SAFETY 0.25

MACGrid2D::MACGrid2D()
{
	dLeft = 0.0;
	dBottom = 0.0;

	mNumCellWidth = 0;
	mNumCellHeight = 0;
}

MACGrid2D::MACGrid2D(const ApplicationData& inData)
{
	InitializeGrid(inData);
	InitializeGridPressure();
}

void MACGrid2D::InitializeGrid(const ApplicationData& inData)
{
	dLeft = inData.GetGridLeft();
	dBottom = inData.GetGridBottom();

	mNumCellWidth = inData.GetNumGridCellsWidth();
	mNumCellHeight = inData.GetNumGridCellsHeight();

	mNumCells = inData.GetNumGridCells();

	mCellSize = inData.GetGridCellSize();
	mInvCellSize = 1 / mCellSize;

	double halfCell = mCellSize * 0.5;

	mCellCenters.reserve(mNumCells);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		double centerX = dLeft + (x * mCellSize) + halfCell;
		double centerY = dBottom + (y * mCellSize) + halfCell;

		mCellCenters.push_back(glm::dvec2(centerX, centerY));
	}

	mCellPressures.assign(mNumCells, 0.0);
	mCellXVelocities.assign(mNumCells, 0.0);
	mCellYVelocities.assign(mNumCells, 0.0);
	mCellXVelocitiesDiff.assign(mNumCells, 0.0);
	mCellYVelocitiesDiff.assign(mNumCells, 0.0);
	mCellDivergence.assign(mNumCells, 0.0);
	mCellMassX.assign(mNumCells, 0.0);
	mCellMassY.assign(mNumCells, 0.0);
	
	mIntXVelocities.assign(mNumCells, 0.0);
	mIntYVelocities.assign(mNumCells, 0.0);

	mDensity = inData.GetFluidDensity();
	mInvDensity = 1.0 / mDensity;

	//Initialize cell types and solid boundary.
	mCellType.assign(mNumCells, CellType::eNONE);

	// TO DO: initialize solid boundary and initial fluid cells from file.
	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (x < 1 || x > mNumCellWidth - 2 || y < 1 || y > mNumCellHeight - 2)
		{
			mCellType[index] = CellType::eSOLID;
		}

		if (x > 10 && x < 26 && y > 20 && y < 70)
		{
			mCellType[index] = CellType::eFLUID;
		}
	}
}

void MACGrid2D::InitializeGridPressure()
{
	for (int index = 0; index < mNumCells; index++)
	{
		if (mCellType[index] == CellType::eFLUID)
		{
			mCellPressures[index] = 1.0;
		}
	}
}

void MACGrid2D::UpdateCellTypesFromParticles(const std::vector<glm::dvec2>& inParticlePositions)
{
	tbb::parallel_for(0, mNumCells, 1, [&](int cIndex)
	{
		if (mCellType[cIndex] != CellType::eSOLID)
		{
			mCellType[cIndex] = CellType::eAIR;
		}
	});

	// For each particle, find the closest cell and mark it as a fluid.
	tbb::parallel_for(0, (int)inParticlePositions.size(), 1, [&](int pIndex)
	{
		int cellIndex = GetClosestCell(inParticlePositions[pIndex]);

		if (cellIndex >= 0 && mCellType[cellIndex] != CellType::eSOLID)
		{
			mCellType[cellIndex] = CellType::eFLUID;
		}
	});
}

void MACGrid2D::Update(ApplicationData& inOutData)
{
	double deltaTime = inOutData.GetDeltaTime();

	Advect(inOutData);

	ApplyForces(deltaTime);

	Project(inOutData);
}

void MACGrid2D::Advect(ApplicationData& inOutData)
{
	double deltaTime = inOutData.GetDeltaTime();

	AdvectCellVelocity(deltaTime);
}

void MACGrid2D::Project(ApplicationData& inOutData)
{
	double deltaTime = inOutData.GetDeltaTime();

	CalculateCellDivergence();

	//UpdateCellPressure(deltaTime, 200);
	UpdateCellPressureSparse(deltaTime, 200);

	UpdateCellVelocity(deltaTime);
}

void MACGrid2D::AdvectCellVelocity(double deltaTime)
{
	tbb::parallel_for(0, mNumCells, 1, [&](int index)
	{
		if (mCellType[index] == CellType::eFLUID)
		{
			int x, y;
			std::tie(x, y) = GetXYFromIndex(index);

			glm::dvec2 cellVelocity = { mCellXVelocities[index], mCellYVelocities[index] };

			//Runge Kutta
			// First take a half step backwards over half the timestep to find the average velocity
			glm::dvec2 halfPrevPosition = mCellCenters[index] - cellVelocity * 0.5 * deltaTime;

			int midStepIndex = GetClosestCell(halfPrevPosition);

			int midStepX, midStepY;
			std::tie(midStepX, midStepY) = GetXYFromIndex(midStepIndex);

			int midStepNeighbourRight = GetIndexFromXY(midStepX + 1, midStepY);
			int midStepNeighbourTop = GetIndexFromXY(midStepX, midStepY + 1);

			// Compute the avgerage velocity by interpolating the velocity of the mid step position.
			glm::dvec2 diff = mCellCenters[midStepIndex] - halfPrevPosition;

			glm::dvec2 weight = (diff * mInvCellSize) + 0.5;

			double avgStepXVelocity = (mCellXVelocities[midStepIndex] * weight.x) + (mCellXVelocities[midStepNeighbourRight] * (1 - weight.x));
			double avgStepYVelocity = (mCellYVelocities[midStepIndex] * weight.y) + (mCellYVelocities[midStepNeighbourTop] * (1 - weight.y));

			glm::dvec2 avgStepVelocity = { avgStepXVelocity, avgStepYVelocity };

			// Now that we have the average velocity of the step, we can work out the original position at the beginning of the step.
			glm::dvec2 prevPosition = mCellCenters[index] - avgStepVelocity * deltaTime;

			int prevCellIndex = GetClosestCell(prevPosition);

			int prevX, prevY;
			std::tie(prevX, prevY) = GetXYFromIndex(prevCellIndex);

			int prevNeighbourRight = GetIndexFromXY(prevX + 1, prevY);
			int prevNeighbourTop = GetIndexFromXY(prevX, prevY + 1);

			// Recalculate our weights for interpolation
			diff = mCellCenters[prevCellIndex] - prevPosition;

			weight = (diff * mInvCellSize) + 0.5;

			double prevXVelocity = (mCellXVelocities[prevCellIndex] * weight.x) + (mCellXVelocities[prevNeighbourRight] * (1 - weight.x));
			double prevYVelocity = (mCellYVelocities[prevCellIndex] * weight.y) + (mCellYVelocities[prevNeighbourTop] * (1 - weight.y));

			mIntXVelocities[index] = prevXVelocity;
			mIntYVelocities[index] = prevYVelocity;
		}
		else
		{
			mIntXVelocities[index] = 0.0;
			mIntYVelocities[index] = 0.0;
		}
	});
}

int MACGrid2D::GetClosestCell(const glm::dvec2& inPos) const
{
	int x = static_cast<int>(round((inPos.x - dLeft - (mCellSize * 0.5)) * mInvCellSize));
	int y = static_cast<int>(round((inPos.y - dBottom - (mCellSize * 0.5)) * mInvCellSize));

	int approxIndex = GetIndexFromXY(x, y);

	if (approxIndex >= mNumCells || approxIndex < 0)
	{
		approxIndex = -1;
	}

	return approxIndex;
}

void MACGrid2D::ApplyForces(double deltaTime)
{
	double gravity = -9.8 * deltaTime;

	for (int index = 0; index < mNumCells; index++)
	{
		if (mCellType[index] != CellType::eSOLID)
		{
			mIntYVelocities[index] += gravity;
		}
	}
}

void MACGrid2D::UpdateCellVelocity(double deltaTime)
{
	double scale = (double)deltaTime * mInvCellSize * mInvDensity;

	double solidXVel = 0.0;
	double solidYVel = 0.0;

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (x > 0)
		{
			int neighbourLeft = GetIndexFromXY(x - 1, y);

			if (mCellType[neighbourLeft] == CellType::eFLUID || mCellType[index] == CellType::eFLUID)
			{
				if (mCellType[neighbourLeft] == CellType::eSOLID || mCellType[index] == CellType::eSOLID)
				{
					mCellXVelocities[index] = solidXVel;
				}
				else
				{
					mCellXVelocities[index] = mIntXVelocities[index] - scale * (mCellPressures[index] - mCellPressures[neighbourLeft]);
				}
			}
		}
		if (y > 0)
		{
			int neighbourBottom = GetIndexFromXY(x, y - 1);

			if (mCellType[neighbourBottom] == CellType::eFLUID || mCellType[index] == CellType::eFLUID)
			{
				if (mCellType[neighbourBottom] == CellType::eSOLID || mCellType[index] == CellType::eSOLID)
				{
					mCellYVelocities[index] = solidYVel;
				}
				else
				{
					mCellYVelocities[index] = mIntYVelocities[index] - scale * (mCellPressures[index] - mCellPressures[neighbourBottom]);
				}
			}
		}
	}
}

void MACGrid2D::ExtrapolateVelocityField(bool extrapolateIntVelocities)
{
	int maxSearchDepth = mNumCells;
	std::vector<int> marker;
	marker.assign(mNumCells, maxSearchDepth);

	std::vector<int> wavefrontIndices;
	wavefrontIndices.reserve(mNumCells);

	// Mark known velocities as zero.
	for (int cellIndex = 0; cellIndex < mNumCells; cellIndex++)
	{
		if (mCellType[cellIndex] == CellType::eFLUID)
		{
			marker[cellIndex] = 0;
		}
	}

	// Initialise first wave for search.
	for (int cellIndex = 0; cellIndex < mNumCells; cellIndex++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(cellIndex);

		int numFluidNeighbours = 0;

		if (marker[cellIndex] == 0)
		{
			continue;
		}

		if (x > 0)
		{
			int left = GetIndexFromXY(x - 1, y);
			if (mCellType[left] == CellType::eFLUID)
			{
				numFluidNeighbours++;
			}
		}
		if (x < mNumCellWidth - 1)
		{
			int right = GetIndexFromXY(x + 1, y);
			if (mCellType[right] == CellType::eFLUID)
			{
				numFluidNeighbours++;
			}
		}
		if (y > 0)
		{
			int bottom = GetIndexFromXY(x, y - 1);
			if (mCellType[bottom] == CellType::eFLUID)
			{
				numFluidNeighbours++;
			}
		}
		if (y < mNumCellHeight - 1)
		{
			int top = GetIndexFromXY(x, y + 1);
			if (mCellType[top] == CellType::eFLUID)
			{
				numFluidNeighbours++;
			}
		}

		// If at least one neighbour is a fluid, add the cell to the first wavefront.
		if (numFluidNeighbours > 0)
		{
			marker[cellIndex] = 1;
			wavefrontIndices.push_back(cellIndex);
		}
	}

	int iteration = 0;

	// Breadth first search extrapolation.
	while (iteration < wavefrontIndices.size())
	{
		int cellIndex = wavefrontIndices[iteration];
		int cellMarker = marker[cellIndex];

		int x, y;
		std::tie(x, y) = GetXYFromIndex(cellIndex);

		int numSearchedNeighbours = 0;
		double sumXVel = 0;
		double sumYVel = 0;

		if (x > 0)
		{
			int left = GetIndexFromXY(x - 1, y);

			if (marker[left] < cellMarker)
			{
				if (extrapolateIntVelocities)
				{
					sumXVel += mIntXVelocities[left];
					sumYVel += mIntYVelocities[left];
				}
				else
				{
					sumXVel += mCellXVelocities[left];
					sumYVel += mCellYVelocities[left];
				}

				++numSearchedNeighbours;
			}
			else if(marker[left] == maxSearchDepth && mCellType[left] != CellType::eSOLID)
			{
				marker[left] = cellMarker + 1;
				wavefrontIndices.push_back(left);
			}
		}

		if (x < mNumCellWidth - 1)
		{
			int right = GetIndexFromXY(x + 1, y);

			if (marker[right] < cellMarker)
			{
				if (extrapolateIntVelocities)
				{
					sumXVel += mIntXVelocities[right];
					sumYVel += mIntYVelocities[right];
				}
				else
				{
					sumXVel += mCellXVelocities[right];
					sumYVel += mCellYVelocities[right];
				}

				++numSearchedNeighbours;
			}
			else if(marker[right] == maxSearchDepth && mCellType[right] != CellType::eSOLID)
			{
				marker[right] = cellMarker + 1;
				wavefrontIndices.push_back(right);
			}
		}

		if (y > 0)
		{
			int bottom = GetIndexFromXY(x, y - 1);

			if (marker[bottom] < cellMarker)
			{
				if (extrapolateIntVelocities)
				{
					sumXVel += mIntXVelocities[bottom];
					sumYVel += mIntYVelocities[bottom];
				}
				else
				{
					sumXVel += mCellXVelocities[bottom];
					sumYVel += mCellYVelocities[bottom];
				}
		
				++numSearchedNeighbours;
			}
			else if(marker[bottom] == maxSearchDepth && mCellType[bottom] != CellType::eSOLID)
			{
				marker[bottom] = cellMarker + 1;
				wavefrontIndices.push_back(bottom);
			}
		}
		
		if (y < mNumCellHeight - 1)
		{
			int top = GetIndexFromXY(x, y + 1);

			if (marker[top] < cellMarker)
			{
				if (extrapolateIntVelocities)
				{
					sumXVel += mIntXVelocities[top];
					sumYVel += mIntYVelocities[top];
				}
				else
				{
					sumXVel += mCellXVelocities[top];
					sumYVel += mCellYVelocities[top];
				}

				++numSearchedNeighbours;
			}
			else if(marker[top] == maxSearchDepth && mCellType[top] != CellType::eSOLID)
			{
				marker[top] = cellMarker + 1;
				wavefrontIndices.push_back(top);
			}
		}

		if (numSearchedNeighbours > 0)
		{
			if (extrapolateIntVelocities)
			{
				SetIntXVelocity(cellIndex, sumXVel / numSearchedNeighbours);
				SetIntYVelocity(cellIndex, sumYVel / numSearchedNeighbours);
			}
			else
			{
				SetCellXVelocity(cellIndex, sumXVel / numSearchedNeighbours);
				SetCellYVelocity(cellIndex, sumYVel / numSearchedNeighbours);
			}
		}

		++iteration;
	}
}

void MACGrid2D::CalculateCellDivergence()
{
	double scale = -mInvCellSize;

	double solidXVel = 0.0;
	double solidYVel = 0.0;

	// Calculate negative divergence of each cell.
	for (int index = 0; index < mNumCells; index++)
	{
		mCellDivergence[index] = 0.0;

		if (mCellType[index] == CellType::eFLUID)
		{
			int x, y;
			std::tie(x, y) = GetXYFromIndex(index);

			double divergence = 0.0;

			int neighbourRight = GetIndexFromXY(x + 1, y);
			if (mCellType[neighbourRight] != CellType::eSOLID)
			{
				divergence += mIntXVelocities[neighbourRight] - mIntXVelocities[index];
			}
			
			int neighbourTop = GetIndexFromXY(x, y + 1);
			if (mCellType[neighbourTop] != CellType::eSOLID)
			{
				divergence += mIntYVelocities[neighbourTop] - mIntYVelocities[index];
			}
			
			divergence *= scale;

			mCellDivergence[index] = divergence;
		}
	}
	
	// Update divergence to account for solid cells.
	for (int index = 0; index < mNumCells; index++)
	{
		if (mCellType[index] == CellType::eFLUID)
		{
			int x, y;
			std::tie(x, y) = GetXYFromIndex(index);

			int neighbourLeft = GetIndexFromXY(x - 1, y);

			if (mCellType[neighbourLeft] == CellType::eSOLID)
			{
				mCellDivergence[index] += scale * (mIntXVelocities[index] - solidXVel);
			}

			int neighbourRight = GetIndexFromXY(x + 1, y);

			if (mCellType[neighbourRight] == CellType::eSOLID)
			{
				mCellDivergence[index] -= scale * (mIntXVelocities[neighbourRight] - solidXVel);
			}

			int neighbourBottom = GetIndexFromXY(x, y - 1);

			if (mCellType[neighbourBottom] == CellType::eSOLID)
			{
				mCellDivergence[index] += scale * (mIntYVelocities[index] - solidYVel);
			}

			int neighbourTop = GetIndexFromXY(x, y + 1);

			if (mCellType[neighbourTop] == CellType::eSOLID)
			{
				mCellDivergence[index] -= scale * (mIntYVelocities[neighbourTop] - solidYVel);
			}
		}
	}
}

void MACGrid2D::UpdateCellPressure(double deltaTime, int maxIterations)
{
	std::vector<double> Adiagonal;
	std::vector<double> Ax;
	std::vector<double> Ay;

	Adiagonal.assign(mNumCells, 0.0);
	Ax.assign(mNumCells, 0.0);
	Ay.assign(mNumCells, 0.0);

	InitializeLinearSystem(deltaTime, Adiagonal, Ax, Ay);

	// Preconditioned Conjugate Gradient.
	Eigen::VectorXd newPressure(mNumCells);
	newPressure.fill(0.0);

	Eigen::VectorXd residuals(mNumCells);
	residuals.fill(0.0);
	// TO DO: make cell divergence an eigen vector.
	for (int i = 0; i < mNumCells; i++)
	{
		residuals(i) = mCellDivergence[i];
	}

	double maxResidual = -1.0;

	for (int index = 0; index < mNumCells; index++)
	{
		if (abs(residuals[index]) > maxResidual)
		{
			maxResidual = abs(residuals[index]);
		}
	}

	if (maxResidual < TOLERANCE)
	{
		// Reached convergence;
		std::cout << "Zero cell divergence. Skipping pressure update. \n";
		return;
	}

	Eigen::VectorXd z(mNumCells);
	z.fill(0.0);

	std::vector<double> precon;
	precon.assign(mNumCells, 0.0);
	CalculatePreconditioner(precon, Adiagonal, Ax, Ay);

	ApplyPreconditioner(z, residuals, precon, Ax, Ay);

	Eigen::VectorXd search = z;

	double sigma = z.dot(residuals);

	int iteration;

	for (iteration = 0; iteration < maxIterations; ++iteration)
	{
		ApplyA(deltaTime, z, search, Adiagonal, Ax, Ay);

		double phi = z.dot(search);
		double alpha = sigma / phi;

		// Update pressure and residual vectors.
		residuals -= z * alpha;
		newPressure += search * alpha;

		// Calculate maximum residual error.
		maxResidual = -1.0;

		for (int index = 0; index < mNumCells; index++)
		{
			if (abs(residuals[index]) > maxResidual)
			{
				maxResidual = abs(residuals[index]);
			}
		}

		if (maxResidual < TOLERANCE)
		{
			// Reached convergence;
			break;
		}

		ApplyPreconditioner(z, residuals, precon, Ax, Ay);

		// SigmaNew = dot product of z and residuals.
		double sigmaNew = z.dot(residuals);
		double beta = sigmaNew / sigma;

		// Update search vector for next iteration.
		search *= beta;
		search += z;

		sigma = sigmaNew;
	}
	std::cout << "------------------------------\n";
	if (iteration >= maxIterations)
	{
		std::cout << "WARNING: MAX NUMBER OF ITERATIONS REACHED IN PRESSURE SOLVE: ";
	}

	std::cout << "NUM PRESSURE SOLVE ITERATIONS: " << iteration << "\n";
	std::cout << "PRESSURE SOLVE ERROR: " << maxResidual << "\n";

	// Set pressure
	for (int i = 0; i < mNumCells; i++)
	{
		mCellPressures[i] = newPressure[i];
	}
}

void MACGrid2D::InitializeLinearSystem(double deltaTime, std::vector<double>& inDiag, std::vector<double>& inX, std::vector<double>& inY)
{
	double scale = (double)deltaTime * mInvCellSize * mInvCellSize * mInvDensity;

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y; 
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			// Handle left neighbour
			if (x > 0)
			{
				int neighbourLeft = GetIndexFromXY(x - 1, y);

				if (mCellType[neighbourLeft] == CellType::eFLUID)
				{
					inDiag[index] += scale;
				}
			}
			// Handle right neighbour
			if (x < mNumCellWidth - 1)
			{
				int neighbourRight = GetIndexFromXY(x + 1, y);

				if (mCellType[neighbourRight] == CellType::eFLUID)
				{
					inDiag[index] += scale;
					inX[index] -= scale;
				}
				else if (mCellType[neighbourRight] != CellType::eSOLID)
				{
					inDiag[index] += scale;
				}
			}
			// Handle bottom neighbour
			if (y > 0)
			{
				int neighbourBottom = GetIndexFromXY(x, y - 1);

				if (mCellType[neighbourBottom] == CellType::eFLUID)
				{
					inDiag[index] += scale;
				}
			}
			// Handle top neighbour
			if (y < mNumCellHeight - 1)
			{
				int neighbourTop = GetIndexFromXY(x, y + 1);

				if (mCellType[neighbourTop] == CellType::eFLUID)
				{
					inDiag[index] += scale;
					inY[index] -= scale;
				}
				else if (mCellType[neighbourTop] != CellType::eSOLID)
				{
					inDiag[index] += scale;
				}
			}
		}
	}
}

void MACGrid2D::UpdateCellPressureSparse(double deltaTime, int maxIterations)
{
	Eigen::SparseMatrix<double> A(mNumCells, mNumCells);

	InitializeLinearSystemSparse(deltaTime, A);

	Eigen::VectorXd divergence(mNumCells);
	Eigen::VectorXd pressure(mNumCells);

	for (int i = 0; i < mNumCells; i++)
	{
		divergence[i] = mCellDivergence[i];
	}

	// Use full matrix rather than just upper as in RBridson as this achieves better performance in Eigen.
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;

	solver.setMaxIterations(maxIterations);
	solver.setTolerance(TOLERANCE);

	solver.compute(A);

	pressure = solver.solve(divergence);

	for (int i = 0; i < mNumCells; i++)
	{
		mCellPressures[i] = pressure[i];
	}

	std::cout << "NUM PRESSURE SOLVE ITERATIONS: " << solver.iterations() << "\n";
	std::cout << "PRESSURE SOLVE ERROR: " << solver.error() << "\n";
}

void MACGrid2D::InitializeLinearSystemSparse(double deltaTime, Eigen::SparseMatrix<double>& A)
{
	double scale = (double)deltaTime * mInvCellSize * mInvCellSize * mInvDensity;

	std::vector<Eigen::Triplet<double>> coefficients; // TO DO: reserve memory. Num fluid cells?

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{

			int numFluidNeighbours = 0;

			// Handle left neighbour
			if (x > 0)
			{
				int neighbourLeft = GetIndexFromXY(x - 1, y);

				if (mCellType[neighbourLeft] == CellType::eFLUID)
				{
					++numFluidNeighbours;
					coefficients.push_back(Eigen::Triplet<double>(neighbourLeft, index, -scale));
				}
				else if (mCellType[neighbourLeft] != CellType::eSOLID)
				{
					++numFluidNeighbours;
				}
			}
			// Handle right neighbour
			if (x < mNumCellWidth - 1)
			{
				int neighbourRight = GetIndexFromXY(x + 1, y);

				if (mCellType[neighbourRight] == CellType::eFLUID)
				{
					++numFluidNeighbours;
					coefficients.push_back(Eigen::Triplet<double>(neighbourRight, index, -scale));
				}
				else if (mCellType[neighbourRight] != CellType::eSOLID)
				{
					++numFluidNeighbours;
				}
			}
			// Handle bottom neighbour
			if (y > 0)
			{
				int neighbourBottom = GetIndexFromXY(x, y - 1);

				if (mCellType[neighbourBottom] == CellType::eFLUID)
				{
					++numFluidNeighbours;
					coefficients.push_back(Eigen::Triplet<double>(neighbourBottom, index, -scale));
				}
				else if (mCellType[neighbourBottom] != CellType::eSOLID)
				{
					++numFluidNeighbours;
				}
			}
			// Handle top neighbour
			if (y < mNumCellHeight - 1)
			{
				int neighbourTop = GetIndexFromXY(x, y + 1);

				if (mCellType[neighbourTop] == CellType::eFLUID)
				{
					++numFluidNeighbours;
					coefficients.push_back(Eigen::Triplet<double>(neighbourTop, index, -scale));
				}
				else if (mCellType[neighbourTop] != CellType::eSOLID)
				{
					++numFluidNeighbours;
				}
			}


			coefficients.push_back(Eigen::Triplet<double>(index, index, scale * numFluidNeighbours));
		}
	}

	A.setFromTriplets(coefficients.begin(), coefficients.end());
}

void MACGrid2D::ApplyA(double deltaTime, Eigen::VectorXd& outResult, const Eigen::VectorXd& inVec, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY)
{
	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		double value = 0.0;

		if (x > 0)
		{
			int neighbourLeft = GetIndexFromXY(x - 1, y);
			value += inVec[neighbourLeft] * inX[neighbourLeft];
		}

		if (x < mNumCellWidth - 1)
		{
			int neighbourRight = GetIndexFromXY(x + 1, y);  
			value += inVec[neighbourRight] * inX[index];

		}

		if (y > 0)
		{
			int neighbourBottom = GetIndexFromXY(x, y - 1);
			value += inVec[neighbourBottom] * inY[neighbourBottom];
		}

		if (y < mNumCellHeight - 1)
		{
			int neighbourTop = GetIndexFromXY(x, y + 1);
			value += inVec[neighbourTop] * inY[index];
		}

  		value += inDiag[index] * inVec[index];
		outResult[index] = value;
	}
}

void MACGrid2D::CalculatePreconditioner(std::vector<double>& inOutPrecon, const std::vector<double>& inDiag, const std::vector<double>& inX, const std::vector<double>& inY)
{
	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			double Axi = 0.0; // Aplusi_iminus1
			double Axj = 0.0; // Aplusi_jminus1

			double Ayi = 0.0; // Aplusj_iminus1
			double Ayj = 0.0; // Aplusj_jminus1

			double preconi = 0.0; // precon_iminus1
			double preconj = 0.0; // precon_jminus1

			if (x > 0)
			{
				int neighbourLeft = GetIndexFromXY(x - 1, y);

				Axi = inX[neighbourLeft];
				Ayi = inY[neighbourLeft];

				preconi = inOutPrecon[neighbourLeft];
			}

			if (y > 0)
			{
				int neighbourBottom = GetIndexFromXY(x, y - 1);

				Axj = inX[neighbourBottom];
				Ayj = inY[neighbourBottom];

				preconj = inOutPrecon[neighbourBottom];
			}

			double a = Axi * preconi;
			double b = Ayj * preconj;

			double termOne = a * a + b * b;

			double d = Axi * Ayi * preconi * preconi;
			double e = Ayj * Axj * preconj * preconj;

			double termTwo = d + e;

			double newPrecon = inDiag[index] - termOne - PRECON_TUNER * termTwo;

			if (newPrecon < PRECON_SAFETY * inDiag[index])
			{
				newPrecon = inDiag[index];
			}

			inOutPrecon[index] = 1.0 / sqrt(newPrecon);
		}
	}
}

void MACGrid2D::ApplyPreconditioner(Eigen::VectorXd& outResult, const Eigen::VectorXd& inResidual, const std::vector<double>& inPrecon, const std::vector<double>& inX, const std::vector<double>& inY)
{
	std::vector<double> intermediate;  // q
	intermediate.assign(mNumCells, 0.f);

	for (int index = 0; index < mNumCells; index++)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{

			double Axi = 0.0; // Aplusi_iminus1
			double Ayj = 0.0; // Aplusj_jminus1

			double preconi = 0.0; // precon_iminus1
			double preconj = 0.0; // precon_jminus1

			double intermediatei = 0.0; // qminusi
			double intermediatej = 0.0; // qminusj

			int neighbourLeft = GetIndexFromXY(x - 1, y);

			Axi = inX[neighbourLeft];
			preconi = inPrecon[neighbourLeft];
			intermediatei = intermediate[neighbourLeft];

			int neighbourBottom = GetIndexFromXY(x, y - 1);

			Ayj = inY[neighbourBottom];
			preconj = inPrecon[neighbourBottom];
			intermediatej = intermediate[neighbourBottom];

			double t = inResidual[index] - Axi * preconi * intermediatei
										- Ayj * preconj * intermediatej;

			intermediate[index] = t * inPrecon[index];
		}
	}

	for (int index = mNumCells - 1; index >= 0; index--)
	{
		int x, y;
		std::tie(x, y) = GetXYFromIndex(index);

		if (mCellType[index] == CellType::eFLUID)
		{
			double zi = 0.0;
			double zj = 0.0;

			int neighbourRight = GetIndexFromXY(x + 1, y);

			zi = outResult[neighbourRight];

			int neighbourTop = GetIndexFromXY(x, y + 1);

			zj = outResult[neighbourTop];

			double t = intermediate[index] - inX[index] * inPrecon[index] * zi
										  - inY[index] * inPrecon[index] * zj;

			outResult[index] = t * inPrecon[index];
		}
	}
}

const glm::dvec2& MACGrid2D::GetCellVelocity(int index) const
{                                                        
	int x, y;
	std::tie(x, y) = GetXYFromIndex(index);

	if (x >= mNumCellWidth - 1 || y >= mNumCellHeight - 1)
	{
		return glm::dvec2(0.0, 0.0);
	}

	int neighbourRight = GetIndexFromXY(x + 1, y);
	int neighbourTop = GetIndexFromXY(x, y + 1);

	double xVelocity = GetCellXVelocity(index) + GetCellXVelocity(neighbourRight);
	double yVelocity = GetCellYVelocity(index) + GetCellYVelocity(neighbourTop);

	return glm::dvec2(xVelocity, yVelocity);
}

void MACGrid2D::SavePreviousVelocities()
{
	mCellXVelocitiesPrev = mIntXVelocities;
	mCellYVelocitiesPrev = mIntYVelocities;
}

void MACGrid2D::CalculateVelocityChange()
{
	for (int cellIndex = 0; cellIndex < mNumCells; cellIndex++)
	{
		mCellXVelocitiesDiff[cellIndex] = mCellXVelocities[cellIndex] - mCellXVelocitiesPrev[cellIndex];
		mCellYVelocitiesDiff[cellIndex] = mCellYVelocities[cellIndex] - mCellYVelocitiesPrev[cellIndex];
	}
}

const CellType MACGrid2D::GetCellTypeFromPosition(const glm::dvec2& inPos) const
{
	int index = GetClosestCell(inPos);

	return GetCellType(index);
}

std::tuple<int, int> MACGrid2D::GetXYFromIndex(int index) const
{
	int x = 0;
	int y = 0;

	y = index % mNumCellHeight;
	x = static_cast<int>(round(index / mNumCellHeight)) % mNumCellWidth;

	return std::tuple<int, int>(x, y);
}

int MACGrid2D::GetIndexFromXY(int X, int Y) const
{
	return  Y + X * mNumCellHeight;
}

void MACGrid2D::UpdateApplicationData(ApplicationData& inOutData)
{
	inOutData.SetCellTypes(mCellType);
	inOutData.SetCellCenters2D(mCellCenters);
}