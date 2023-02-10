#include "Interpolation.h"

#include <iostream>

double Interpolation::CubicInterpolate(const double(&points)[4], double intervalSize, double x)
{
	double xSqr = x * x;
	double xCube = x * x * x;

	double h00 = 2.0 * xCube - 3.0 * xSqr + 1.0;
	double h10 = xCube - 2.0 * xSqr + x;
	double h01 = -2.0 * xCube + 3.0 * xSqr;
	double h11 = xCube - xSqr;

//	double m1 = FiniteDifference(points[1], points[2], intervalSize);
//	double m2 = FiniteDifference(points[2], points[3], intervalSize);

	double m1 = CardinalSpline(points[0], points[2], intervalSize, 0.5);
	double m2 = CardinalSpline(points[1], points[3], intervalSize, 0.5);

	return (h00 * points[1]) + (h10 * intervalSize * m1) + (h01 * points[2]) + (h11 * intervalSize * m2);
}

double Interpolation::BicubicInterpolate(const double(&points)[4][4], double intervalSize, glm::dvec2 x)
{
	double p0 = CubicInterpolate(points[0], intervalSize, x.x);
	double p1 = CubicInterpolate(points[1], intervalSize, x.x);
	double p2 = CubicInterpolate(points[2], intervalSize, x.x);
	double p3 = CubicInterpolate(points[3], intervalSize, x.x);

	double newPoints[4] = {p0, p1, p2, p3};

	return CubicInterpolate(newPoints, intervalSize, x.y);
}

double Interpolation::TricubicInterpolate()
{
	return 0.0;
}

double Interpolation::NCubicInterpolate()
{
	return 0.0;
}

double Interpolation::NBicubicInterpolate()
{
	return 0.0;
}

double Interpolation::NTricubicInterpolate()
{
	return 0.0;
}

double Interpolation::FiniteDifference(double point0,double point1, double intervalSize)
{
	double p0p1 = point1 - point0;
	double p1p0 = point0 - point1;

	return ((p0p1 / intervalSize) - (p1p0 / intervalSize)) * 0.5;
}

double Interpolation::CardinalSpline(double point0, double point2, double intervalSize, float tension)
{
	if (tension < 0 || tension > 1)
	{
		// Tension must be in range [0,1].
		std::cout << "WARNING: Tension value for cardinal spline not in range [0,1] \n";
		return 0.0;
	}

	double p0p2 = point2 - point0;

	return (1.0 - tension) * (p0p2 / (2.0 * intervalSize));
}

double Interpolation::CatmullRomSpline(double point0, double point2, double intervalSize)
{
	// Equivalent to cardinal spline with tension 0.5.
	return CardinalSpline(point0, point2, intervalSize, 0.5);
}

double Interpolation::NFiniteDifference(double point0, double point1)
{
	return point1 - point0;
}

double Interpolation::NCardinalSpline(double point0, double point2, float tension)
{
	if (tension < 0 || tension > 1)
	{
		// Tension must be in range [0,1].
		std::cout << "WARNING: Tension value for normalized cardinal spline not in range [0,1] \n";
		return 0.0;
	}

	double p0p2 = point2 - point0;

	// Interval range is [0,1] so interval between two points is [0,2].
	return (1.0 - tension) * p0p2 * 0.5;
}

double Interpolation::NCatmullRomSpline(double point0, double point2)
{
	// Equivalent to cardinal spline with tension 0.5.
	return NCardinalSpline(point0, point2, 0.5);
}