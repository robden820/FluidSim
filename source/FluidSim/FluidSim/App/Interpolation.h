#pragma once

#include "glm/vec2.hpp"
#include "glm/vec3.hpp"

// Interpolation helper class.
namespace Interpolation
{
	// Used for non-normalized interpolate interval - [x,y]. Interval size is (y - x)
	double CubicInterpolate(const double (&points)[4], double intervalSize, double x);
	double BicubicInterpolate(const double (&points)[4][4], double intervalSize, glm::dvec2 x);
	double TricubicInterpolate();

	double FiniteDifference(double point0, double point1, double intervalSize); // Returns tangent gradient at point 0.
	double CardinalSpline(double point0, double point2, double intervalSize, double tension);   // Returns tangent gradient at point 1.
	double CatmullRomSpline(double point0, double point2, double intervalSize); // Returns tangent gradient at point 1.

	// Used for normalized interpolation interval - [0,1].
	double NCubicInterpolate();
	double NBicubicInterpolate();
	double NTricubicInterpolate();

	double NFiniteDifference(double point0, double point1); // Returns tangent gradient at point 0.
	double NCardinalSpline(double point0, double point2, double tension);   // Returns tangent gradient at point 1.
	double NCatmullRomSpline(double point0, double point2); // Returns tangent gradient at point 1.
}