#include "Sphere.h"

#define GOLDEN_RATIO 1.61803398875

float* Sphere::GetVertices()
{
	float vertices[] =
	{
		-0.5f, -0.5f, -0.5f, // left bottom back
		0.5f, -0.5f, -0.5f,  // right bottom back
		-0.5f, 0.5f, -0.5f,  // left top back
		0.5f, 0.5f, -0.5f,   // right top back
		-0.5f, -0.5f, 0.5f, // left bottom front
		0.5f, -0.5f, 0.5f,  // right bottom front
		-0.5f, 0.5f, 0.5f,  // left top front
		0.5f, 0.5f, 0.5f,   // right top front

	};

	return vertices;
}
