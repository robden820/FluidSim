#pragma once

#include <vector>
#include <memory>
#include "glm/vec3.hpp"

class Sphere
{
	public:
		Sphere();
		~Sphere() = default;

		const std::vector<float>& GetVertices() { return mVertices; }
		const std::vector<unsigned int>& GetIndices() {return mIndices; }

		void SetRadius(float inRadius);

		int mNumVertices;

	private:

		std::vector<float> mVertices;
		std::vector<unsigned int> mIndices;

		float mRadius;
};