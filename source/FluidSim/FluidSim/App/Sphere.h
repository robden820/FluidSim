#pragma once

#include <vector>
#include <memory>
#include "glm/vec3.hpp"

class Sphere
{
	public:
		Sphere();
		~Sphere() = default;

		std::shared_ptr<std::vector<float>> GetVertices() { return std::make_shared<std::vector<float>>(mVertices); }
		std::shared_ptr<std::vector<unsigned int>> GetIndices() {return std::make_shared<std::vector<unsigned int>>(mIndices); }

		void SetRadius(float inRadius);

		int mNumVertices;

	private:

		std::vector<float> mVertices;
		std::vector<unsigned int> mIndices;

		float mRadius;
};