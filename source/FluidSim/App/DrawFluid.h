#pragma once

#include <memory>

#include "Shader.h"
#include "Fluid.h"
#include "glm/vec3.hpp"

class DrawFluid
{
	public:
		DrawFluid() = default;
		~DrawFluid() = default;

		DrawFluid(std::shared_ptr<Fluid> inFluid);



		void DrawParticles();
		void DrawGrid();

	private:
		std::vector<glm::vec3> mParticlePoints;
		std::vector<glm::vec3> mGridPoints;

		Shader* mShader;
};

	Shader* mShader;
private:
	DebugDraw(const DebugDraw&);
	DebugDraw& operator=(const DebugDraw&);
public:
	DebugDraw();
	DebugDraw(unsigned int size);
	~DebugDraw();

	unsigned int Size();
	void Resize(unsigned int newSize);
	vec3& operator[](unsigned int index);
	void Push(const vec3& v);

	void FromPose(Pose& pose);

	void UpdateOpenGLBuffers();
	void Draw(DebugDrawMode mode, const vec3& color, const matrix4& mvp);