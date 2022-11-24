#pragma once

#include <memory>

#include "Camera.h"
#include "Shader.h"
#include "DrawFluid.h"
#include "VoxelFluid.h"

#include "Sphere.h"

class Application
{
public:
	Application(Camera& inCamera, Shader& inShader);
	~Application() = default;

	void Initialize();
	void Update(float deltaTime);
	void Render(float aspectRatio);
	void Shutdown();

	void SetCamera(Camera& inCamera) { mCamera = inCamera; }
	void SetShader(Shader& inShader);

	void SetScreenWidth(float inWidth) { mScreenWidth = inWidth; }
	void SetScreenHeight(float inHeight) { mScreenHeight = inHeight; }

private:

	void InitGLObjects();

	Camera& mCamera;
	Shader& mShader;

	float mScreenWidth;
	float mScreenHeight;

	DrawFluid mDrawFluid;
	VoxelFluid mVoxelFluid;

	Fluid mFluid;

	Sphere mSphere; // Contains vertex & index information to draw particles.

	unsigned int VBO;
	unsigned int VAO;
	unsigned int EBO;
};
