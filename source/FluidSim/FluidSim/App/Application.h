#pragma once

#include <memory>

#include "Camera.h"
#include "Shader.h"
#include "DrawFluid3D.h"
#include "DrawFluid2D.h"
#include "VoxelFluid3D.h"
#include "VoxelFluid2D.h"

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

	bool m3Dsimulation;

	DrawFluid3D mDrawFluid3D;
	VoxelFluid3D mVoxelFluid3D;

	DrawFluid2D mDrawFluid2D;
	VoxelFluid2D mVoxelFluid2D;

	std::unique_ptr<Fluid> mFluid;

	Sphere mSphere; // Contains vertex & index information to draw particles.

	unsigned int VBO;
	unsigned int VAO;
	unsigned int EBO;
};
