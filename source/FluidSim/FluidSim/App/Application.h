#pragma once

#include <memory>

#include "Camera.h"
#include "Shader.h"
#include "VoxelFluid3D.h"
#include "VoxelFluid2D.h"

#include "Fluid3D.h"
#include "Fluid2D.h"

#include "Simulation2D.h"

#include "ApplicationData.h"

#include "Sphere.h"

class Application
{
public:
	Application(Camera& inCamera, Shader& inShader);
	~Application() = default;

	void Initialize();
	void Update();
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

	ApplicationData mApplicationData;

	float mScreenWidth;
	float mScreenHeight;

	bool m3Dsimulation;

	VoxelFluid2D mVoxelFluid2D;
	VoxelFluid3D mVoxelFluid3D;

	std::unique_ptr<Fluid> mFluid;

	Simulation2D* mSimulation;

	Sphere mSphere; // Contains vertex & index information to draw particles.

	unsigned int VBO;
	unsigned int VAO;
	unsigned int EBO;
};
