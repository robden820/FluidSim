#pragma once

#include <memory>

#include "Camera.h"
#include "Shader.h"
#include "DrawFluid.h"

#include "Sphere.h"

class Application
{
public:
	Application() = default;
	~Application() = default;

	void Initialize();
	void Update(float deltaTime);
	void Render(float aspectRatio);
	void Shutdown();

	void SetCamera(std::shared_ptr<Camera> inCamera) { mCamera = inCamera; }
	void SetShader(std::shared_ptr<Shader> inShader);

	void SetScreenWidth(float inWidth) { mScreenWidth = inWidth; }
	void SetScreenHeight(float inHeight) { mScreenHeight = inHeight; }

private:

	void InitGLObjects();

	std::shared_ptr<Camera> mCamera;
	std::shared_ptr<Shader> mShader;

	float mScreenWidth;
	float mScreenHeight;

	DrawFluid mDrawFluid;

	Fluid mFluid;

	Sphere mSphere; // Contains vertex & index information to draw particles.

	unsigned int VBO;
	unsigned int VAO;
	unsigned int EBO;
};
