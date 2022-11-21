#include "Application.h"

void Application::Initialize()
{
	Fluid fluid(10, 1);
	mFluid = fluid;

	DrawFluid drawFluid(std::make_shared<Fluid>(mFluid));
	mDrawFluid = drawFluid;
}

void Application::SetShader(std::shared_ptr<Shader> inShader)
{
	mShader = inShader;
	mShader.get()->Use();
}

void Application::Update(float deltaTime)
{
}

void Application::Render(float inAspectRatio)
{
	glm::mat4 projection = glm::perspective(glm::radians(mCamera.get()->Zoom), mScreenWidth / mScreenHeight, 0.1f, 100.0f);
	glm::mat4 view = glm::lookAt(mCamera.get()->Position, mCamera.get()->Position + mCamera.get()->Forward, mCamera.get()->Up);
	glm::mat4 mvp = projection * view; // No model


	mShader.get()->SetMatrix("view", view);
	mShader.get()->SetMatrix("projection", projection);

	for (int p = 0; p < mDrawFluid.mParticlePoints.size(); p++)
	{
		glm::mat4 model = glm::mat4(1.0f);
		model = glm::translate(model, mDrawFluid.mParticlePoints[p]);

		mShader->SetMatrix("model", model);

		glDrawArrays(GL_TRIANGLES, 0, 36);
	}
}

void Application::Shutdown()
{
}

