#include "Application.h"

void Application::Initialize()
{
	Fluid fluid(100);
	mFluid = fluid;

	DrawFluid drawFluid(std::make_shared<Fluid>(mFluid));
	mDrawFluid = drawFluid;

	Sphere s;
	mSphere = s;
	mSphere.SetRadius(0.1f);

	InitGLObjects();
}

void Application::InitGLObjects()
{
	// Generate a buffer object with a unique ID.
	// We also need to assign a vertex array object, and vertex attribute calls can be stored inside the VAO. This saves us having to make lots of calls the set up our vertex attributes.
	// This is used in conjunction with our VBO.
	
	glGenVertexArrays(1, &VAO); // Can generate multiple array and buffer objects at the same time.
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	// Bind our vertex array object first.
	glBindVertexArray(VAO);

	// Bind our buffer object to the array buffer type. This is the type needed for a vertex buffer object.
	glBindBuffer(GL_ARRAY_BUFFER, VBO);

	// Now any buffer calls we make to GL_ARRAY_BUFFER will be used with the currently bound buffer, which is our VBO.
	// This function copies our vertex data into the buffer.
	// GL_STATIC_DRAW means we are going to set the data once and use it multiple times.
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * mSphere.GetVertices().get()->size(), mSphere.GetVertices()->data(), GL_STATIC_DRAW);

	// Similarly for our element buffer object and our indices
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * mSphere.GetVertices().get()->size(), mSphere.GetIndices()->data(), GL_STATIC_DRAW);

	// We can specify any input we want to a vertex shader, so we need to tell openGl how to interpret our vertex data.
	// - The first parameter specifies which vertex attribute in the shader we want to configure. In the shader we used layout(location=0), so we pass in 0 to match that attribute.
	// - The second parameter is the size of our attribute. Each vertex is a vec3, so the size of the vertex attribute is 3.
	// - The third parameter is the type of data we are using.
	// - The fourth specifies if we want to normalize the data.
	// - The fifth argument is the stide of the data, and is the space between consecutive pieces of data. Since each vertex is separated by 3 floats, we note that as the stride here.
	// - The final argument is the offset of where the data begins in the buffer, which for us is 0.

	// Deals with out position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	// Deal with our color attribute
//	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0); // Our color data is offset by 3 floats.
//	glEnableVertexAttribArray(1);

	// We can now safely unbind our VAO and VBO.
	glBindBuffer(GL_ARRAY_BUFFER, 0);
//	glBindVertexArray(0);
	// - WARNING - we cannot unbind the EBO while the VAO is active as the bound EBO is stored in the VAO.
	//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	// Bind our vertex array object.
	//glBindVertexArray(VAO);
}

void Application::SetShader(std::shared_ptr<Shader> inShader)
{
	mShader = inShader;
	mShader.get()->Use();
}

void Application::Update(float deltaTime)
{
	mFluid.StepSimulation(deltaTime);

	mDrawFluid.FromFluid(std::make_shared<Fluid>(mFluid));
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
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &EBO);
}

