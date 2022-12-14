#include "Application.h"

#include "GLFW/glfw3.h"

Application::Application(Camera& inCamera, Shader& inShader)
	: mCamera(inCamera), mShader(inShader)
{
	mShader.Use();

	m3Dsimulation = false;
}

void Application::Initialize()
{
	if (m3Dsimulation)
	{
		float start = glfwGetTime();

		Fluid3D fluid(1000);
		mFluid = std::make_unique<Fluid3D>(fluid);

		std::cout << "Initializing fluid: " << glfwGetTime() - start << "\n";
		start = glfwGetTime();

		DrawFluid3D drawFluid(dynamic_cast<Fluid3D&>(*mFluid.get()));
		mDrawFluid3D = drawFluid;

		std::cout << "Initializing Draw fluid: " << glfwGetTime() - start << "\n";
		start = glfwGetTime();

		// Needs to be set to 10/MACGrid resolution
		// TO DO: fix all initialization values.
		VoxelFluid3D voxelFluid(dynamic_cast<Fluid3D&>(*mFluid.get()));
		mVoxelFluid3D = voxelFluid;

		std::cout << "Initializing voxel grid: " << glfwGetTime() - start << "\n";

		Sphere s;
		mSphere = s;

		start = glfwGetTime();

		InitGLObjects();

		std::cout << "Initializing GL objects: " << glfwGetTime() - start << "\n";
	}
	else
	{
		float start = glfwGetTime();

		Fluid2D fluid(100);
		mFluid = std::make_unique<Fluid2D>(fluid);

		std::cout << "Initializing fluid: " << glfwGetTime() - start << "\n";
		start = glfwGetTime();

		DrawFluid2D drawFluid(dynamic_cast<Fluid2D&>(*mFluid.get()));
		mDrawFluid2D = drawFluid;

		std::cout << "Initializing Draw fluid: " << glfwGetTime() - start << "\n";
		start = glfwGetTime();

		// Needs to be set to 10/MACGrid resolution
		// TO DO: fix all initialization values.
		VoxelFluid2D voxelFluid(dynamic_cast<Fluid2D&>(*mFluid.get()));
		mVoxelFluid2D = voxelFluid;

		std::cout << "Initializing voxel grid: " << glfwGetTime() - start << "\n";


		Sphere s;
		mSphere = s;

		start = glfwGetTime();

		InitGLObjects();

		std::cout << "Initializing GL objects: " << glfwGetTime() - start << "\n";
	}
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
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * mSphere.GetVertices().size(), mSphere.GetVertices().data(), GL_STATIC_DRAW);

	// Similarly for our element buffer object and our indices
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * mSphere.GetVertices().size(), mSphere.GetIndices().data(), GL_STATIC_DRAW);

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

void Application::SetShader(Shader& inShader)
{
	mShader = inShader;
	mShader.Use();
}

void Application::Update(float deltaTime)
{
	mFluid->StepSimulation(deltaTime);

	if (m3Dsimulation)
	{
		mDrawFluid3D.FromFluid(dynamic_cast<Fluid3D&>(*mFluid.get()));
		mVoxelFluid3D.UpdateVoxelStates(dynamic_cast<Fluid3D&>(*mFluid.get()));
	}
	else
	{
		mDrawFluid2D.FromFluid(dynamic_cast<Fluid2D&>(*mFluid.get()));
		mVoxelFluid2D.UpdateVoxelStates(dynamic_cast<Fluid2D&>(*mFluid.get()));
	}
}

void Application::Render(float inAspectRatio)
{
	glm::mat4 projection = glm::perspective(glm::radians(mCamera.Zoom), mScreenWidth / mScreenHeight, 0.1f, 100.0f);
	glm::mat4 view = glm::lookAt(mCamera.Position, mCamera.Position + mCamera.Forward, mCamera.Up);
	glm::mat4 mvp = projection * view; // No model


	mShader.SetMatrix("view", view);
	mShader.SetMatrix("projection", projection);

	mShader.SetVector("color", glm::vec3(1.0f, 0.0f, 0.0f));
	mShader.SetFloat("alpha", 1.0f);
	
	if (m3Dsimulation)
	{
		for (int p = 0; p < mDrawFluid3D.mParticlePoints.size(); p++)
		{
			glm::mat4 model = glm::mat4(1.0f);
			model = glm::translate(model, mDrawFluid3D.mParticlePoints[p]);
			model = glm::scale(model, glm::vec3(0.1f, 0.1f, 0.1f));

			mShader.SetMatrix("model", model);

			glDrawArrays(GL_TRIANGLES, 0, 36);
		}

		for (int v = 0; v < mVoxelFluid3D.GetVoxelCenters().size(); v++)
		{
			glm::mat4 model = glm::mat4(1.0f);
			model = glm::translate(model, mVoxelFluid3D.GetVoxelCenter(v));
			model = glm::scale(model, glm::vec3(0.5f));

			mShader.SetMatrix("model", model);

			if (mVoxelFluid3D.GetVoxelState(v) == VoxelFluid::VoxelState::eFLUID)
			{
				mShader.SetVector("color", glm::vec3(0.0f, 0.0f, 1.0f));
				mShader.SetFloat("alpha", 1.0f);
				glDrawArrays(GL_LINES, 0, 36);
			}
			else if (mVoxelFluid3D.GetVoxelState(v) == VoxelFluid::VoxelState::eSOLID)
			{
				mShader.SetVector("color", glm::vec3(1.0f, 1.0f, 0.0f));
				mShader.SetFloat("alpha", 0.1f);
				glDrawArrays(GL_LINES, 0, 36);
			}
		}
	}
	else
	{
		for (int p = 0; p < mDrawFluid2D.mParticlePoints.size(); p++)
		{
			glm::mat4 model = glm::mat4(1.0f);

			glm::vec3 vec = { mDrawFluid2D.mParticlePoints[p].x, mDrawFluid2D.mParticlePoints[p].y, 0.0f };

			model = glm::translate(model, vec);
			model = glm::scale(model, glm::vec3(0.1f, 0.1f, 0.1f));

			mShader.SetMatrix("model", model);

			glDrawArrays(GL_TRIANGLES, 0, 36);
		}

		for (int v = 0; v < mVoxelFluid2D.GetVoxelCenters().size(); v++)
		{
			glm::mat4 model = glm::mat4(1.0f);
			
			glm::vec3 vec = { mVoxelFluid2D.GetVoxelCenter(v).x, mVoxelFluid2D.GetVoxelCenter(v).y, 0.0f };

			model = glm::translate(model, vec);
			model = glm::scale(model, glm::vec3(0.5f));

			mShader.SetMatrix("model", model);

			if (mVoxelFluid2D.GetVoxelState(v) == VoxelFluid::VoxelState::eFLUID)
			{
				mShader.SetVector("color", glm::vec3(0.0f, 0.0f, 1.0f));
				mShader.SetFloat("alpha", 1.0f);
				glDrawArrays(GL_LINES, 0, 36);
			}
			else if (mVoxelFluid2D.GetVoxelState(v) == VoxelFluid::VoxelState::eSOLID)
			{
				mShader.SetVector("color", glm::vec3(1.0f, 1.0f, 0.0f));
				mShader.SetFloat("alpha", 0.1f);
				glDrawArrays(GL_LINES, 0, 36);
			}
		}
	}
}

void Application::Shutdown()
{
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &EBO);
}

