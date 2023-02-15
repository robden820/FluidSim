#include "Application.h"

#include "GLFW/glfw3.h"

Application::Application(Camera& inCamera, Shader& inShader)
	: mCamera(inCamera), mShader(inShader)
{
	mShader.Use();

	m3Dsimulation = false;

	ApplicationData newData;
	mApplicationData = newData;

	// Initialise simulation data.
	mApplicationData.SetDeltaTime(0.03);
	mApplicationData.SetFluidDensity(1000.0);
	mApplicationData.SetFLIPBlend(0.95);

	// Set MACGrid data
	mApplicationData.SetGridLeft(-10.0);
	mApplicationData.SetGridBottom(-10.0);

	mApplicationData.SetNumGridCellsWidth(70);
	mApplicationData.SetNumGridCellsHeight(80);
	mApplicationData.UpdateNumGridCells();

	mApplicationData.SetGridCellSize(0.2);
}

void Application::Initialize()
{
	if (m3Dsimulation)
	{
		double start = glfwGetTime();

		Fluid3D fluid(mApplicationData);
		mFluid = std::make_unique<Fluid3D>(fluid);

		std::cout << "Initializing fluid: " << glfwGetTime() - start << "\n";
		start = glfwGetTime();

		std::cout << "Initializing Draw fluid: " << glfwGetTime() - start << "\n";
		start = glfwGetTime();

		// Needs to be set to 10/MACGrid resolution
		// TO DO: fix all initialization values.
		VoxelFluid3D voxelFluid(mApplicationData);
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
		double start = glfwGetTime();

//		Fluid2D fluid(mApplicationData);
//		mFluid = std::make_unique<Fluid2D>(fluid);

		Simulation2D sim(mApplicationData);
		mSimulation = sim;

		std::cout << "Initializing fluid: " << glfwGetTime() - start << "\n";
		start = glfwGetTime();

		std::cout << "Initializing Draw fluid: " << glfwGetTime() - start << "\n";
		start = glfwGetTime();

		VoxelFluid2D voxelFluid(mApplicationData);
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
	// - The fifth argument is the stide of the data, and is the space between consecutive pieces of data. Since each vertex is separated by 3 doubles, we note that as the stride here.
	// - The final argument is the offset of where the data begins in the buffer, which for us is 0.

	// Deals with out position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

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

void Application::Update()
{
	//mFluid->Update(mApplicationData);
	mSimulation.StepSimulation(mApplicationData);
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

	double scale = mApplicationData.GetGridCellSize();
	
	if (m3Dsimulation)
	{
		for (int p = 0; p < mApplicationData.GetNumParticles(); p++)
		{
			glm::mat4 model = glm::mat4(1.0f);
			glm::vec3 particlePosition = mApplicationData.Get3DParticlePosition(p);

			model = glm::translate(model, particlePosition);
			model = glm::scale(model, glm::vec3(scale * 0.25f));

			mShader.SetMatrix("model", model);

			glDrawArrays(GL_TRIANGLES, 0, 36);
		}

		for (int v = 0; v < mVoxelFluid3D.GetVoxelCenters().size(); v++)
		{
			glm::mat4 model = glm::mat4(1.0f);
			model = glm::translate(model, mVoxelFluid3D.GetVoxelCenter(v));
			model = glm::scale(model, glm::vec3(scale));

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
		for (int p = 0; p < mApplicationData.GetNumParticles(); p++)
		{
			glm::mat4 model = glm::mat4(1.0f);

			glm::vec2 vec = mApplicationData.Get2DParticlePosition(p);
			glm::vec3 particlePosition = { vec.x, vec.y, 0.0f };

			model = glm::translate(model, particlePosition);
			model = glm::scale(model, glm::vec3(scale * 0.5f));

			mShader.SetMatrix("model", model);

			glDrawArrays(GL_TRIANGLES, 0, 36);
		}

		for (int v = 0; v < mVoxelFluid2D.GetVoxelCenters().size(); v++)
		{
			glm::mat4 model = glm::mat4(1.0f);	
			glm::vec3 vec = { mVoxelFluid2D.GetVoxelCenter(v).x, mVoxelFluid2D.GetVoxelCenter(v).y, 0.0f };

			model = glm::translate(model, vec);
			model = glm::scale(model, glm::vec3(scale));

			mShader.SetMatrix("model", model);

			if (mApplicationData.GetCellType(v) == CellType::eFLUID)
			{
				mShader.SetVector("color", glm::vec3(0.0f, 0.0f, 1.0f));
				mShader.SetFloat("alpha", 1.0f);
				glDrawArrays(GL_LINES, 0, 36);
			}
			else if (mApplicationData.GetCellType(v) == CellType::eSOLID)
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

