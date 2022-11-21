// FluidSim.cpp : Defines the entry point for the application.
//

#include "FluidSim.h"

#include <iostream>
#include <memory>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Application.h"
#include "Camera.h"
#include "Shader.h"



void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void mouse_callback(GLFWwindow* window, double xPos, double yPos);
void scroll_callback(GLFWwindow* window, double xOffset, double yOffset);

// Screen settings
const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 600;

// Camera
Camera camera;
// Timing 
float deltaTime = 0.0f;   // Time between previous and current frame
float lastFrame = 0.0f;   // Time of the last frame.

// Mouse input
float lastX = SCREEN_WIDTH / 2;
float lastY = SCREEN_HEIGHT / 2;
bool firstMouse = true;

Application* gApplication = 0;


int main()
{

	gApplication = new Application();

	// Initialize the glfw library
	glfwInit();

	// Set window hints, these are effectively window settings.
	// Version major and minor specify the client API version that the context must be compatible with. We are using version 3.3.
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	// OpenGL profile specifies which profile we want to use, in this case the core profile.
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Create our window object.

	GLFWwindow* window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "LearnOpenGL", NULL, NULL);
	// If we failed to create the window, bail.
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}

	// Make the context of our window the main context.
	glfwMakeContextCurrent(window);

	// GLAD manages function pointers for OpenGL so we need to initialize it.
	// // We pass GLAD the function to load the address of the OpenGL function pointers.
	// If we fail to load GLAD, bail.
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}

	// Tell OpenGL the size of the rendering window so it knows how to display data with respect to the window.
	// First two parameters set the location of the lower left corner of the window, then its the width and height parameters.
	glViewport(0, 0, 800, 600);

	// Register our framesize function with GLFW so that when the window gets resized, so does the viewport.
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	// Tell glfw we want our cursor to be captured for input. We also want the scroll wheel as well.
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetCursorPosCallback(window, mouse_callback);
	glfwSetScrollCallback(window, scroll_callback);

	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	/* ~~~~~ Set up our triangle vertices ~~~~~ */

	float vertices[] = {
	-0.5f, -0.5f, -0.5f,  0.0f, 0.0f, 0.0f,
	 0.5f, -0.5f, -0.5f,  1.0f, 0.0f, 0.0f,
	 0.5f,  0.5f, -0.5f,  1.0f, 1.0f, 0.0f,
	 0.5f,  0.5f, -0.5f,  1.0f, 1.0f, 0.0f,
	-0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 0.0f,
	-0.5f, -0.5f, -0.5f,  0.0f, 0.0f, 0.0f,

	-0.5f, -0.5f,  0.5f,  0.0f, 0.0f, 1.0f,
	 0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 1.0f,
	 0.5f,  0.5f,  0.5f,  1.0f, 1.0f, 1.0f,
	 0.5f,  0.5f,  0.5f,  1.0f, 1.0f, 1.0f,
	-0.5f,  0.5f,  0.5f,  0.0f, 1.0f, 1.0f,
	-0.5f, -0.5f,  0.5f,  0.0f, 0.0f, 1.0f,

	-0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
	-0.5f,  0.5f, -0.5f,  1.0f, 1.0f, 0.0f,
	-0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 0.0f,
	-0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 0.0f,
	-0.5f, -0.5f,  0.5f,  0.0f, 0.0f, 0.0f,
	-0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 0.0f,

	 0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 1.0f,
	 0.5f,  0.5f, -0.5f,  1.0f, 1.0f, 1.0f,
	 0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 1.0f,
	 0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 1.0f,
	 0.5f, -0.5f,  0.5f,  0.0f, 0.0f, 1.0f,
	 0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 1.0f,

	-0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 0.0f,
	 0.5f, -0.5f, -0.5f,  1.0f, 1.0f, 0.0f,
	 0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
	 0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
	-0.5f, -0.5f,  0.5f,  0.0f, 0.0f, 0.0f,
	-0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 0.0f,

	-0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 1.0f,
	 0.5f,  0.5f, -0.5f,  1.0f, 1.0f, 1.0f,
	 0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 1.0f,
	 0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 1.0f,
	-0.5f,  0.5f,  0.5f,  0.0f, 0.0f, 1.0f,
	-0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 1.0f
	};

	float verticesWithColor[] = {
		// positions         // colors
		 0.5f, -0.5f, 0.0f,  1.0f, 0.0f, 0.0f,   // bottom right
		-0.5f, -0.5f, 0.0f,  0.0f, 1.0f, 0.0f,   // bottom left
		 0.0f,  0.5f, 0.0f,  0.0f, 0.0f, 1.0f    // top 
	};

	unsigned int indices[] = {
		0, 1, 2,   // Triangle one
		0, 2, 3    // Triangle two
	};

	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	/* ~~~~~ Here we can set up our camera ~~~~~ */
	Camera cam(glm::vec3(0.0f, 0.0f, 3.0f));
	camera = cam;

	/* ~~~~~ Set up our shaders ~~~~~ */
	Shader shader("../FluidSim/App/Shaders/vertex.txt", "../FluidSim/App/Shaders/fragment.txt");

	/* ~~~~~ Set those in the application ~~~~~ */
	gApplication->SetShader(std::make_shared<Shader>(shader));
	gApplication->SetCamera(std::make_shared<Camera>(camera));

	gApplication->SetScreenWidth((float)SCREEN_WIDTH);
	gApplication->SetScreenHeight((float)SCREEN_HEIGHT);

	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	/* ~~~~~ Set up our vertex array and buffer objects ~~~~~ */

	// Generate a buffer object with a unique ID.
	// We also need to assign a vertex array object, and vertex attribute calls can be stored inside the VAO. This saves us having to make lots of calls the set up our vertex attributes.
	// This is used in conjunction with our VBO.
	unsigned int VBO;
	unsigned int VAO;
	unsigned int EBO;
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
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	// Similarly for our element buffer object and our indices
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	// We can specify any input we want to a vertex shader, so we need to tell openGl how to interpret our vertex data.
	// - The first parameter specifies which vertex attribute in the shader we want to configure. In the shader we used layout(location=0), so we pass in 0 to match that attribute.
	// - The second parameter is the size of our attribute. Each vertex is a vec3, so the size of the vertex attribute is 3.
	// - The third parameter is the type of data we are using.
	// - The fourth specifies if we want to normalize the data.
	// - The fifth argument is the stide of the data, and is the space between consecutive pieces of data. Since each vertex is separated by 3 floats, we note that as the stride here.
	// - The final argument is the offset of where the data begins in the buffer, which for us is 0.

	// Deals with out position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	// Deal with our color attribute
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float))); // Our color data is offset by 3 floats.
	glEnableVertexAttribArray(1);

	// We can now safely unbind our VAO and VBO.
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	// - WARNING - we cannot unbind the EBO while the VAO is active as the bound EBO is stored in the VAO.
	//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	// We want to enable depth testing using the z-buffer.
	glEnable(GL_DEPTH_TEST);
	// Enable face culling. Do we want this? TO_DO
	//glEnable(GL_CULL_FACE);


	gApplication->Initialize();

	// The render loop. Checks if the widnow has been told to close.
	while (!glfwWindowShouldClose(window))
	{
		processInput(window);
		// Any rendering commands go here.
		// ~~~ RENDERING ~~~ //

		// Specify our clear color.
		glClearColor(0.2f, 0.3, 0.3f, 1.0f);
		// Clears the specified buffer with the color we just set. We also want to clear our depth buffer.
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Bind our vertex array object.
		glBindVertexArray(VAO);

		// Render the application
		if (gApplication != 0)
		{
			// Do we need to do this every frame?
			glBindVertexArray(VAO);
			float aspect = (float)SCREEN_WIDTH / (float)SCREEN_HEIGHT;
			gApplication->Render(aspect);
		}

		// ~~~~~~~~~~~~~~~~~ //

		// Swaps the color buffer that we use to render to the output buffer.
		// See Double Buffering.
		glfwSwapBuffers(window);

		// Check for any events like mouse inputs.
		glfwPollEvents();
	}

	// De-allocate all of our resources.
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &EBO);

	// Properly clean up all of GLFW's allocated resources.
	glfwTerminate();

	// Clear application resources.
	if (gApplication != 0)
	{
		std::cout << "Expected application to be null on exit\n";
		delete gApplication;
	}
	return 0;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void processInput(GLFWwindow* window)
{
	// Check if we are pressing the escape key
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
	{
		// Alert glfw that we should close the window.
		glfwSetWindowShouldClose(window, true);
	}

	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		camera.ProcessKeyboard(Camera_Movement::FORWARD, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		camera.ProcessKeyboard(Camera_Movement::BACKWARD, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		camera.ProcessKeyboard(Camera_Movement::LEFT, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		camera.ProcessKeyboard(Camera_Movement::RIGHT, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		camera.ProcessKeyboard(Camera_Movement::UP, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
		camera.ProcessKeyboard(Camera_Movement::DOWN, deltaTime);
}

void mouse_callback(GLFWwindow* window, double xPosIn, double yPosIn)
{
	float xpos = static_cast<float>(xPosIn);
	float ypos = static_cast<float>(yPosIn);

	if (firstMouse)
	{
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	float xoffset = xpos - lastX;
	float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

	lastX = xpos;
	lastY = ypos;

	camera.ProcessMouseMovement(xoffset, yoffset);
}

void scroll_callback(GLFWwindow* window, double xOffset, double yOffset)
{
	camera.ProcessMouseScroll(static_cast<float>(yOffset));
}