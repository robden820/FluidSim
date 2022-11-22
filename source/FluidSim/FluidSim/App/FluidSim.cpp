// FluidSim.cpp : Defines the entry point for the application.
//

#include "FluidSim.h"

#include <iostream>
#include <memory>
#include <chrono>

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
//	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetCursorPosCallback(window, mouse_callback);
	glfwSetScrollCallback(window, scroll_callback);

	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	/* ~~~~~ Here we can set up our camera ~~~~~ */
	Camera cam(glm::vec3(0.0f, 0.0f, 40.0f));
	camera = cam;

	/* ~~~~~ Set up our shaders ~~~~~ */
	Shader shader("../FluidSim/App/Shaders/vertex.txt", "../FluidSim/App/Shaders/fragment.txt");

	/* ~~~~~ Set those in the application ~~~~~ */
	gApplication->SetShader(std::make_shared<Shader>(shader));
	gApplication->SetCamera(std::make_shared<Camera>(camera));

	gApplication->SetScreenWidth((float)SCREEN_WIDTH);
	gApplication->SetScreenHeight((float)SCREEN_HEIGHT);

	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

	// We want to enable depth testing using the z-buffer.
	glEnable(GL_DEPTH_TEST);
	// Enable face culling. Do we want this? TO_DO
	//glEnable(GL_CULL_FACE);


	gApplication->Initialize();

	// Timing 
	float deltaTime = 0.0f;   // Time between previous and current frame
	float lastFrame = 0.0f;   // Time of the last frame.

	// The render loop. Checks if the widnow has been told to close.
	while (!glfwWindowShouldClose(window))
	{
		//processInput(window);
		// Any rendering commands go here.
		// ~~~ RENDERING ~~~ //

		// Specify our clear color.
		glClearColor(0.2f, 0.3, 0.3f, 1.0f);
		// Clears the specified buffer with the color we just set. We also want to clear our depth buffer.
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		float currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		// Render the application
		if (gApplication != 0)
		{
			gApplication->Update(deltaTime);
			// Do we need to do this every frame?
//			glBindVertexArray(VAO);
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
	gApplication->Shutdown();

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