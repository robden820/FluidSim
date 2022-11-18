#include "Camera.h"

#include <iostream>

Camera::Camera(glm::vec3 position, glm::vec3 worldUp, glm::vec3 forward, float yaw, float pitch)
{
	Position = position;
	WorldUp = worldUp;
	Forward = forward;

	Yaw = yaw;
	Pitch = pitch;

	MovementSpeed = SPEED;
	MouseSensitivity = SENSITIVITY;
	Zoom = ZOOM;

	UpdateVectors();
}

Camera::Camera(float posX, float posY, float posZ, float worldUpX, float worldUpY, float worldUpZ, float forwardX, float forwardY, float forwardZ, float yaw, float pitch)
{
	Position = glm::vec3(posX, posY, posZ);
	WorldUp = glm::vec3(worldUpX, worldUpY, worldUpZ);
	Forward = glm::vec3(forwardX, forwardY, forwardZ);

	Yaw = yaw;
	Pitch = pitch;

	MovementSpeed = SPEED;
	MouseSensitivity = SENSITIVITY;
	Zoom = ZOOM;

	UpdateVectors();
}

void Camera::ProcessKeyboard(Camera_Movement direction, float deltaTime)
{
	float velocity = MovementSpeed * deltaTime;
	if (direction == Camera_Movement::FORWARD)
		Position += Forward * velocity;
	if (direction == Camera_Movement::BACKWARD)
		Position -= Forward * velocity;
	if (direction == Camera_Movement::LEFT)
		Position -= Right * velocity;
	if (direction == Camera_Movement::RIGHT)
		Position += Right * velocity;
	if (direction == Camera_Movement::UP)
		Position += Up * velocity;
	if (direction == Camera_Movement::DOWN)
		Position -= Up * velocity;
}

void Camera::ProcessMouseMovement(float xOffset, float yOffset, GLboolean constrainPitch)
{
	xOffset *= MouseSensitivity;
	yOffset *= MouseSensitivity;

	Yaw += xOffset;
	Pitch += yOffset;

	//	std::cout << "YAW: " << Yaw << " || PITCH: " << Pitch << std::endl;
	//	std::cout << "xOff: " << xOffset << " || yOff: " << yOffset << std::endl;

		// make sure that when pitch is out of bounds, screen doesn't get flipped
	if (constrainPitch)
	{
		if (Pitch > 89.0f)
			Pitch = 89.0f;
		if (Pitch < -89.0f)
			Pitch = -89.0f;
	}

	// update Front, Right and Up Vectors using the updated Euler angles
	UpdateVectors();
}

void Camera::ProcessMouseScroll(float yoffset)
{
	Zoom -= (float)yoffset;
	if (Zoom < 1.0f)
		Zoom = 1.0f;
	if (Zoom > 45.0f)
		Zoom = 45.0f;
}

void Camera::UpdateVectors()
{
	// calculate the new Front vector
	glm::vec3 forward;
	forward.x = cos(glm::radians(Yaw)) * cos(glm::radians(Pitch));
	forward.y = sin(glm::radians(Pitch));
	forward.z = sin(glm::radians(Yaw)) * cos(glm::radians(Pitch));
	Forward = glm::normalize(forward);

	//	std::cout << "YAW: " << Yaw << " || PITCH: " << Pitch << std::endl;
	//	std::cout << "FORWARD: (" << Forward.x << ", " << Forward.y << ", " << Forward.z << ")" << std::endl;

		// also re-calculate the Right and Up vector
	Right = glm::normalize(glm::cross(Forward, WorldUp));  // normalize the vectors, because their length gets closer to 0 the more you look up or down which results in slower movement.
	Up = glm::normalize(glm::cross(Right, Forward));
}