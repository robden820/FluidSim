#pragma once

#include <glad/glad.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

class Shader
{
public:
	/**
	*	Constructor will read and build the shader.
	*	These are the filepaths to the source code for each shader.
	*/
	Shader(const char* vertexPath, const char* fragmentPath);
	~Shader();

	/* Use and activate the shader */
	void Use();

	/* Uniform utility functions */
	void SetBool(const std::string& name, bool value) const;
	void SetInt(const std::string& name, int value) const;
	void SetFloat(const std::string& name, float value) const;
	void SetVector(const std::string& name, glm::vec3 value) const;
	void SetMatrix(const std::string& name, glm::mat4 value) const;

	// program ID
	unsigned int programID;
};
