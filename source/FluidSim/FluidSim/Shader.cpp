#include "Shader.h"

Shader::Shader(const char* vertexPath, const char* fragmentPath)
{
	// Retrieve each shaders sourcecode from the provided filepaths.
	std::string vertexSourceCode;
	std::ifstream vertexShaderFile;

	std::string fragmentSourceCode;
	std::ifstream fragmentShaderFile;

	// Check our filestreams for any issues. Throw any exceptions.
	vertexShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	fragmentShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try
	{
		// Open each file.
		vertexShaderFile.open(vertexPath);
		fragmentShaderFile.open(fragmentPath);

		// Read each files buffer contents into a stream.
		std::stringstream vertexShaderStream;
		std::stringstream fragmentShaderStream;

		vertexShaderStream << vertexShaderFile.rdbuf();
		fragmentShaderStream << fragmentShaderFile.rdbuf();

		// Close each file.
		vertexShaderFile.close();
		fragmentShaderFile.close();

		// Convert the data from each stream into a string.
		vertexSourceCode = vertexShaderStream.str();
		fragmentSourceCode = fragmentShaderStream.str();
	}
	catch (std::ifstream::failure exp)
	{
		std::cout << "ERROR:: SHADER FILE NOT SUCCESSFULLY READ:\n" << exp.what() << std::endl;
	}

	const char* vertexShaderCode = vertexSourceCode.c_str();
	const char* fragmentShaderCode = fragmentSourceCode.c_str();

	// Compile the shaders.

	int success;
	char infoLog[512];

	// Handle our vertex shader first.

	// OpenGl needs to dynamically compile the shader code at run-time.
	// Create a shader object, which we can use to store the vertex shader.
	unsigned int vertexShader;
	vertexShader = glCreateShader(GL_VERTEX_SHADER);

	// Pass over our shader source code to the newly created vertex shader to be compiled.
	// This specifies we are passing one string to be used as our source code.
	glShaderSource(vertexShader, 1, &vertexShaderCode, NULL);
	glCompileShader(vertexShader);

	// check for shader compile errors
	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
		std::cout << "ERROR:: VERTEX SHADER COMPILATION_FAILED:\n" << infoLog << std::endl;
	}

	// Handle the fragment shader second.
		// Similarly for the fragment shader.
	unsigned int fragmentShader;
	fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShader, 1, &fragmentShaderCode, NULL);
	glCompileShader(fragmentShader);

	// check for shader compile errors
	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
		std::cout << "ERROR:: FRAGMENT SHADER COMPILATION_FAILED:\n" << infoLog << std::endl;
	}

	// We now need to link these shaders into a shader program.
	// We can have multiple different shader programs if we want to use multiple different shaders.
	programID = glCreateProgram();

	// Attach both of our shaders to the shader program. OpenGL will handle linking the output of the vertex shader to the input of the fragment shader.
	glAttachShader(programID, vertexShader);
	glAttachShader(programID, fragmentShader);
	glLinkProgram(programID);

	// check for linking errors
	glGetProgramiv(programID, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(programID, 512, NULL, infoLog);
		std::cout << "ERROR:: SHADER PROGRAM LINKING_FAILED:\n" << infoLog << std::endl;
	}

	// Once we have linked the shaders into the program, we no longer need them.
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
}

Shader::~Shader()
{
	glDeleteProgram(programID);
}

void Shader::Use()
{
	// We can then activate this shader program.
	glUseProgram(programID);
}

void Shader::SetBool(const std::string& name, bool value) const
{
	// Find our uniform in the shader program and set it's value.
	glUniform1i(glGetUniformLocation(programID, name.c_str()), (int)value);
}

void Shader::SetInt(const std::string& name, int value) const
{
	glUniform1i(glGetUniformLocation(programID, name.c_str()), value);
}

void Shader::SetFloat(const std::string& name, float value) const
{
	glUniform1f(glGetUniformLocation(programID, name.c_str()), value);
}

void Shader::SetMatrix(const std::string& name, glm::mat4 value) const
{
	glUniformMatrix4fv(glGetUniformLocation(programID, name.c_str()), 1, GL_FALSE, glm::value_ptr(value));   // Specify that we are sending a single matrix and that we aren't using the transpose.
}