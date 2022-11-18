# CMakeList.txt : CMake project for FluidSim, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
FILE(GLOB SOURCEFILES *.cpp)
FILE(GLOB HEADERFILES *.h)

add_executable(FluidSim ${HEADERFILES} ${SOURCEFILES})

target_sources(FluidSim PRIVATE "ThirdParty/glad/glad.c")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET FluidSim PROPERTY CXX_STANDARD 20)
endif()


include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/")

target_link_libraries(FluidSim "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/GLFW/glfw3.lib")
# TODO: Add tests and install targets if needed.
