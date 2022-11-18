add_subdirectory(App)

target_sources(FluidSim PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/glad/glad.c")

target_include_directories(FluidSim PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/")

target_link_libraries(FluidSim "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/GLFW/glfw3.lib")
