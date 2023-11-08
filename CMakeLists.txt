cmake_minimum_required(VERSION 3.0)

project(lqr)

# Set the C++ standard (you can adjust this to your preferred version)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define the source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Define the header files
file(GLOB_RECURSE HEADERS "include/*.h")

# Add the executable
add_executable(lqr ${SOURCES} ${HEADERS})

# Include directories for header files
target_include_directories(lqr PRIVATE "include")
