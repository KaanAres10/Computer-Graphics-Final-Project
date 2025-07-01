cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
project(ComputerGraphics C CXX)
enable_testing()

#option(USE_PREBUILT_INTERSECT "Enable using prebuilt intersection library" ON)
option(USE_PREBUILT "Enable using prebuilt libraries" ON)

if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/framework")
	# Create framework library and include CMake scripts (compiler warnings, sanitizers and static analyzers).
	add_subdirectory("framework")
	set(REFERENCE_MODE OFF)
elseif(PROJECT_IS_TOP_LEVEL)
	add_subdirectory("../../../framework/" "${CMAKE_CURRENT_BINARY_DIR}/framework/")
else()
	# During development the framework lives in parent folder.
	set(REFERENCE_MODE ON)
endif()

find_package(OpenGL REQUIRED)
find_package(OpenMP REQUIRED)

add_library(Bachelor_FinalProjectLib
	"src/scene.cpp"
	"src/draw.cpp"
	"src/screen.cpp"
	"src/light.cpp"
	"src/config.cpp"
	"src/texture.cpp"
	"src/shading.cpp"
	"src/interpolate.cpp"
	"src/recursive.cpp"
	"src/render.cpp"
	"src/extra.cpp"
	"src/verification.cpp"
	"src/bvh.cpp"
	"src/splines.cpp")

target_include_directories(Bachelor_FinalProjectLib PUBLIC "src")
target_link_libraries(Bachelor_FinalProjectLib PUBLIC CGFramework OpenGL::GLU OpenMP::OpenMP_CXX)
target_compile_features(Bachelor_FinalProjectLib PUBLIC cxx_std_20)
enable_sanitizers(Bachelor_FinalProjectLib)
set_project_warnings(Bachelor_FinalProjectLib)

if (USE_PREBUILT)
	add_library(Bachelor_Intersect_Import STATIC IMPORTED)
	if (WIN32)
		set_target_properties(Bachelor_Intersect_Import PROPERTIES IMPORTED_LOCATION_DEBUG "${CMAKE_CURRENT_SOURCE_DIR}/prebuilt/Provided_x64_debug.lib")
		set_target_properties(Bachelor_Intersect_Import PROPERTIES IMPORTED_LOCATION_RELWITHDEBINFO "${CMAKE_CURRENT_SOURCE_DIR}/prebuilt/Provided_x64_relwithdebinfo.lib")
		set_target_properties(Bachelor_Intersect_Import PROPERTIES IMPORTED_LOCATION_RELEASE "${CMAKE_CURRENT_SOURCE_DIR}/prebuilt/Provided_x64_release.lib")
	elseif (APPLE)
		set_target_properties(Bachelor_Intersect_Import PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/prebuilt/libProvided_mac_arm64.a")
	else()
		set_target_properties(Bachelor_Intersect_Import PROPERTIES IMPORTED_LOCATION_DEBUG "${CMAKE_CURRENT_SOURCE_DIR}/prebuilt/libProvided_linux_x64_debug.a")
		set_target_properties(Bachelor_Intersect_Import PROPERTIES IMPORTED_LOCATION_RELEASE "${CMAKE_CURRENT_SOURCE_DIR}/prebuilt/libProvided_linux_x64_release.a")
		set_target_properties(Bachelor_Intersect_Import PROPERTIES IMPORTED_LOCATION_RELWITHDEBINFO "${CMAKE_CURRENT_SOURCE_DIR}/prebuilt/libProvided_linux_x64_relwithdebinfo.a")
	endif()


#	set_target_properties(Bachelor_Intersect_Import PROPERTIES IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake-build-debug/Bachelor/final_project2_ray_tracing/CodeProvided/libProvided.a")
	target_link_libraries(Bachelor_FinalProjectLib PUBLIC Bachelor_Intersect_Import)
else()
	target_sources(Bachelor_FinalProjectLib PRIVATE "src/intersect.cpp")
endif()

target_compile_definitions(Bachelor_FinalProjectLib PUBLIC
	"-DDATA_DIR=\"${CMAKE_CURRENT_LIST_DIR}/data/\"")

add_executable(Bachelor_FinalProject "src/main.cpp")
target_link_libraries(Bachelor_FinalProject PUBLIC Bachelor_FinalProjectLib)
target_compile_features(Bachelor_FinalProject PUBLIC cxx_std_20)
enable_sanitizers(Bachelor_FinalProject)
set_project_warnings(Bachelor_FinalProject)

# Add test files
if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/tests/")
	add_subdirectory("tests")
endif()
if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/grading_tests/")
	add_subdirectory("grading_tests")
endif()