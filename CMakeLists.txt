cmake_minimum_required(VERSION 3.24)

project(MeshOffsetting
	VERSION 0.1.0
	LANGUAGES CXX
)

include(FetchContent)

set(MESHOFFSET_SOURCES
	MeshOffsetting/main.cpp
	MeshOffsetting/src/SSVH.cpp
	MeshOffsetting/src/FCPWClosestPointQuery.cpp
	MeshOffsetting/src/CMeshO.cpp
	MeshOffsetting/src/UniformGrid.cpp
	MeshOffsetting/src/NarrowBandSDF.cpp
	MeshOffsetting/src/SlicerSDF.cpp
	MeshOffsetting/Test/Test.cpp
)

set(MESHOFFSET_HEADERS
	MeshOffsetting/include/SSVH.h
	MeshOffsetting/include/FCPWClosestPointQuery.h
	MeshOffsetting/include/base_type.h
	MeshOffsetting/include/CMeshO.h
	MeshOffsetting/include/Tolerance.h
	MeshOffsetting/include/UniformGrid.h
	MeshOffsetting/include/NarrowBandSDF.h
	MeshOffsetting/include/SlicerSDF.h
	MeshOffsetting/FCPW/FCPWPointTriangleDistance.h
	MeshOffsetting/Test/Test.h
)

add_executable(MeshOffset
	${MESHOFFSET_SOURCES}
	${MESHOFFSET_HEADERS}
)

source_group(
	TREE "${CMAKE_CURRENT_SOURCE_DIR}/MeshOffsetting"
	FILES
		${MESHOFFSET_SOURCES}
		${MESHOFFSET_HEADERS}
)

target_include_directories(MeshOffset
	PRIVATE
		"${CMAKE_CURRENT_SOURCE_DIR}/MeshOffsetting"
		"${CMAKE_CURRENT_SOURCE_DIR}/MeshOffsetting/include"
		"${CMAKE_CURRENT_SOURCE_DIR}/MeshOffsetting/FCPW"
)

target_compile_features(MeshOffset PRIVATE cxx_std_20)

target_compile_definitions(MeshOffset
	PRIVATE
		MESHLAB_SCALAR=float
)

set(VCGLIB_ROOT
""
CACHE PATH
"VCGLib source directory or installation prefix"
)

find_path(VCGLIB_INCLUDE_DIR
	NAMES vcg/complex/complex.h

	HINTS
		"${VCGLIB_ROOT}"
		ENV VCGLIB_ROOT
		"D:/vcglib"
		"${CMAKE_CURRENT_SOURCE_DIR}/third_party/vcglib"
		
	PATH_SUFFIXES
	include
	
	DOC "Directory containing the vcg folder"
)

add_library(meshoffset_vcglib INTERFACE)

add_library(
	MeshOffsetting::VCGLib
	ALIAS meshoffset_vcglib
)

if(VCGLIB_INCLUDE_DIR)
	message(STATUS "Using local VCGLib: ${VCGLIB_INCLUDE_DIR}")
	set(VCGLIB_PLY_SOURCE "${VCGLIB_INCLUDE_DIR}/wrap/ply/plylib.cpp")
	
	target_include_directories(meshoffset_vcglib
		SYSTEM INTERFACE
		"${VCGLIB_INCLUDE_DIR}"
	)
	
	if(EXISTS "${VCGLIB_INCLUDE_DIR}/eigenlib/Eigen/Core")
		target_include_directories(meshoffset_vcglib
		SYSTEM INTERFACE
			"${VCGLIB_INCLUDE_DIR}/eigenlib"
		)
	else()
		find_package(Eigen3 CONFIG REQUIRED)
		target_link_libraries(meshoffset_vcglib
		INTERFACE Eigen3::Eigen
		)
endif()

else()
	message(STATUS "Local VCGLib not found; downloading it")
	
	set(VCG_HEADER_ONLY ON CACHE BOOL "" FORCE)
	set(VCG_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
	set(VCG_ALLOW_BUNDLED_EIGEN ON CACHE BOOL "" FORCE)
	FetchContent_Declare(
        vcglib

        GIT_REPOSITORY
            https://github.com/cnr-isti-vclab/vcglib.git

        GIT_TAG
            2025.07

        GIT_SHALLOW
            TRUE
    )

    FetchContent_MakeAvailable(vcglib)
	set(VCGLIB_PLY_SOURCE "${vcglib_SOURCE_DIR}/wrap/ply/plylib.cpp")

    target_link_libraries(meshoffset_vcglib
        INTERFACE vcglib
	)
endif()

set(FCPW_USE_ENOKI OFF CACHE BOOL "" FORCE)
set(FCPW_USE_EIGHT_WIDE_BRANCHING OFF CACHE BOOL "" FORCE)
set(FCPW_ENABLE_GPU_SUPPORT OFF CACHE BOOL "" FORCE)
set(FCPW_BUILD_BINDINGS OFF CACHE BOOL "" FORCE)
set(FCPW_BUILD_DEMO OFF CACHE BOOL "" FORCE)
set(FCPW_BUILD_TESTS OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
	meshoffset_fcpw_dependency
	GIT_REPOSITORY https://github.com/rohan-sawhney/fcpw.git
	GIT_TAG v1.2.0
	GIT_SHALLOW TRUE
	GIT_SUBMODULES ""
	UPDATE_DISCONNECTED TRUE
)

FetchContent_MakeAvailable(meshoffset_fcpw_dependency)

target_sources(MeshOffset PRIVATE "${VCGLIB_PLY_SOURCE}")
source_group("Third Party" FILES "${VCGLIB_PLY_SOURCE}")

target_link_libraries(MeshOffset
	PRIVATE
		MeshOffsetting::VCGLib
		fcpw
)

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
	target_link_libraries(MeshOffset PRIVATE OpenMP::OpenMP_CXX)
endif()
