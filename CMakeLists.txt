cmake_minimum_required(VERSION 3.24)

project(MeshOffsetting
	VERSION 0.1.0
	LANGUAGES CXX
)

set(MESHOFFSET_SOURCES
	MeshOffsetting/main.cpp
	MeshOffsetting/src/SSVH.cpp
	MeshOffsetting/src/CMeshO.cpp
	MeshOffsetting/src/UniformGrid.cpp
	MeshOffsetting/src/NarrowBandSDF.cpp
)

set(MESHOFFSET_HEADERS
	MeshOffsetting/include/SSVH.h
	MeshOffsetting/include/base_type.h
	MeshOffsetting/include/CMeshO.h
	MeshOffsetting/include/Tolerance.h
	MeshOffsetting/include/UniformGrid.h
	MeshOffsetting/include/NarrowBandSDF.h
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
		"${CMAKE_CURRENT_SOURCE_DIR}/MeshOffsetting/include"
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
	
	include(FetchContent)
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

    target_link_libraries(meshoffset_vcglib
        INTERFACE vcglib
	)
endif()

target_link_libraries(MeshOffset
	PRIVATE MeshOffsetting::VCGLib
)
