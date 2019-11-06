cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(PhotonProject)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall")
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

find_package(Geant4 REQUIRED)
find_package(Eigen3 3.0.1 REQUIRED)

option(BUILD_HDF5 "Use HDF5" ON)
if(BUILD_GP)
	add_definitions(-DUSEHDF5)
	find_package(HDF5 REQUIRED COMPONENTS C CXX)
endif()

option(BUILD_GP "Use GP" OFF)
if(BUILD_GP)
	add_definitions(-DUSEGP)
	find_package(PkgConfig REQUIRED)
	pkg_check_modules(libgp REQUIRED libgp)
endif()


file(GLOB source ${CMAKE_CURRENT_SOURCE_DIR}/Source/Tools/*.cpp
	   	 ${CMAKE_CURRENT_SOURCE_DIR}/Source/Processes/*.cpp
		 ${CMAKE_CURRENT_SOURCE_DIR}/Source/Fields/*.cpp)
file(GLOB header ${CMAKE_CURRENT_SOURCE_DIR}/Source/Tools/*.hh
           	 ${CMAKE_CURRENT_SOURCE_DIR}/Source/Processes/*.hh
		 ${CMAKE_CURRENT_SOURCE_DIR}/Source/Fields/*.hh)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Source/Tools
		    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Processes
		    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Fields
		    ${CMAKE_CURRENT_SOURCE_DIR}/Source/GaussianProcess)

# Add the base required libraries
if(BUILD_GP)
	add_compile_definitions(USEGP)
	file(GLOB source_gp ${CMAKE_CURRENT_SOURCE_DIR}/Source/GaussianProcess/*.cpp)
	file(GLOB header_gp ${CMAKE_CURRENT_SOURCE_DIR}/Source/GaussianProcess/*.hh)
	add_library(PhotonProcess ${source} ${header} ${source_gp} ${header_gp})
	target_include_directories(PhotonProcess PUBLIC ${Geant4_INCLUDE_DIRS} ${libgp_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
	target_link_libraries(PhotonProcess ${Geant4_LIBRARIES} Eigen3::Eigen ${libgp_LIBRARIES} ${HDF5_LIBRARIES})
	add_executable(test test.cpp ${source} ${header} ${source_gp} ${header_gp})
	target_include_directories(test PUBLIC ${Geant4_INCLUDE_DIRS} ${libgp_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
        target_link_libraries(test ${Geant4_LIBRARIES} Eigen3::Eigen ${libgp_LIBRARIES} ${HDF5_LIBRARIES})
else()
	add_library(PhotonProcess ${source} ${header})
        target_include_directories(PhotonProcess PUBLIC ${libgp_INCLUDE_DIRS})
        target_link_libraries(PhotonProcess ${Geant4_LIBRARIES} Eigen3::Eigen)
endif()
