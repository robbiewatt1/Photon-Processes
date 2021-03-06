cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(PhotonProcess)

if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. "
	    		"Remove CMakeCache.txt and CMakeFiles and try again.")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

find_package(Geant4 REQUIRED ui_all)
find_package(Eigen3 3.0.1 REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)

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
file(GLOB data  ${CMAKE_CURRENT_SOURCE_DIR}/DataTables/*.h5)
file(GLOB bin  ${CMAKE_CURRENT_SOURCE_DIR}/Photon-Process.sh)


include(${Geant4_USE_FILE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Source/Tools
		    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Processes
		    ${CMAKE_CURRENT_SOURCE_DIR}/Source/Fields
		    ${CMAKE_CURRENT_SOURCE_DIR}/Source/GaussianProcess)

if(BUILD_GP)
	add_compile_definitions(USEGP)
	file(GLOB source_gp ${CMAKE_CURRENT_SOURCE_DIR}/Source/GaussianProcess/*.cpp)
	file(GLOB header_gp ${CMAKE_CURRENT_SOURCE_DIR}/Source/GaussianProcess/*.hh)
endif()

add_library(PhotonProcess SHARED ${source} ${header} ${source_gp} ${header_gp})
target_include_directories(PhotonProcess PUBLIC $<INSTALL_INTERFACE:include/PhotonProcess>
						${libgp_INCLUDE_DIRS}
						${Geant4_INCLUDE_DIRS}
						${HDF5_INCLUDE_DIRS})
target_link_libraries(PhotonProcess PUBLIC Eigen3::Eigen 
					   ${libgp_LIBRARIES}
					   ${HDF5_LIBRARIES}
					   ${Geant4_LIBRARIES})
install(TARGETS PhotonProcess EXPORT PhotonProcess-targets DESTINATION lib/PhotonProcess)
install(EXPORT PhotonProcess-targets DESTINATION share/PhotonProcess/cmake)
install(FILES ${header} ${header_gp}  DESTINATION include/PhotonProcess)

include(CMakePackageConfigHelpers)
configure_package_config_file(
	"PhotonProcessConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/PhotonProcessConfig.cmake"
	INSTALL_DESTINATION "share/PhotonProcess/cmake"
)
if(BUILD_GP)
	File(APPEND ${CMAKE_CURRENT_BINARY_DIR}/PhotonProcessConfig.cmake "add_definitions(-DUSEGP)")
endif()
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/PhotonProcessConfig.cmake" DESTINATION share/PhotonProcess/cmake)
install(FILES ${data} DESTINATION share/PhotonProcess/DataTables)
install(FILES ${bin} DESTINATION bin)
