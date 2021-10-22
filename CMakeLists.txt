cmake_minimum_required(VERSION 3.8)
PROJECT( VascuSynth )

#FIND_PACKAGE( ITK )
#IF( ITK_FOUND )
#	MESSAGE(STATUS "  - Found ITK library")
#ELSE()
#	MESSAGE(ERROR "ITK library not found")
#ENDIF (ITK_FOUND )


MESSAGE(STATUS "Trying to find HDF5 libraries (minimum version 1.8.13) with C++ support...")
find_package(HDF5 1.8.13 COMPONENTS C CXX REQUIRED)
IF(HDF5_FOUND)
	MESSAGE(STATUS "  - Found HDF5 library")
	SET(MOD_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
	MESSAGE(STATUS "Include Directories: " ${HDF5_INCLUDE_DIRS})
	MESSAGE(STATUS "Libraries: " ${HDF5_LIBRARIES})
ENDIF()

SET(CMAKE_CXX_STANDARD 11)


ADD_EXECUTABLE(VascuSynth VascuSynth.cpp SupplyMap.cpp OxygenationMap.cpp NodeTable.cpp VascularTree.cpp TreeDrawer.cpp)

TARGET_LINK_LIBRARIES( VascuSynth ${HDF5_LIBRARIES})
#INCLUDE(${HDF5_INCLUDE_DIRS})
