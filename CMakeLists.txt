cmake_minimum_required(VERSION 3.8)
PROJECT( VascuSynth )

FIND_PACKAGE( ITK )
IF( ITK_FOUND )
	MESSAGE(STATUS "  - Found ITK library")
ELSE()
	MESSAGE(ERROR "ITK library not found")
ENDIF (ITK_FOUND )


MESSAGE(STATUS "Trying to find HDF5 libraries (minimum version 1.8.13) with C++ support...")
IF(${VRN_USE_HDF5_VERSION} MATCHES "1.8")

	MESSAGE(STATUS "Trying to find HDF5 libraries (in between versions 1.8.13 and 1.10.0) with C++ support...")

	# First: Try to find the needed headers and libraries in some hard coded locations corresponding to
	# a definite 1.8 library. This is useful if HDF5 version 1.10 is the default on a system, but
	# version 1.8 can be installed in a non-standard location. This is similar to the procedure for ffmpeg
	# packages (at least on Arch Linux).
	FIND_PATH(HDF5_INCLUDE_DIR_MANUAL H5Cpp.h
		PATHS
		/usr/include/hdf5_18/
		/usr/local/include/hdf5_18/
		/opt/local/include/hdf5_18/
		NO_DEFAULT_PATH
	)

	FIND_LIBRARY(HDF5_LIB_C_MANUAL
		NAMES hdf5
		PATHS
		/usr/lib/hdf5_18/
		/usr/local/lib/hdf5_18/
		/opt/local/lib/hdf5_18/
		NO_DEFAULT_PATH
	)

	FIND_LIBRARY(HDF5_LIB_CPP_MANUAL
		NAMES hdf5_cpp
		PATHS
		/usr/lib/hdf5_18/
		/usr/local/lib/hdf5_18/
		/opt/local/lib/hdf5_18/
		NO_DEFAULT_PATH
	)

	IF(HDF5_INCLUDE_DIR_MANUAL AND HDF5_LIB_C_MANUAL AND HDF5_LIB_CPP_MANUAL)
		SET(HDF5_LIBRARIES
			${HDF5_LIB_C_MANUAL}
			${HDF5_LIB_CPP_MANUAL}
		)
		SET(HDF5_INCLUDE_DIRS
			${HDF5_INCLUDE_DIR_MANUAL}
		)

		SET(HDF5_VERSION 1.8)
		SET(HDF5_FOUND TRUE)
	ENDIF()
ENDIF()

IF(HDF5_FOUND)
	MESSAGE(STATUS "  - Found HDF5 library")
	IF(HDF5_VERSION VERSION_LESS 1.10.0)
		SET(MOD_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
		MESSAGE(STATUS "Include Directories: " ${HDF5_INCLUDE_DIRS})
		MESSAGE(STATUS "Libraries: " ${HDF5_LIBRARIES})
	ELSE()
		MESSAGE(FATAL_ERROR "HDF5 version is newer than 1.8 and may cause incompatibilites with other applications and installations. Disable the HDF5-Module or downgrade your HDF5 installation.")
	ENDIF()
ENDIF()


ADD_EXECUTABLE(VascuSynth VascuSynth.cpp SupplyMap.cpp OxygenationMap.cpp NodeTable.cpp VascularTree.cpp TreeDrawer.cpp)

TARGET_LINK_LIBRARIES( VascuSynth ${ITK_LIBRARIES} ${HDF5_LIBRARIES})
INCLUDE( ${ITK_USE_FILE} ${HDF5_INCLUDE_DIRS})
