cmake_minimum_required(VERSION 3.24)
project(src)

set(CMAKE_CXX_STANDARD 17)

# find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(HDF5 REQUIRED COMPONENTS C CXX HL REQUIRED)


link_directories( ${HDF5_HL_LIBRARIES} )
include_directories ( ${HDF5_INCLUDE_DIR} )

# link_directories( ${HDF5_LIBRARY_DIRS} )
# include_directories ( ${HDF5_HL_LIBRARIES} )

#set(src -lhdf5_hl -lhdf5 -lhdf5_cpp)

#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O2 -std=c++17 -Wall -g -lstdc++fs -lpthread -fopenmp -lhdf5_hl -lhdf5 -lhdf5_cpp")

#SET(MY_VARIABLE "option_value" CACHE STRING "../tst/dynamics/parfile.par")
set_target_properties(${TARGET_NAME} PROPERTIES VS_DEBUGGER_COMMAND_ARGUMENTS "../tst/dynamics/parfile.par" )

add_executable(src src/main.cpp src/composition.h src/model_ejecta.h src/old/ejecta_grb.h src/old/ejecta_radiation.h src/initial_data.h src/blastwave.h src/radiation.h)

target_link_libraries(${PROJECT_NAME} stdc++fs HDF5::HDF5)

#cmake_minimum_required(VERSION 3.24)
#project(PyBlastAfterglowMag)
#
#set(CMAKE_CXX_STANDARD 17)
#
#add_executable(PyBlastAfterglowMag main.cpp)
