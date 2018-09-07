cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

if (NOT ParaView_BINARY_DIR)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif()

include(ParaViewPlugins)

# define a warning option
option(ENABLE_SUPPORT_WARNING "Check this to disable the warnings at the filters loading" ON)
mark_as_advanced(FORCE ENABLE_SUPPORT_WARNING)
if(${ENABLE_SUPPORT_WARNING})
	add_definitions(-DPRINT_SUPPORT_WARNING)
endif()

find_package(yaml-cpp REQUIRED)
find_package(Boost REQUIRED)

include_directories(${YAML_CPP_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})
set(LIBS ${LIBS} ${YAML_CPP_LIBRARIES} ${Boost_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

# create a paraview plugin containing server manager xml and the server
# manager classes to build
# this plugin can be loaded on the server side
set_source_files_properties(
IzarHelpers
WRAP_EXCLUDE
)

ADD_PARAVIEW_PLUGIN(Izar "0.1"
  SERVER_MANAGER_XML Izar.xml
  SERVER_MANAGER_SOURCES 	vtkIzarAddFieldData.cpp
				vtkIzarAddFieldDataPG.cpp
				vtkIzarAddFieldDataCpT.cpp
				vtkIzarAddFieldDataGeneric.cpp
				vtkIzarAddFieldDataByBlock.cpp
				vtkIzarAddCpTPolynomial.cpp
				vtkIzarExtractSubsetFromFile.cpp
				vtkIzarFilterBlocksByFieldData.cpp
				vtkIzarReconstructFields.cpp
				vtkIzarCpTReconstructFields.cpp
				vtkIzarChangeFrame.cpp
				vtkIzarGenericXRCut.cpp
				vtkIzarConeCut.cpp
				vtkIzarXRCutFromCSV.cpp
				vtkIzarRotateSteady.cpp
				vtkIzarRotationalDuplication.cpp
				vtkIzarChorochronicDuplication.cpp
				vtkIzarRotate.cpp
				vtkIzarCylindricalProjection.cpp
				vtkIzarAxialProjection.cpp
				vtkIzarRestoreCoordinates.cpp
				vtkIzarRemoveEmptyBlocks.cpp
				vtkIzarExtractCylindricalComponents.cpp
				vtkIzarFilterPolygons.cpp
				vtkIzarRThetaCalculator.cpp			
 SERVER_SOURCES			IzarHelpers.cpp)

set_property(TARGET Izar PROPERTY CXX_STANDARD 11)
target_link_libraries(Izar LINK_PRIVATE ${LIBS})
