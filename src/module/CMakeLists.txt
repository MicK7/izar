
# define a warning option
#option(ENABLE_SUPPORT_WARNING "Check this to disable the warnings at the filters loading" ON)
#mark_as_advanced(FORCE ENABLE_SUPPORT_WARNING)
#if(${ENABLE_SUPPORT_WARNING})
#	add_definitions(-DPRINT_SUPPORT_WARNING)
#endif()

set(classes
	vtkIzarAddCpTPolynomial
	vtkIzarAddFieldDataByBlock
	vtkIzarAddFieldDataCpT
	vtkIzarAddFieldDataGeneric
	vtkIzarAddFieldData
	vtkIzarAddFieldDataPG
	vtkIzarAxialProjection
	vtkIzarChangeFrame
	vtkIzarChorochronicDuplication
	vtkIzarConeCut
	vtkIzarCpTReconstructFields
	vtkIzarCylindricalProjection
	vtkIzarExtractCylindricalComponents
	vtkIzarExtractSubsetFromFile
	vtkIzarFilterBlocksByFieldData
	vtkIzarFilterPolygons
	vtkIzarGenericXRCut
	vtkIzarReconstructFields
	vtkIzarRemoveEmptyBlocks
	vtkIzarRestoreCoordinates
	vtkIzarRotate
	vtkIzarRotateSteady
	vtkIzarRotationalDuplication
	vtkIzarRThetaCalculator
	vtkIzarXRCutFromCSV
)

set(SERVER_MANAGER_XML_files
	Izar.xml)

set(SERVER_MANAGER_SOURCES_files
	vtkIzarAddFieldData.cpp
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
)

set(SERVER_SOURCES_files
	IzarHelpers.cpp)

set(SERVER_MANAGER_HEADERS_files
	IzarDefines.h
	vtkIzarAddCpTPolynomial.h
	vtkIzarAddFieldDataByBlock.h
	vtkIzarAddFieldDataCpT.h
	vtkIzarAddFieldDataGeneric.h
	vtkIzarAddFieldData.h
	vtkIzarAddFieldDataPG.h
	vtkIzarAxialProjection.h
	vtkIzarChangeFrame.h
	vtkIzarChorochronicDuplication.h
	vtkIzarConeCut.h
	vtkIzarCpTReconstructFields.h
	vtkIzarCylindricalProjection.h
	vtkIzarExtractCylindricalComponents.h
	vtkIzarExtractSubsetFromFile.h
	vtkIzarFilterBlocksByFieldData.h
	vtkIzarFilterPolygons.h
	vtkIzarGenericXRCut.h
	vtkIzarReconstructFields.h
	vtkIzarRemoveEmptyBlocks.h
	vtkIzarRestoreCoordinates.h
	vtkIzarRotate.h
	vtkIzarRotateSteady.h
	vtkIzarRotationalDuplication.h
	vtkIzarRThetaCalculator.h
	vtkIzarXRCutFromCSV.h
)

set(private_headers
	IzarHelpers.h
)

set_source_files_properties(
IzarHelpers
WRAP_EXCLUDE
)

vtk_module_add_module(Izar::IzarModule
  SOURCES ${SERVER_MANAGER_SOURCES_files} ${SERVER_SOURCES_files}
  PRIVATE_HEADERS ${private_headers}
  HEADERS ${SERVER_MANAGER_HEADERS_files}
)

vtk_module_link(Izar::IzarModule
  PRIVATE
  ${LIBS}
)

paraview_add_server_manager_xmls(
  XMLS ${SERVER_MANAGER_XML_files}
)

