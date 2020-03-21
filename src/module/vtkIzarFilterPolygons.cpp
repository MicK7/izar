
#include "vtkIzarFilterPolygons.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkMath.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"

vtkStandardNewMacro(vtkIzarFilterPolygons)

vtkIzarFilterPolygons::vtkIzarFilterPolygons():
	vtkPolyDataAlgorithm()
{
	IZAR_WARNING
}

vtkIzarFilterPolygons::~vtkIzarFilterPolygons()
{
}

void vtkIzarFilterPolygons::PrintSelf(ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarFilterPolygons\n";
}

int vtkIzarFilterPolygons::FillInputPortInformation(
	int vtkNotUsed(port), vtkInformation* info)
{
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	return 1;
}

int vtkIzarFilterPolygons::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
{
	vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
	if(!inInfo)
	{
		return 0;
	}
	vtkPolyData* inData = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!inData)
	{
		return 0;
	}
	
	vtkInformation* outInfo = outVector->GetInformationObject(0);
	if(!outInfo)
	{
		return 0;
	}
	vtkPolyData* outData = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!outData)
	{
		return 0;
	}
	outData->Initialize();
	vtkPoints* pts = inData->GetPoints();
	if(!pts)
	{
		vtkErrorMacro("The dataset does not have any point");
		return 0;
	}
	outData->SetPoints(pts);
	outData->SetPolys(inData->GetPolys());
	
	// Add field and point data as is
	outData->SetFieldData(inData->GetFieldData());
	for(int ii = 0 ; ii < inData->GetPointData()->GetNumberOfArrays() ; ii++)
	{
		outData->GetPointData()->AddArray(inData->GetPointData()->GetArray(ii));
	}
	
	// Add cell data
	int nVerts = inData->GetNumberOfVerts();
	int nLines = inData->GetNumberOfLines();
	int nPolys = inData->GetNumberOfPolys();
	int iBegin = nVerts + nLines;
	int iEnd = iBegin + nPolys;
	for(int ii = 0 ; ii < inData->GetCellData()->GetNumberOfArrays() ; ii++)
	{
		vtkDataArray* arOrig = inData->GetCellData()->GetArray(ii);
		vtkDataArray* arNew = arOrig->NewInstance();
		arNew->SetName(arOrig->GetName());
		arNew->SetNumberOfComponents(arOrig->GetNumberOfComponents());
		arNew->SetNumberOfTuples(nPolys);
		arNew->InsertTuples(0, nPolys, iBegin, arOrig);
		outData->GetCellData()->AddArray(arNew);
		arNew->Delete();
	}
	
	return 1;
}
	
	
