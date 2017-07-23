
#include "vtkIzarAxialProjection.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"

vtkStandardNewMacro(vtkIzarAxialProjection)

vtkIzarAxialProjection::vtkIzarAxialProjection():
	vtkPointSetAlgorithm()
{
	IZAR_WARNING
}

vtkIzarAxialProjection::~vtkIzarAxialProjection()
{
}

void vtkIzarAxialProjection::PrintSelf(ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarAxialProjection\n";
}

int vtkIzarAxialProjection::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
{
	vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
	if(!inInfo)
	{
		return 0;
	}
	vtkPointSet* inData = vtkPointSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!inData)
	{
		return 0;
	}
	
	vtkInformation* outInfo = outVector->GetInformationObject(0);
	if(!outInfo)
	{
		return 0;
	}
	vtkPointSet* outData = vtkPointSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!outData)
	{
		return 0;
	}
	outData->ShallowCopy(inData);
	vtkPoints* oldPts = outData->GetPoints();
	if(!oldPts)
	{
		vtkErrorMacro("The dataset does not have any point");
		return 0;
	}
	vtkPoints* newPts = vtkPoints::New();
	newPts->SetDataType(oldPts->GetDataType());
	int N = oldPts->GetNumberOfPoints();
	newPts->SetNumberOfPoints(N);
	
	switch(oldPts->GetDataType())
	{
		vtkTemplateMacro((this->CallSMPAxialProjection<VTK_TT>(oldPts->GetData(), newPts->GetData(), N)));
		default:
			return 0;
			break;
	}
	
	// Perform a deep copy of the original points. We do not want to modify this data array because it is already used by the input data
	vtkDataArray* savedCoords = oldPts->GetData()->NewInstance();
	savedCoords->DeepCopy(oldPts->GetData());
	savedCoords->SetName(SAVED_COORDS);
	if(outData->GetPointData()->HasArray(SAVED_COORDS))
	{
		outData->GetPointData()->RemoveArray(SAVED_COORDS);
	}
	outData->GetPointData()->AddArray(savedCoords);
	savedCoords->Delete();
	outData->SetPoints(newPts);
	newPts->Delete();
	
	return 1;
}
	
	
