
#include "vtkIzarRestoreCoordinates.h"
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

vtkStandardNewMacro(vtkIzarRestoreCoordinates)

vtkIzarRestoreCoordinates::vtkIzarRestoreCoordinates():
	vtkPointSetAlgorithm()
{
	IZAR_WARNING
}

vtkIzarRestoreCoordinates::~vtkIzarRestoreCoordinates()
{
}

void vtkIzarRestoreCoordinates::PrintSelf(ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarRestoreCoordinates\n";
}

int vtkIzarRestoreCoordinates::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
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
	vtkPoints* pts = outData->GetPoints();
	if(!pts)
	{
		vtkErrorMacro("The dataset does not have any point");
		return 0;
	}
	vtkDataArray* toRestore = outData->GetPointData()->GetArray(SAVED_COORDS);
	if(!toRestore)
	{
		vtkErrorMacro("The dataset does not have any points to restore (" << SAVED_COORDS << ")");
		return 0;
	}
	pts->SetData(toRestore);
	outData->GetPointData()->RemoveArray(SAVED_COORDS);
	
	return 1;
}
	
	
