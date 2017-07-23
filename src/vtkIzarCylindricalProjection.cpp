
#include "vtkIzarCylindricalProjection.h"
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

vtkStandardNewMacro(vtkIzarCylindricalProjection)

vtkIzarCylindricalProjection::vtkIzarCylindricalProjection():
	vtkPointSetAlgorithm(), ThetaRef(0.), Radius(1.0), FlattenNormals(0)
{
	IZAR_WARNING
}

vtkIzarCylindricalProjection::~vtkIzarCylindricalProjection()
{
}

void vtkIzarCylindricalProjection::PrintSelf(ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarCylindricalProjection\n";
}

int vtkIzarCylindricalProjection::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
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
		vtkTemplateMacro((this->CallSMPCylindricalProjection<VTK_TT>(oldPts->GetData(), newPts->GetData(), N)));
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
	
	// If the user request it, flatten the normals
	if(this->FlattenNormals)
	{
		vtkDoubleArray* normalsPts = vtkDoubleArray::New();
		
		normalsPts->SetNumberOfComponents(3);
		normalsPts->SetNumberOfTuples(N);
		normalsPts->FillComponent(0, 0.);
		normalsPts->FillComponent(1, 0.);
		normalsPts->FillComponent(2, 1.);
		normalsPts->SetName("Normals");
		outData->GetPointData()->SetNormals(normalsPts);
		normalsPts->Delete();
		
		vtkDoubleArray* normalsCells = vtkDoubleArray::New();
		normalsCells->SetNumberOfComponents(3);
		normalsCells->SetNumberOfTuples(outData->GetNumberOfCells());
		normalsCells->FillComponent(0, 0.);
		normalsCells->FillComponent(1, 0.);
		normalsCells->FillComponent(2, 1.);
		normalsCells->SetName("Normals");
		outData->GetCellData()->SetNormals(normalsCells);
		normalsCells->Delete();
	}
	
	return 1;
}
	
	
