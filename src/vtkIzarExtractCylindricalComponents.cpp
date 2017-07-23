#include "vtkIzarExtractCylindricalComponents.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointSet.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkSetGet.h"

#include <sstream>

vtkStandardNewMacro(vtkIzarExtractCylindricalComponents)

vtkIzarExtractCylindricalComponents::vtkIzarExtractCylindricalComponents():
	VectorName("")
{
	IZAR_WARNING
}

vtkIzarExtractCylindricalComponents::~vtkIzarExtractCylindricalComponents()
{

}

void vtkIzarExtractCylindricalComponents::PrintSelf(ostream& os, vtkIndent indent)
{
    os << indent << "vtkIzarExtractCylindricalComponents";
}

int vtkIzarExtractCylindricalComponents::RequestData(vtkInformation* request,
	vtkInformationVector** inVector, vtkInformationVector* outVector)
{
	vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
	if(!inInfo)
	{
		return 0;
	}
	vtkPointSet* inData = vtkPointSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!inData)
	{
		return 0;
	}
	
	vtkInformation* outInfo = outVector->GetInformationObject(0);
	vtkPointSet* outData = vtkPointSet::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!outData)
	{
		return 0;
	}
	
	outData->ShallowCopy(inData);
	
	if(outData->GetNumberOfPoints() == 0) { return 1; } // Do nothing
	
	vtkDataArray* coords = outData->GetPoints()->GetData();
	vtkDataArray* inputVector = outData->GetPointData()->GetArray(this->VectorName.c_str());
	
	if(!coords) { vtkErrorMacro("Cannot find coordinates array"); return 0; }
	if(coords->GetNumberOfComponents() != 3) { vtkErrorMacro("Cannot find vector coordinates"); return 0; }
	if((!inputVector) || (inputVector->GetNumberOfComponents() != 3)) { vtkErrorMacro("Cannot find vector field '" << this->VectorName << "'"); return 0; }
	
	int N = coords->GetNumberOfTuples();
	if(inputVector->GetNumberOfTuples() != N)
	{
		vtkErrorMacro("Number of tuples mismatch");
		return 0;
	}
	
	vtkDataArray* outputR = inputVector->NewInstance();
	outputR->SetNumberOfComponents(1);
	outputR->SetNumberOfTuples(N);
	vtkDataArray* outputTheta = inputVector->NewInstance();
	outputTheta->SetNumberOfComponents(1);
	outputTheta->SetNumberOfTuples(N);
	
	switch(vtkTemplate2PackMacro(coords->GetDataType(), inputVector->GetDataType()))
	{
		vtkTemplate2Macro((this->CallSMPOp<VTK_T1, VTK_T2>(
			coords, inputVector, outputR, outputTheta, N)));
		default:
			vtkErrorMacro("Should not happen");
			return 0;
			break;
	}
	
	std::string nameR = this->VectorName + "_r";
	outputR->SetName(nameR.c_str());
	std::string nameT = this->VectorName + "_theta";
	outputTheta->SetName(nameT.c_str());
	
	if(outData->GetPointData()->HasArray(nameR.c_str()))
	{
		outData->GetPointData()->RemoveArray(nameR.c_str());
	}
	if(outData->GetPointData()->HasArray(nameT.c_str()))
	{
		outData->GetPointData()->RemoveArray(nameT.c_str());
	}
	outData->GetPointData()->AddArray(outputR);
	outData->GetPointData()->AddArray(outputTheta);
	
	outputR->Delete();
	outputTheta->Delete();
	
	return 1;
}
