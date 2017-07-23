
#include "vtkIzarRotateSteady.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkPointSet.h"
#include "vtkMath.h"


vtkStandardNewMacro(vtkIzarRotateSteady)

vtkIzarRotateSteady::vtkIzarRotateSteady():
	vtkPointSetAlgorithm(), Mode(0), AngleInDegrees(0.), AngleInRadians(0.), NumberOfPassagesToShift(0)
{
	IZAR_WARNING
}

vtkIzarRotateSteady::~vtkIzarRotateSteady()
{
}

void vtkIzarRotateSteady::PrintSelf(ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarRotateSteady\n";
}

int vtkIzarRotateSteady::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
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
	outData->Initialize();
	
	double actualAngle = 0.0; // Angle in radians
	if(this->Mode == 0)
	{
		// Mode : degrees
		actualAngle = this->AngleInDegrees * 2.*vtkMath::Pi() / 360.;
	}
	else if(this->Mode == 1)
	{
		// Mode : radians
		actualAngle = this->AngleInRadians;
	}
	else if(this->Mode == 2)
	{
		// Mode : shift a given number of passages
		int zSector = 0;
		bool isOk = Izar::GetIntFieldData(inData, ZSECTOR, zSector);
		if(!isOk)
		{
			vtkErrorMacro("Error: could not find the FieldData named " << ZSECTOR << " in the dataset");
			return 0;
		}
		actualAngle = ((double)this->NumberOfPassagesToShift) * 2.*vtkMath::Pi() / ((double)zSector);
	}
	else
	{
		vtkErrorMacro("Error: Bad Mode, should not happen");
		return 0;
	}
	
	vtkPointSet* tmpRes = vtkPointSet::SafeDownCast(Izar::NewRotatedDataObject(inData, actualAngle));
	outData->ShallowCopy(tmpRes);
	tmpRes->Delete();
	return 1;
}
