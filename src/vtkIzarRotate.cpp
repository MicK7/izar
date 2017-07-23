
#include "vtkIzarRotate.h"
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
#include "vtkStreamingDemandDrivenPipeline.h"


vtkStandardNewMacro(vtkIzarRotate)

vtkIzarRotate::vtkIzarRotate():
	RotationSpeed(0.),
	RotationSpeedUnit(0),
	RotationSpeedDim(0.),
	TimeUnitInSeconds(1.),
	AngleAtInitialTime(0.),
	AreTimestepsProvided(true),
	CurrentTimestep(0.)
{
	IZAR_WARNING
}

vtkIzarRotate::~vtkIzarRotate()
{
}

void vtkIzarRotate::PrintSelf(ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarRotate\n";
}

int vtkIzarRotate::ExecuteInformation(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
{
    vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outVector->GetInformationObject(0);

    if(inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
    {
        int nTimesteps = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        double* timesteps = new double[nTimesteps];
        inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), timesteps);
        outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), timesteps, nTimesteps);
        delete[] timesteps;
        this->AreTimestepsProvided = true;
    }
    else
    {
        this->AreTimestepsProvided = false;
    }
    // Compute internal data for the rotation
    if(this->RotationSpeedUnit) // RPM
    {
		this->RotationSpeedDim = this->RotationSpeed * 2.*vtkMath::Pi() / 60. * this->TimeUnitInSeconds;
	}
	else // rad/s
	{
		this->RotationSpeedDim = this->RotationSpeed * this->TimeUnitInSeconds;
	}
    return 1;
}

int vtkIzarRotate::ComputeInputUpdateExtent(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
{
    vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outVector->GetInformationObject(0);

    // We simply forward the requested timestep to the input
    if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
    {
        this->CurrentTimestep = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
        inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->CurrentTimestep);
    }
    else
    {
        this->CurrentTimestep = 0;
    }
    return 1;
}

int vtkIzarRotate::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
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
	double angle = this->AngleAtInitialTime + this->RotationSpeedDim * this->CurrentTimestep;
	vtkDataObject* newData = Izar::NewRotatedDataObject(inData, angle);
	vtkPointSet* newPtSet = vtkPointSet::SafeDownCast(newData);
	if(!newPtSet)
	{
		vtkErrorMacro("Should not happen");
		newData->Delete();
		return 0;
	}
	else
	{
		outData->ShallowCopy(newPtSet);
		newData->Delete();
		return 1;
	}
}
