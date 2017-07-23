
#include "vtkIzarChorochronicDuplication.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkMath.h"

#include <sstream>
#include <algorithm>

vtkStandardNewMacro(vtkIzarChorochronicDuplication)

vtkIzarChorochronicDuplication::vtkIzarChorochronicDuplication():
	ModeRotationSpeed(0.0),
	RotationSpeedUnit(0),
	NumberOfDiameters(1),
	ZSector(1),
	IndexMin(0),
	IndexMax(1),
	TimeUnitInSeconds(1.0),
	RotationSpeedDim(0.0),
	TimestepsPeriod(1.0),
	DeltaTChoro(1.0),
	DeltaTheta(1.0),
	TemporalPeriod(1.0),
	NTimestepsPerPeriod(1),
	NTimestepsPerChoroPeriod(1),
	TimeTolerance(1.e-9),
	CurrentRequestIndex(-1),
	TimestepsAreOK(true),
	InfoErrorMsg("")
{
	IZAR_WARNING
}

vtkIzarChorochronicDuplication::~vtkIzarChorochronicDuplication()
{
}

void vtkIzarChorochronicDuplication::PrintSelf(ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarChorochronicDuplication\n";
}


int vtkIzarChorochronicDuplication::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
	return 1;
}

int vtkIzarChorochronicDuplication::RequestDataObject(
  vtkInformation* request,
  vtkInformationVector** inputVector ,
  vtkInformationVector* outputVector)
{
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	if(!outInfo)
	{
		return 0;
	}
	vtkMultiBlockDataSet* outData = vtkMultiBlockDataSet::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!outData)
	{
		vtkMultiBlockDataSet* newOutData = vtkMultiBlockDataSet::New();
		outInfo->Set(vtkDataObject::DATA_OBJECT(), newOutData);
		newOutData->Delete();
	}
	return 1;
}

int vtkIzarChorochronicDuplication::RequestInformation(
	vtkInformation* request,
	vtkInformationVector** inputVector ,
	vtkInformationVector* outputVector)
{
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	// Get the input timesteps and save them in the InputTimesteps vector
	if(!inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
	{
		vtkErrorMacro("The input dataset does not contain temporal data");
		return 0;
	}
	int n = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	double* timesteps = inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
	this->InputTimesteps.resize(n);
	for(int ii = 0; ii < n; ii++)
	{
		this->InputTimesteps[ii] = timesteps[ii];
	}
	
	// Determines if this is suitable for a chorochronic reconstruction
	// We have to make sure that no RequestData is made if this fails. However,
	// this can fail (eg at the first time step), and therefore we have to
	// deal with this case without sending error messages
	// To do so, we set an error flag if this fails, and when requestdata is
	// executed, this flag is checked and the corresponding error message
	// is sent if necessary
	this->SetAndCheckTemporalData();
	
	// Set the output information time steps
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	if(!outInfo) { return 0; }
	outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), this->OutputTimesteps.data(), this->OutputTimesteps.size());
	
	this->CurrentRequestIndex = -1; // To indicate that we are not currently looping
	
	return 1;
}

int vtkIzarChorochronicDuplication::RequestUpdateExtent(vtkInformation* request,
	vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	if(this->CurrentRequestIndex == -1) // We are not looping yet
	{
		// Check that the input has provided us valid timesteps
		if(!this->TimestepsAreOK)
		{
			vtkErrorMacro("Error in requestinformation" << this->InfoErrorMsg);
			return 0;
		}
		
		vtkInformation* outInfo = outputVector->GetInformationObject(0);
		if(!outInfo) { return 0; }
		if(!outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
		{
			vtkErrorMacro("The output does not request any time step...");
			return 0;
		}
		double currentTimestep = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
		this->InitializeRequest(currentTimestep);
		this->CurrentRequestIndex = 0;
	}
	if((this->CurrentRequestIndex < 0) || (this->CurrentRequestIndex >= this->TimestepsRequests.size()))
	{
		vtkErrorMacro("Should not happen");
		return 0;
	}
	
	// We are now requesting from the input the timestep at index this->CurrentRequestIndex
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	if(!inInfo) { this->CurrentRequestIndex == -1; return 0; }
	inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->TimestepsRequests[this->CurrentRequestIndex]);
	return 1;
}

int vtkIzarChorochronicDuplication::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
{
	vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
	if(!inInfo)
	{
		return 0;
	}
	vtkDataObject* inData = inInfo->Get(vtkDataObject::DATA_OBJECT());
	if(!inData)
	{
		return 0;
	}
	
	vtkInformation* outInfo = outVector->GetInformationObject(0);
	if(!outInfo)
	{
		return 0;
	}
	vtkMultiBlockDataSet* outData = vtkMultiBlockDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!outData)
	{
		return 0;
	}
	
	// If this is the first RequestData, initializes the result
	if(this->CurrentRequestIndex == 0)
	{
		int nBlocks = this->IndexMax - this->IndexMin + 1;
		outData->Initialize();
		outData->SetNumberOfBlocks(nBlocks);
		for(int ii = 0; ii < nBlocks; ii++)
		{
			vtkInformation* meta = outData->GetMetaData(ii);
			std::ostringstream ss;
			ss << "Duplication_" << this->IndexMin + ii;
			meta->Set(vtkMultiBlockDataSet::NAME(), ss.str().c_str());
		}
	}
	// The input data contains the dataset at timestep this->TimestepsRequest[this->CurrentRequestIndex]
	// Find the duplication(s) which match the requested timestep
	for(auto it = std::begin(this->DuplicationsTimesteps), end = std::end(this->DuplicationsTimesteps); it != end; it++)
	{
		if(this->TimestepsRequests[this->CurrentRequestIndex] == it->second)
		{
			double duplicationAngle = this->DeltaTheta * double(it->first);
			vtkDataObject* curDuplication = Izar::NewRotatedDataObject(inData, duplicationAngle);
			outData->SetBlock(it->first - this->IndexMin, curDuplication);
			curDuplication->Delete();
		}
	}
	// If this was the last RequestData, tell the pipeline not to continue
	if(this->CurrentRequestIndex == this->TimestepsRequests.size() -1)
	{
		request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
		this->CurrentRequestIndex = -1;
	}
	else
	{
		this->CurrentRequestIndex++;
		request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
	}
	return 1;
}

void vtkIzarChorochronicDuplication::SetAndCheckTemporalData()
{
	/*
	 * Called after receiving the input temporal information
	 */
	// First : order the timesteps
	std::sort(std::begin(this->InputTimesteps), std::end(this->InputTimesteps));
	// This is the exepected period between two timesteps, we will check later on if this is really the case
	this->TimestepsPeriod = (this->InputTimesteps.back() - this->InputTimesteps.front()) / double(this->InputTimesteps.size()-1);
	// But first we have to compute the tolerance so that we can perform comparisons between time steps
	this->TimeTolerance = this->TimestepsPeriod * 1.e-3;
	
	// Check that the provided time steps are evenly spaced
	for(int ii = 0; ii < this->InputTimesteps.size() ; ii++)
	{
		if(! this->TimestepsAreEqual(this->InputTimesteps[ii], this->InputTimesteps.front() + (double(ii) * this->TimestepsPeriod)))
		{
			std::ostringstream ss;
			ss << "The provided timesteps are not evenly spaced. Timesteps at position " << ii << " is " << this->InputTimesteps[ii] << " but should be " << (this->InputTimesteps.front() + (double(ii) * this->TimestepsPeriod));
			this->InfoErrorMsg = ss.str();
			this->TimestepsAreOK = false;
			return;
		}
	}
	
	// Temporal periodicity
	if(this->RotationSpeedUnit)
	{
		this->RotationSpeedDim = this->ModeRotationSpeed * 2.*vtkMath::Pi() / 60. * this->TimeUnitInSeconds;
	}
	else
	{
		this->RotationSpeedDim = this->ModeRotationSpeed * this->TimeUnitInSeconds;
	}
	this->TemporalPeriod = 2.*vtkMath::Pi() / (fabs(this->RotationSpeedDim) * double(this->NumberOfDiameters));
	
	this->NTimestepsPerPeriod = int(std::round(this->TemporalPeriod / this->TimestepsPeriod));
	// We check that the number of timesteps per period is an integer
	if(!this->TimestepsAreEqual(this->TemporalPeriod, double(this->NTimestepsPerPeriod)*this->TimestepsPeriod))
	{
		this->InfoErrorMsg = "The time step does not divide the temporal period";
		this->TimestepsAreOK = false;
		return;
	}
	
	// We then check that we have enough time steps to cover at least one period
	// We support providing several periods of the oscillation
	// (it might be useful to check for the convergence of the computation)
	double inputPeriod = this->InputTimesteps.back() - this->InputTimesteps.front() + this->TimestepsPeriod;
	if( ! (this->TimestepsAreEqual(inputPeriod, this->TemporalPeriod) || (inputPeriod > this->TemporalPeriod)))
	{
		std::ostringstream ss;
		ss << "The period provided is too short (" << inputPeriod << " time units while we need at least " << this->TemporalPeriod << " units)";
		this->InfoErrorMsg = ss.str();
		this->TimestepsAreOK = false;
		return;
	}
	
	// We select the output time steps
	this->OutputTimesteps.resize(this->NTimestepsPerPeriod);
	std::copy(std::end(this->InputTimesteps) - this->NTimestepsPerPeriod, std::end(this->InputTimesteps), std::begin(this->OutputTimesteps));
	
	// Chorochronic periodicity
	this->DeltaTheta = 2.*vtkMath::Pi() / double(this->ZSector);
	this->DeltaTChoro = 2.*vtkMath::Pi() / (this->RotationSpeedDim * double(this->ZSector));
	this->NTimestepsPerChoroPeriod = int(std::round(this->DeltaTChoro / this->TimestepsPeriod));
	
	// Check that the number of timesteps per chorochornic period is really an integer
	if(!this->TimestepsAreEqual(this->DeltaTChoro, double(this->NTimestepsPerChoroPeriod)*this->TimestepsPeriod))
	{
		this->InfoErrorMsg = "The time step does not divide the chorochronic period";
		this->TimestepsAreOK = false;
		return;
	}
	
	// This may not be the best place to do that, but we have to make sure that
	// this check is done before entering the request update extent / request data loop
	if(this->IndexMax < this->IndexMin)
	{
		this->InfoErrorMsg = "IndexMax should be larger than IndexMin";
		this->TimestepsAreOK = false;
		return;
	}
	
	this->InfoErrorMsg = "";
	this->TimestepsAreOK = true;
}

void vtkIzarChorochronicDuplication::InitializeRequest(double timestep)
{
	int indexInPeriod = Izar::PositiveModulo(vtkMath::Round((timestep - this->OutputTimesteps[0]) / this->TimestepsPeriod), this->NTimestepsPerPeriod);
	std::set<int> indicesTraversed;
	this->DuplicationsTimesteps.clear();
	this->TimestepsRequests.clear();
	for(int ii = this->IndexMin; ii <= this->IndexMax; ii++)
	{
		int indexRotated = Izar::PositiveModulo(indexInPeriod - ii * this->NTimestepsPerChoroPeriod, this->NTimestepsPerPeriod);
		this->DuplicationsTimesteps.insert(std::pair<int, double>(ii, this->OutputTimesteps[indexRotated]));
		if(indicesTraversed.find(indexRotated) == indicesTraversed.end())
		{
			indicesTraversed.insert(indexRotated);
			this->TimestepsRequests.push_back(this->OutputTimesteps[indexRotated]);
		}
	}
}
