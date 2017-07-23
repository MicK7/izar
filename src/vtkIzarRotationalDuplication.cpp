
#include "vtkIzarRotationalDuplication.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMath.h"

#include <sstream>

vtkStandardNewMacro(vtkIzarRotationalDuplication)

vtkIzarRotationalDuplication::vtkIzarRotationalDuplication():
	IndexMin(0), IndexMax(0), ForceSectorPeriodicity(0), ZSector(1)
{
	IZAR_WARNING
}

vtkIzarRotationalDuplication::~vtkIzarRotationalDuplication()
{
}

void vtkIzarRotationalDuplication::PrintSelf(ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarRotationalDuplication\n";
}


int vtkIzarRotationalDuplication::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
	return 1;
}

int vtkIzarRotationalDuplication::RequestDataObject(
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

int vtkIzarRotationalDuplication::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
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
	std::cout << "Initializing the out data structure" << std::endl;
	
	// Initialize the structure
	outData->Initialize();
	if(this->IndexMin > this->IndexMax)
	{
		vtkErrorMacro("IndexMin should be smaller or equal than IndexMax");
		return 0;
	}
	int nBlocks = this->IndexMax - this->IndexMin + 1;
	outData->SetNumberOfBlocks(nBlocks);
	std::string baseName("Duplication_");
	for(int ii = 0; ii < nBlocks; ii++)
	{
		std::ostringstream os;
		os << baseName << (this->IndexMin + ii);
		vtkInformation* info = outData->GetMetaData(ii);
		info->Set(vtkCompositeDataSet::NAME(), os.str().c_str());
	}
	
	if(this->ForceSectorPeriodicity)
	{
		// We use the ZSector given by the user
		for(int ii = 0; ii < nBlocks; ii++)
		{
			double angle = ((double)(this->IndexMin + ii)) * 2.*vtkMath::Pi() / ((double)this->ZSector);
			vtkDataObject* newDO = Izar::NewRotatedDataObject(inData, angle);
			outData->SetBlock(ii, newDO);
			newDO->Delete();
		}
	}
	else
	{
		// We use the ZSector given in the FieldData
		RotationTransformFromZSectorGenerator generator;
		for(int ii = 0; ii < nBlocks; ii++)
		{
			generator.rotationIndex = this->IndexMin + ii;
			vtkDataObject* newDO = Izar::NewTransformedDataObject(inData, generator);
			outData->SetBlock(ii, newDO);
			newDO->Delete();
		}
	}
	
	return 1;
}
