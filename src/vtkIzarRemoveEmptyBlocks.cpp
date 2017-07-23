
#include "vtkIzarRemoveEmptyBlocks.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataSet.h"
#include "vtkDataSet.h"

#include <sstream>
#include <stdexcept>

vtkStandardNewMacro(vtkIzarRemoveEmptyBlocks)

vtkIzarRemoveEmptyBlocks::vtkIzarRemoveEmptyBlocks() : vtkMultiBlockDataSetAlgorithm()
{
	IZAR_WARNING
}

vtkIzarRemoveEmptyBlocks::~vtkIzarRemoveEmptyBlocks()
{
}

void vtkIzarRemoveEmptyBlocks::PrintSelf(ostream& os, vtkIndent indent)
{
	os << "vtkIzarRemoveEmptyBlocks";
}

int vtkIzarRemoveEmptyBlocks::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
	return 1;
}

int vtkIzarRemoveEmptyBlocks::RequestDataObject(
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
 
int vtkIzarRemoveEmptyBlocks::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
{
	vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
	if(!inInfo)
	{
		return 0;
	}
	vtkMultiBlockDataSet* inData = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
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
	
	outData->Initialize();
	vtkIzarRemoveEmptyBlocks::FillMultiBlockWithNonEmptyNodes(inData, outData);
	
	return 1;
}

void vtkIzarRemoveEmptyBlocks::FillMultiBlockWithNonEmptyNodes(vtkMultiBlockDataSet* orig, vtkMultiBlockDataSet* toFill)
{
	// toFill is assumed to be empty
	int n = orig->GetNumberOfBlocks();
	int idCurrentBlock = 0;
	for(int ii = 0; ii < n; ii++)
	{
		vtkDataObject* block = orig->GetBlock(ii);
		if(block != NULL)
		{
			// If it is a multiblock
			vtkMultiBlockDataSet* mbChildOrig = vtkMultiBlockDataSet::SafeDownCast(block);
			if(mbChildOrig)
			{
				vtkMultiBlockDataSet* mbChildFilled = vtkMultiBlockDataSet::New();
				vtkIzarRemoveEmptyBlocks::FillMultiBlockWithNonEmptyNodes(mbChildOrig, mbChildFilled);
				if(mbChildFilled->GetNumberOfBlocks() != 0)
				{
					toFill->SetBlock(idCurrentBlock, mbChildFilled);
					if(orig->HasMetaData(ii))
					{
						vtkInformation* newMetaData = toFill->GetMetaData(idCurrentBlock);
						vtkInformation* oldMetaData = orig->GetMetaData(ii);
						newMetaData->Copy(oldMetaData);
					}
					idCurrentBlock++;
				}
				
				mbChildFilled->Delete();
			}
			else
			{
				// If it is a vtkDataSet
				vtkDataSet* dsChildOrig = vtkDataSet::SafeDownCast(block);
				if(dsChildOrig && (dsChildOrig->GetNumberOfPoints() != 0))
				{
					toFill->SetBlock(idCurrentBlock, dsChildOrig);
					if(orig->HasMetaData(ii))
					{
						vtkInformation* newMetaData = toFill->GetMetaData(idCurrentBlock);
						vtkInformation* oldMetaData = orig->GetMetaData(ii);
						newMetaData->Copy(oldMetaData);
					}
					idCurrentBlock++;
				}
			}
		}		
	}
}

