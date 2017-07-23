
#include "vtkIzarFilterBlocksByFieldData.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataSet.h"
#include "vtkDataSet.h"
#include "vtkDataObjectTreeIterator.h"
#include "vtkFieldData.h"
#include "vtkDataArray.h"
#include "vtkIntArray.h"

#include <sstream>
#include <stdexcept>
#include <cmath>

vtkStandardNewMacro(vtkIzarFilterBlocksByFieldData)

vtkIzarFilterBlocksByFieldData::vtkIzarFilterBlocksByFieldData() : vtkMultiBlockDataSetAlgorithm(),
	UseRange(0), Value(0.0), Tolerance(0.0), RangeMin(0.0), RangeMax(0.0)
{
	IZAR_WARNING
}

vtkIzarFilterBlocksByFieldData::~vtkIzarFilterBlocksByFieldData()
{
}

void vtkIzarFilterBlocksByFieldData::PrintSelf(ostream& os, vtkIndent indent)
{
	os << "vtkIzarFilterBlocksByFieldData";
}

int vtkIzarFilterBlocksByFieldData::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
	return 1;
}

int vtkIzarFilterBlocksByFieldData::RequestDataObject(
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
 
int vtkIzarFilterBlocksByFieldData::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
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
	outData->CopyStructure(inData);
	
	// Algorithm: copy the structure, and shallow copy the blocks which respect the criterion
	vtkDataObjectTreeIterator* iterIn = inData->NewTreeIterator();
	vtkDataObjectTreeIterator* iterOut = outData->NewTreeIterator();
	iterIn->SkipEmptyNodesOff();
	iterIn->TraverseSubTreeOn();
	iterIn->VisitOnlyLeavesOn();
	iterOut->SkipEmptyNodesOff();
	iterOut->TraverseSubTreeOn();
	iterOut->VisitOnlyLeavesOn();
	iterIn->InitTraversal();
	iterOut->InitTraversal();
	
	while(!iterIn->IsDoneWithTraversal())
	{
		vtkDataObject* inObj = inData->GetDataSet(iterIn);
		if(inObj)
		{
			vtkDataSet* inDs = vtkDataSet::SafeDownCast(inObj);
			if(this->MatchCriterion(inDs))
			{
				vtkDataSet* outDs = inDs->NewInstance();
				outDs->ShallowCopy(inDs);
				outData->SetDataSet(iterOut, outDs);
				outDs->Delete();
			}
		}
		iterIn->GoToNextItem();
		iterOut->GoToNextItem();
	}
	iterIn->Delete();
	iterOut->Delete();
	
	// Remove the empty nodes
	if(!this->KeepEmptyNodes)
	{
		this->RemoveEmptyNodes(outData);
	}
	
	return 1;
}

bool vtkIzarFilterBlocksByFieldData::RemoveEmptyNodes(vtkMultiBlockDataSet* data)
{
	bool allRemoved = true;
	int nBks = data->GetNumberOfBlocks();
	for(int ii = nBks-1; ii >= 0; ii--)
	{
		vtkDataObject* child = data->GetBlock(ii);
		if(child == NULL)
		{
			data->RemoveBlock(ii);
		}
		else
		{
			vtkMultiBlockDataSet* mbChild = vtkMultiBlockDataSet::SafeDownCast(child);
			if(mbChild)
			{
				if(this->RemoveEmptyNodes(mbChild))
				{
					data->RemoveBlock(ii);
				}
				else
				{
					allRemoved = false;
				}
			}
			else
			{
				allRemoved = false;
			}
		}
	}
	return allRemoved;
}

bool vtkIzarFilterBlocksByFieldData::MatchCriterion(vtkDataSet* block)
{
	vtkFieldData* fData = block->GetFieldData();
	vtkDataArray* ar = fData->GetArray(this->FieldName.c_str());
	if(!ar) { return false; }
	if(ar->GetNumberOfTuples() == 0) { return false; }
	if(this->UseRange)
	{
		double valueField = ar->GetTuple(0)[0];
		return (valueField >= this->RangeMin) && (valueField <= this->RangeMax);
	}
	else
	{
		// If the fielddata is an int, we don't use the tolerance and we
		// round the value
		vtkIntArray* intAr = vtkIntArray::SafeDownCast(ar);
		if(intAr != NULL)
		{
			int valueRounded = static_cast<int>(std::lround(this->Value));
			int valueField = (reinterpret_cast<int*>(intAr->GetVoidPointer(0)))[0];
			return valueRounded == valueField;
		}
		else
		{
			double valueField = ar->GetTuple(0)[0];
			return std::fabs(this->Value - valueField) <= this->Tolerance;
		}
	}
}
