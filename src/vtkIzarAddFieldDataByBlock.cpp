
#include "vtkIzarAddFieldDataByBlock.h"
#include "vtkIzarAddFieldData.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataSet.h"
#include "vtkDataSet.h"

#include <sstream>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

vtkStandardNewMacro(vtkIzarAddFieldDataByBlock)

vtkIzarAddFieldDataByBlock::vtkIzarAddFieldDataByBlock() : vtkMultiBlockDataSetAlgorithm()
{
	IZAR_WARNING
}

vtkIzarAddFieldDataByBlock::~vtkIzarAddFieldDataByBlock()
{
	this->ClearAddFieldDataObjects();
}

void vtkIzarAddFieldDataByBlock::PrintSelf(ostream& os, vtkIndent indent)
{
	os << "vtkIzarAddFieldDataByBlock: FileName = " << this->FileName;
}

int vtkIzarAddFieldDataByBlock::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
	return 1;
}

int vtkIzarAddFieldDataByBlock::RequestDataObject(
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
 
int vtkIzarAddFieldDataByBlock::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
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
	
	this->ClearAddFieldDataObjects();
	
	// Open and parse the YAML file
	try
	{
		YAML::Node node = YAML::LoadFile(this->FileName);
		this->IterateOverLeaves(inData, outData, node, std::string("/"));
	}
	catch(const YAML::BadFile& e)
	{
		vtkErrorMacro("Could not open JSON/YAML file " << this->FileName);
		outData->Initialize();
		this->ClearAddFieldDataObjects();
		return 0;
	}
	catch(const Izar::TreeTraversalException& e)
	{
		vtkErrorMacro("Error while adding FieldData to individual blocks:\n" << e.what());
		outData->Initialize();
		this->ClearAddFieldDataObjects();
		return 0;
	}
	
	return 1;
}

void vtkIzarAddFieldDataByBlock::IterateOverLeaves(vtkMultiBlockDataSet* in, vtkMultiBlockDataSet* out, YAML::Node& currentYamlLevel, const std::string& curPath)
{
	// Recursive function which performs the operation at the leaves and iterates over the children if we are not at a leaf
	if(!currentYamlLevel.IsMap())
	{
		throw Izar::TreeTraversalException(std::string("At '") + curPath + std::string("': node in YAML file is not a Map"));
	}
	int nChildren = in->GetNumberOfBlocks();
	for(int ii = 0; ii < nChildren; ii++)
	{
		// Retrieve the name of the block
		if(!in->HasMetaData(ii))
		{
			std::stringstream err;
			err << "At '" << curPath << "': VTK node " << ii << " has no MetaData";
			throw Izar::TreeTraversalException(err.str());
		}
		vtkInformation* meta = in->GetMetaData(ii);
		if(!meta->Has(vtkCompositeDataSet::NAME()))
		{
			std::stringstream err;
			err << "At '" << curPath << "': VTK node " << ii << " has no NAME MetaData";
			throw Izar::TreeTraversalException(err.str());
		}
		std::string nodeName = meta->Get(vtkCompositeDataSet::NAME());
		std::stringstream nextPath;
		nextPath << curPath << "/" << nodeName;
		std::string path = nextPath.str();
		
		// Check that the next YAML node exists
		if(!currentYamlLevel[nodeName])
		{
			std::stringstream err;
			err << "At '" << curPath << "': node in YAML file has no child named " << nodeName;
			throw Izar::TreeTraversalException(err.str());
		}
		YAML::Node nextYaml = currentYamlLevel[nodeName];
		vtkDataObject* bkIn = in->GetBlock(ii);
		if(bkIn == NULL)
		{
			// If empty leaf, do nothing (bkOut is already NULL)
		}
		else if(bkIn->IsA("vtkMultiBlockDataSet"))
		{
			// If the child is a multiblock, apply the same function
			vtkMultiBlockDataSet* mbIn = vtkMultiBlockDataSet::SafeDownCast(bkIn);
			vtkMultiBlockDataSet* mbOut = vtkMultiBlockDataSet::SafeDownCast(out->GetBlock(ii));
			this->IterateOverLeaves(mbIn, mbOut, nextYaml, path);
		}
		else if(bkIn->IsA("vtkDataSet"))
		{
			// If the child is a monoblock data set, add the FieldData
			if(!nextYaml.IsMap())
			{
				std::stringstream err;
				err << "At '" << path << "': node in YAML file is not a Map";
				throw Izar::TreeTraversalException(err.str());
			}
			// If the child is not a multiblock, add the FieldData
			vtkIzarAddFieldData* addFld = vtkIzarAddFieldData::New();
			this->AddFieldDataObjects.push_back(addFld);
			addFld->SetInputData(vtkDataSet::SafeDownCast(bkIn));
			for(YAML::const_iterator it = nextYaml.begin(); it != nextYaml.end(); it++)
			{
				std::string fldName = it->first.as<std::string>();
				// Try to read int, else try to read double
				try
				{
					int val = it->second.as<int>();
					addFld->AddIntData(fldName, val);
				}
				catch(const YAML::BadConversion& e)
				{
					try
					{
						double val = it->second.as<double>();
						addFld->AddDoubleData(fldName, val);
					}
					catch(const YAML::BadConversion& e)
					{
						std::stringstream err;
						err << "At '" << path <<"': key in YAML file '" << fldName << "' cannot be cast into an int or a double";
						throw Izar::TreeTraversalException(err.str());
					}
				}
			}
			addFld->Update();
			out->SetBlock(ii, addFld->GetOutput());
		}
		else
		{
			std::stringstream err;
			err << "At '" << path << "': Block is neither a MultiBlockDataSet or a vtkDataSet (it is a " << bkIn->GetClassName() << ")";
			throw Izar::TreeTraversalException(err.str());
		}
	}
}

void vtkIzarAddFieldDataByBlock::ClearAddFieldDataObjects()
{
	for(int ii = 0; ii < this->AddFieldDataObjects.size(); ii++)
	{
		this->AddFieldDataObjects[ii]->Delete();
	}
	this->AddFieldDataObjects.clear();
}
