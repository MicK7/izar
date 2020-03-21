
#include "vtkIzarExtractSubsetFromFile.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataSet.h"
#include "vtkStructuredGrid.h"
#include "vtkExtractGrid.h"

#include <sstream>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

vtkStandardNewMacro(vtkIzarExtractSubsetFromFile)

vtkIzarExtractSubsetFromFile::vtkIzarExtractSubsetFromFile() : vtkMultiBlockDataSetAlgorithm()
{
	IZAR_WARNING
}

vtkIzarExtractSubsetFromFile::~vtkIzarExtractSubsetFromFile()
{
}

void vtkIzarExtractSubsetFromFile::PrintSelf(ostream& os, vtkIndent indent)
{
	os << "vtkIzarExtractSubsetFromFile: FileName = " << this->FileName;
}

int vtkIzarExtractSubsetFromFile::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
	return 1;
}

int vtkIzarExtractSubsetFromFile::RequestDataObject(
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

int vtkIzarExtractSubsetFromFile::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
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
	try
	{
		YAML::Node topNode = YAML::LoadFile(this->FileName);
		this->FillTree(inData, outData, topNode, std::string("/"));
	}
	catch(const YAML::BadFile& e)
	{
		vtkErrorMacro("Could not open file " << this->FileName);
		return 0;
	}
	catch(const Izar::TreeTraversalException& e)
	{
		vtkErrorMacro("Error while trying to generate the patches : " << e.what());
		outData->Initialize();
		return 0;
	}
	return 1;
}

void vtkIzarExtractSubsetFromFile::FillTree(vtkMultiBlockDataSet* in, vtkMultiBlockDataSet* curNodeOut, const YAML::Node& curYaml, const std::string& curYamlPath)
{
	/* Assumes that:
	 *   - in contains the full input tree
	 *   - curNodeOut contains the multiBlockDataSet to fill at the current tree level
	 *   - curYaml contains the YAML node at the current tree level
	 */
	if(!curYaml.IsMap())
	{
		std::stringstream err;
		err << "Node " << curYamlPath << " in YAML file is not a Map object";
		throw Izar::TreeTraversalException(err.str());
	}
	for(YAML::const_iterator it = curYaml.begin(); it != curYaml.end(); it++)
	{
		std::string yamlNodeName = it->first.as<std::string>();
		std::string nextYamlPath = curYamlPath + std::string("/") + yamlNodeName;
		YAML::Node nextNode = it->second;
		if(!nextNode.IsMap())
		{
			std::stringstream err;
			err << "Node " << nextYamlPath << " is not a Map object";
			throw Izar::TreeTraversalException(err.str());
		}
		if(this->IsPatchNode(nextNode))
		{
			int imin, imax, jmin, jmax, kmin, kmax;
			std::string blockPath;
			int ijk[6];
			this->GetPatchNodeData(nextNode, blockPath, ijk, nextYamlPath);
			//vtkStructuredGrid* block = this->GetBlockFromPath(in, blockPath, blockPath);
			vtkStructuredGrid* block;
			try
			{
				block = Izar::GetStructuredBlockFromPath(in, blockPath);
			}
			catch(const Izar::BlockNotFound& e)
			{
				std::stringstream err;
				err << "At JSON/YAML node " << curYamlPath << ": trying to find block at path " << blockPath << "\n" << e.what();
				throw Izar::TreeTraversalException(err.str());
			}

			vtkStructuredGrid* patch = Izar::NewStructuredGridFromVOI(block, ijk);
			int nBks = curNodeOut->GetNumberOfBlocks();
			curNodeOut->SetNumberOfBlocks(nBks+1);
			curNodeOut->SetBlock(nBks, patch);
			if(patch)
			{
				patch->Delete();
			}
			vtkInformation* metaData = curNodeOut->GetMetaData(nBks);
			metaData->Set(vtkCompositeDataSet::NAME(), yamlNodeName.c_str());
		}
		else if(this->IsBlockNode(nextNode))
		{
			YAML::Node nodePath = nextNode["path"];
			std::string blockPath;
			try
			{
				blockPath = nodePath.as<std::string>();
			}
			catch(const YAML::BadConversion& e)
			{
				std::stringstream err;
				err << "At JSON/YAML node " << curYamlPath << ": node 'path' is not of type string";
				throw Izar::TreeTraversalException(err.str());
			}
			vtkDataObject* block;
			try
			{
				block = Izar::GetBlockFromPath(in, blockPath);
			}
			catch(const Izar::BlockNotFound& e)
			{
				std::stringstream err;
				err << "At JSON/YAML node " << curYamlPath << ": trying to find block at path " << blockPath << "\n" << e.what();
				throw Izar::TreeTraversalException(err.str());
			}
			
			int nBks = curNodeOut->GetNumberOfBlocks();
			curNodeOut->SetNumberOfBlocks(nBks+1);
			curNodeOut->SetBlock(nBks, block);
			vtkInformation* metaData = curNodeOut->GetMetaData(nBks);
			metaData->Set(vtkCompositeDataSet::NAME(), yamlNodeName.c_str());
		}
		else
		{
			vtkMultiBlockDataSet* nextMultiBlockNode = vtkMultiBlockDataSet::New();
			this->FillTree(in, nextMultiBlockNode, nextNode, nextYamlPath);
			int nBks = curNodeOut->GetNumberOfBlocks();
			curNodeOut->SetNumberOfBlocks(nBks+1);
			curNodeOut->SetBlock(nBks, nextMultiBlockNode);
			vtkInformation* metaData = curNodeOut->GetMetaData(nBks);
			metaData->Set(vtkCompositeDataSet::NAME(), yamlNodeName.c_str());
		}
	}
}

bool vtkIzarExtractSubsetFromFile::IsPatchNode(const YAML::Node& node)
{
	return (node.size() == 2) && (node["path"]) && (node["ijk"]);
}

bool vtkIzarExtractSubsetFromFile::IsBlockNode(const YAML::Node& node)
{
	return (node.size() == 1) && (node["path"]);
}

void vtkIzarExtractSubsetFromFile::GetPatchNodeData(const YAML::Node& node,
	std::string& path, int* ijk, const std::string& yamlPath)
{
	if((!node["path"]) || (!node["ijk"]))
	{
		std::stringstream err;
		err << "Bad node format at node " << yamlPath << ": no 'path' or 'ijk' node found";
		throw Izar::TreeTraversalException(err.str());
	}
	YAML::Node nodePath = node["path"];
	YAML::Node nodeIjk = node["ijk"];
	if(!nodePath.IsScalar())
	{
		std::stringstream err;
		err << "At node " << yamlPath << ": 'path' is not a string object";
		throw Izar::TreeTraversalException(err.str());
	}
	if((!nodeIjk.IsSequence()) || (nodeIjk.size() != 6))
	{
		std::stringstream err;
		err << "At node " << yamlPath << ": 'ijk' is not a list object with 6 elements";
		throw Izar::TreeTraversalException(err.str());
	}
	try
	{
		path = nodePath.as<std::string>();
	}
	catch(const YAML::BadConversion& e)
	{
		std::stringstream err;
		err << "At node " << yamlPath << ": 'path' is not a string object";
		throw Izar::TreeTraversalException(err.str());
	}
	for(int ii = 0; ii < 6; ii++)
	{
		try
		{
			ijk[ii] = nodeIjk[ii].as<int>();
		}
		catch(const YAML::BadConversion& e)
		{
			std::stringstream err;
			err << "At node " << yamlPath << ": ijk element number " << ii << " is not an int";
			throw Izar::TreeTraversalException(err.str());
		}
	}
}

