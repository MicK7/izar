#ifndef VTKIZAREXTRACTSUBSETFROMFILE_H
#define VTKIZAREXTRACTSUBSETFROMFILE_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkExtractGrid.h"

#include <string>

#include <yaml-cpp/yaml.h>

class vtkInformation;
class vtkInformationVector;
class vtkStructuredGrid;

/**
 * Extracts subsets from a multiblockdataset of structured grids, or full blocks,
 * and group them according to a hierarchy defined in a JSON or YAML file.
 * Typically used to extract surface patches.
 * 
 * Format of the JSON/YAML file:
 * {
 * "level1": {
 *   "level2": {
 *    ...
 *      "patch1": {"path": "/path/to/the/block",
 * 				   "ijk": [imin, imax, jmin, jmax, kmin, kmax]
 *		},
 * 		"patch2": {...
 *      },
 *      ...
 *    }
 *  }
 */
class VTK_EXPORT vtkIzarExtractSubsetFromFile : public vtkMultiBlockDataSetAlgorithm
{
public:
    static vtkIzarExtractSubsetFromFile* New();
    vtkTypeMacro(vtkIzarExtractSubsetFromFile, vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    
    /**
     * Standard Setter
     */
    void SetFileName(const char* name) { this->FileName = name; this->Modified(); }
    
    /**
     * Standard Getter
     */
    const char* GetFileName() { return this->FileName.c_str(); }
protected:
	vtkIzarExtractSubsetFromFile();
	virtual ~vtkIzarExtractSubsetFromFile();
	
	virtual int RequestDataObject(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
	virtual int FillInputPortInformation(int port, vtkInformation* info);
	
	/**
	 * Recursive traversal of the JSON/YAML file
	 * Called by RequestData
	 */
	void FillTree(vtkMultiBlockDataSet* in, vtkMultiBlockDataSet* curNodeOut, const YAML::Node& curYaml, const std::string& curYamlPath);
	
	/**
	 * Test whether a JSON/YAML node describes a patch
	 */
	bool IsPatchNode(const YAML::Node& node);
	
	/**
	 * Test whether a JSON/YAML node describes a full block
	 */
	bool IsBlockNode(const YAML::Node& node);
	
	/**
	 * From a node describing a patch, extracts the block path and the VOI
	 */
	void GetPatchNodeData(const YAML::Node& node, std::string& path, int* ijk, const std::string& yamlPath);
	
	/**
	 * Name of the JSON/YAML file
	 */
	std::string FileName;
	
private:
    vtkIzarExtractSubsetFromFile operator=(const vtkIzarExtractSubsetFromFile&);
    vtkIzarExtractSubsetFromFile(const vtkIzarExtractSubsetFromFile&);
};
#endif
