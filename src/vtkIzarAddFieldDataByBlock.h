#ifndef VTKIZARADDFIELDDATABYBLOCK_H
#define VTKIZARADDFIELDDATABYBLOCK_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkMultiBlockDataSetAlgorithm.h"

#include <vector>
#include <string>

#include <yaml-cpp/yaml.h>

class vtkInformation;
class vtkInformationVector;
class vtkIzarAddFieldData;

/**
 * Allows the user to add FieldData with different values for each block.
 * The input data set has to be a multiblockdataset. The user has to
 * generate a YAML-compatible file outside PV. The YAML file must follow
 * the structure of the multiblock data set. The keys of the last level
 * of the YAML file are then used to generate the FieldData.
 */
class VTK_EXPORT vtkIzarAddFieldDataByBlock : public vtkMultiBlockDataSetAlgorithm
{
public:
    static vtkIzarAddFieldDataByBlock* New();
    vtkTypeMacro(vtkIzarAddFieldDataByBlock, vtkMultiBlockDataSetAlgorithm);
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
	vtkIzarAddFieldDataByBlock();
	virtual ~vtkIzarAddFieldDataByBlock();
	
	virtual int RequestDataObject(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
	virtual int FillInputPortInformation(int port, vtkInformation* info);
	
	/**
	 * Performs a Delete() on every element in AddFieldDataObjects and
	 * clear this vector
	 */
	void ClearAddFieldDataObjects();
	
	/**
	 * Used by RequestData to go through the elements of a multiblockdataset
	 */ 
	void IterateOverLeaves(vtkMultiBlockDataSet* in, vtkMultiBlockDataSet* out, YAML::Node& currentYamlLevel, const std::string& curPath);
	
	/**
	 * The name of the JSON/YAML file containing the properties to add
	 */
	std::string FileName;
	
	/**
	 * The different vtkIzarAddFieldDataObjects used to add FieldData
	 * to individual blocks
	 */
	std::vector<vtkIzarAddFieldData*> AddFieldDataObjects;
	
private:
    vtkIzarAddFieldDataByBlock operator=(const vtkIzarAddFieldDataByBlock&);
    vtkIzarAddFieldDataByBlock(const vtkIzarAddFieldDataByBlock&);
};

#endif
