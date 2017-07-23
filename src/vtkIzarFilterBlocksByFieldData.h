#ifndef VTKIZARFILTERBLOCKSBYFIELDDATA_H
#define VTKIZARFILTERBLOCKSBYFIELDDATA_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkMultiBlockDataSetAlgorithm.h"

#include <vector>
#include <string>

class vtkInformation;
class vtkInformationVector;
class vtkDataSet;


/**
 * Filters a multi block dataset and keeps only the blocks matching a criteria
 */
class VTK_EXPORT vtkIzarFilterBlocksByFieldData : public vtkMultiBlockDataSetAlgorithm
{
public:
    static vtkIzarFilterBlocksByFieldData* New();
    vtkTypeMacro(vtkIzarFilterBlocksByFieldData, vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    
    void SetFieldName(int pipo1, int pipo2, int pipo3, int pipo4, const char* name) { this->FieldName = name; this->Modified(); }
    const char* GetFieldName() { return this->FieldName.c_str(); }
    
    vtkGetMacro(UseRange, int)
    vtkSetMacro(UseRange, int)
    vtkGetMacro(Value, double)
    vtkSetMacro(Value, double)
    vtkGetMacro(Tolerance, double)
    vtkSetMacro(Tolerance, double)
    vtkGetMacro(KeepEmptyNodes, int)
    vtkSetMacro(KeepEmptyNodes, int)
    
    void SetRange(double valMin, double valMax) { this->RangeMin = valMin; this->RangeMax = valMax; this->Modified(); }
    
protected:
	vtkIzarFilterBlocksByFieldData();
	virtual ~vtkIzarFilterBlocksByFieldData();
	
	virtual int RequestDataObject(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
	virtual int FillInputPortInformation(int port, vtkInformation* info);
	
	/**
	 * Used to test whether a block matches the criteria
	 */
	bool MatchCriterion(vtkDataSet* block);
	
	/**
	 * Recursively removes NULL blocks
	 * Returns true if all the children of data have been removed
	 */
	bool RemoveEmptyNodes(vtkMultiBlockDataSet* data);
	
	/**
	 * Name of the field data array
	 */
	std::string FieldName;
	/**
	 * If true, uses a range instead of a value and a tolerance
	 */
	int UseRange;
	/**
	 * The value, taken into account only if UseRange is false
	 */
	double Value;
	/**
	 * The tolerance, taken into account only if UseRange is false
	 */
	double Tolerance;
	/**
	 * The minimum value of the range, taken into account only if UseRange is true
	 */
	double RangeMin;
	/**
	 * The maximum value of the range, taken into account only if UseRange is true
	 */
	double RangeMax;
	
	/**
	 * Advanced option: by default, the filter removes completely the
	 * blocks that don't match the criteria from the dataset. When this
	 * option is set to true, the blocks are still removed, but the
	 * nodes corresponding to this block now hold a NULL pointer
	 */
	int KeepEmptyNodes;
	
private:
    vtkIzarFilterBlocksByFieldData operator=(const vtkIzarFilterBlocksByFieldData&);
    vtkIzarFilterBlocksByFieldData(const vtkIzarFilterBlocksByFieldData&);
};

#endif
