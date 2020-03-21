#ifndef VTKIZARREMOVEEMPTYBLOCKS_H
#define VTKIZARREMOVEEMPTYBLOCKS_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkMultiBlockDataSetAlgorithm.h"

#include <vector>
#include <string>

class vtkInformation;
class vtkInformationVector;
class vtkIzarAddFieldData;

/**
 * Takes as an input a multiblock, returns the same multiblock but without
 * the nodes with 0 points or NULL data
 */
class VTK_EXPORT vtkIzarRemoveEmptyBlocks : public vtkMultiBlockDataSetAlgorithm
{
public:
    static vtkIzarRemoveEmptyBlocks* New();
    vtkTypeMacro(vtkIzarRemoveEmptyBlocks, vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    
protected:
	vtkIzarRemoveEmptyBlocks();
	virtual ~vtkIzarRemoveEmptyBlocks();
	
	virtual int RequestDataObject(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
	virtual int FillInputPortInformation(int port, vtkInformation* info);
	
	/**
	 * Recursive function to do the job
	 * @param orig original dataset
	 * @param toFill empty dataset, destination of the copy
	 */
	static void FillMultiBlockWithNonEmptyNodes(vtkMultiBlockDataSet* orig, vtkMultiBlockDataSet* toFill);
	
private:
    vtkIzarRemoveEmptyBlocks operator=(const vtkIzarRemoveEmptyBlocks&);
    vtkIzarRemoveEmptyBlocks(const vtkIzarRemoveEmptyBlocks&);
};

#endif
