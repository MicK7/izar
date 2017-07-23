#ifndef VTKIZARROTATIONALDUPLICATION_H
#define VTKIZARROTATIONALDUPLICATION_H

#include "IzarDefines.h"
#include "IzarHelpers.h"

#include "vtkSetGet.h"
#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkTransform.h"

#include <vector>
#include <string>
#include <iostream>

class vtkInformation;
class vtkInformationVector;

/**
 * Duplicates the dataset, using the information stored as ZSECTOR in
 * the FieldData
 */
class VTK_EXPORT vtkIzarRotationalDuplication : public vtkMultiBlockDataSetAlgorithm
{
public:
    static vtkIzarRotationalDuplication* New();
    vtkTypeMacro(vtkIzarRotationalDuplication, vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
	
	vtkGetMacro(IndexMin, int)
	vtkSetMacro(IndexMin, int)
	vtkGetMacro(IndexMax, int)
	vtkSetMacro(IndexMax, int)
	
	vtkGetMacro(ForceSectorPeriodicity, int)
	vtkSetMacro(ForceSectorPeriodicity, int)
	vtkGetMacro(ZSector, int)
	vtkSetMacro(ZSector, int)
    
protected:
	vtkIzarRotationalDuplication();
	virtual ~vtkIzarRotationalDuplication();
	virtual int RequestDataObject(vtkInformation* request,
		vtkInformationVector** inputVector, vtkInformationVector* outputVector);
	virtual int FillInputPortInformation(int port, vtkInformation* info);
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
	/**
	 * The minimal and maximal indices of the passages represented
	 * in the duplicated dataset
	 */
	int IndexMin;
	int IndexMax;
	
	/**
	 * By default, the filter looks for the ZSECTOR FieldData to get
	 * the pitchwise extent of the dataset, and therefore the angle to
	 * rotate when duplicating. Activating this option forces the use
	 * of a user specified pitchwise extent, defined by ZSector.
	 */
	int ForceSectorPeriodicity;
	int ZSector;
	
private:
    vtkIzarRotationalDuplication operator=(const vtkIzarRotationalDuplication&);
    vtkIzarRotationalDuplication(const vtkIzarRotationalDuplication&);
};

/**
 * This functor generates a transform representing a rotation, according
 * to a rotation index given as a parameter and the ZSector FieldData
 */
class RotationTransformFromZSectorGenerator : public Izar::TransformGenerator
{
public:
	RotationTransformFromZSectorGenerator(): rotationIndex(0) { trans = vtkTransform::New(); }
	~RotationTransformFromZSectorGenerator() { trans->Delete(); }
	int rotationIndex;
	inline virtual vtkTransform* operator() (vtkDataSet* data)
	{
		trans->Identity();
		int ZSector;
		if(!Izar::GetIntFieldData(data, ZSECTOR, ZSector))
		{
			std::cout << "The data set does not contain a " << ZSECTOR << " FieldData" << std::endl;
		}
		else
		{
			trans->RotateX(((double)rotationIndex) * 360. / ((double)ZSector));
		}
		return trans;
	}
protected:
	vtkTransform* trans;
};

#endif
