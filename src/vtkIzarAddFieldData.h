#ifndef VTKIZARADDFIELDDATA_H
#define VTKIZARADDFIELDDATA_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkDataSetAlgorithm.h"

#include <vector>
#include <string>

class vtkInformation;
class vtkInformationVector;

/**
 * Add arbitrary FieldData variables to a vtkDataSet
 * 
 * Does not contain any interface to be used as a ParaView plugin. This
 * interface has to be implemented by derived classes.
 */
class VTK_EXPORT vtkIzarAddFieldData : public vtkDataSetAlgorithm
{
public:
    static vtkIzarAddFieldData* New();
    vtkTypeMacro(vtkIzarAddFieldData, vtkDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    
    /**
     * Add an int to the data to add in the FieldData
     */
	void AddIntData(const std::string& name, int data);
	
	/**
	 * Add a double to the data to add in the FieldData
	 */
	void AddDoubleData(const std::string& name, double data);
        
    /**
     * Clear the data names and values previously added
     */
    void clearData();
    
protected:
	vtkIzarAddFieldData();
	virtual ~vtkIzarAddFieldData();
	
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
	/**
	 * Names of the integers to add in the FieldData
	 */
	std::vector<std::string> intDataNames;
	/**
	 * Values of the integers to add in the FieldData
	 */
	std::vector<int> intDataValues;
	
	/**
	 * Names of the doubles to add in the FieldData
	 */
	std::vector<std::string> doubleDataNames;
	/**
	 * Values of the doubles to add in the FieldData
	 */
	std::vector<double> doubleDataValues;
	
private:
    vtkIzarAddFieldData operator=(const vtkIzarAddFieldData&);
    vtkIzarAddFieldData(const vtkIzarAddFieldData&);
};

#endif
