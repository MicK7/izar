#ifndef VTKIZARROTATE_H
#define VTKIZARROTATE_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkPointSetAlgorithm.h"

#include <vector>
#include <string>

class vtkInformation;
class vtkInformationVector;

/**
 * Rotates a dataset
 */
class VTK_EXPORT vtkIzarRotate : public vtkPointSetAlgorithm
{
public:
    static vtkIzarRotate* New();
    vtkTypeMacro(vtkIzarRotate, vtkPointSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    
    vtkGetMacro(RotationSpeed, double)
    vtkSetMacro(RotationSpeed, double)
    vtkGetMacro(RotationSpeedUnit, int)
    vtkSetMacro(RotationSpeedUnit, int)
    vtkGetMacro(TimeUnitInSeconds, double)
    vtkSetMacro(TimeUnitInSeconds, double)
    vtkGetMacro(AngleAtInitialTime, double)
    vtkSetMacro(AngleAtInitialTime, double)
    
protected:
	vtkIzarRotate();
	virtual ~vtkIzarRotate();
	
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	virtual int ExecuteInformation(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	virtual int ComputeInputUpdateExtent(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
	double RotationSpeed;
	/**
	 * 1 : RPM, 0 : rad/s
	 */
	int RotationSpeedUnit;
	double RotationSpeedDim;
	
	double TimeUnitInSeconds;
	
	/**
	 * Angle at time = 0, in radians
	 */
	double AngleAtInitialTime;
	
	/**
	 * This boolean is set to false if the input does not give any
	 * time step
	 */
	bool AreTimestepsProvided;
	double CurrentTimestep;
	
private:
    vtkIzarRotate operator=(const vtkIzarRotate&);
    vtkIzarRotate(const vtkIzarRotate&);
};

#endif
