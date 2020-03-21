#ifndef VTKIZARROTATESTEADY_H
#define VTKIZARROTATESTEADY_H

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
class VTK_EXPORT vtkIzarRotateSteady : public vtkPointSetAlgorithm
{
public:
    static vtkIzarRotateSteady* New();
    vtkTypeMacro(vtkIzarRotateSteady, vtkPointSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    vtkGetMacro(Mode, int)
    vtkSetMacro(Mode, int)
    vtkGetMacro(AngleInDegrees, double)
    vtkSetMacro(AngleInDegrees, double)
    vtkGetMacro(AngleInRadians, double)
    vtkSetMacro(AngleInRadians, double)
    vtkGetMacro(NumberOfPassagesToShift, int)
    vtkSetMacro(NumberOfPassagesToShift, int)
    
protected:
	vtkIzarRotateSteady();
	virtual ~vtkIzarRotateSteady();
	
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
	/**
	 * Selects the mode used to compute the angle :
	 *   0 : Angle in degrees
	 *   1 : Angle in radians
	 *   2 : Number of passages to shift, using the ZSector FieldData
	 */
	int Mode;
	
	double AngleInDegrees;
	double AngleInRadians;
	int NumberOfPassagesToShift;
	
private:
    vtkIzarRotateSteady operator=(const vtkIzarRotateSteady&);
    vtkIzarRotateSteady(const vtkIzarRotateSteady&);
};

#endif
