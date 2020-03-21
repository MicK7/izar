#ifndef VTKIZARRESTORECOORDINATES_H
#define VTKIZARRESTORECOORDINATES_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkPointSetAlgorithm.h"
#include "vtkSMPTools.h"
#include "vtkMath.h"
#include "vtkDataArray.h"

#include <vector>
#include <string>

class vtkInformation;
class vtkInformationVector;

/**
 * Restore coordinates previously stored in the point data
 */
class VTK_EXPORT vtkIzarRestoreCoordinates : public vtkPointSetAlgorithm
{
public:
    static vtkIzarRestoreCoordinates* New();
    vtkTypeMacro(vtkIzarRestoreCoordinates, vtkPointSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    
protected:
	vtkIzarRestoreCoordinates();
	virtual ~vtkIzarRestoreCoordinates();
	
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
private:
    vtkIzarRestoreCoordinates operator=(const vtkIzarRestoreCoordinates&);
    vtkIzarRestoreCoordinates(const vtkIzarRestoreCoordinates&);
};

#endif
