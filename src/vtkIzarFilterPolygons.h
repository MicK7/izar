#ifndef vtkIzarFilterPolygons_H
#define vtkIzarFilterPolygons_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkSMPTools.h"
#include "vtkMath.h"
#include "vtkDataArray.h"

#include <vector>
#include <string>

class vtkInformation;
class vtkInformationVector;

/**
 * A polydata may contain vertices, lines, polygons, and polygon strips
 * This filter takes an input polydata and removes everything except
 * the polygons
 */
class VTK_EXPORT vtkIzarFilterPolygons : public vtkPolyDataAlgorithm
{
public:
    static vtkIzarFilterPolygons* New();
    vtkTypeMacro(vtkIzarFilterPolygons, vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    
protected:
	vtkIzarFilterPolygons();
	virtual ~vtkIzarFilterPolygons();
	
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	virtual int FillInputPortInformation(int port, vtkInformation* info);
	
private:
    vtkIzarFilterPolygons operator=(const vtkIzarFilterPolygons&);
    vtkIzarFilterPolygons(const vtkIzarFilterPolygons&);
};

#endif
