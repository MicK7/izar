#ifndef VTKIZARAXIALPROJECTION_H
#define VTKIZARAXIALPROJECTION_H

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
 * Worker to perform the axial projection
 * Maybe we don't need that...
 */
template<typename CoordsType>
class vtkSMPIzarAxialProjectionOp
{
public:
	CoordsType* oldCoords;
	CoordsType* newCoords;
	void operator() (vtkIdType begin, vtkIdType end)
	{
		for(vtkIdType ii = begin; ii < end; ii++)
		{
			newCoords[ii*3] = 0.;
			newCoords[ii*3+1] = oldCoords[ii*3+1];
			newCoords[ii*3+2] = oldCoords[ii*3+2];
		}
	}
};

/**
 * Projects a dataset
 */
class VTK_EXPORT vtkIzarAxialProjection : public vtkPointSetAlgorithm
{
public:
    static vtkIzarAxialProjection* New();
    vtkTypeMacro(vtkIzarAxialProjection, vtkPointSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    
protected:
	vtkIzarAxialProjection();
	virtual ~vtkIzarAxialProjection();
	
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
	template<typename CoordsType>
	void CallSMPAxialProjection(vtkDataArray* oldCoords, vtkDataArray* newCoords, int nPts)
	{
		vtkSMPIzarAxialProjectionOp<CoordsType> worker;
		worker.oldCoords = reinterpret_cast<CoordsType*>(oldCoords->GetVoidPointer(0));
		worker.newCoords = reinterpret_cast<CoordsType*>(newCoords->GetVoidPointer(0));
		vtkSMPTools::For(0, nPts, worker);
	}
	
private:
    vtkIzarAxialProjection operator=(const vtkIzarAxialProjection&);
    vtkIzarAxialProjection(const vtkIzarAxialProjection&);
};

#endif
