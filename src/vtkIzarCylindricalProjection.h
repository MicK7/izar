#ifndef VTKIZARCYLINDRICALPROJECTION_H
#define VTKIZARCYLINDRICALPROJECTION_H

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

template<typename CoordsType>
class vtkSMPIzarCylindricalProjectionOp
{
public:
	CoordsType* oldCoords;
	CoordsType* newCoords;
	double ThetaRef;
	double Radius;
	void operator() (vtkIdType begin, vtkIdType end)
	{
		for(vtkIdType ii = begin; ii < end; ii++)
		{
			newCoords[ii*3] = oldCoords[ii*3];
			double theta = fmod(atan2(oldCoords[ii*3+2], oldCoords[ii*3+1]) - ThetaRef + vtkMath::Pi(), 2.*vtkMath::Pi());
			if(theta < 0) { theta += 2.*vtkMath::Pi(); }
			theta = theta + ThetaRef - vtkMath::Pi();
			newCoords[ii*3+1] =  - Radius * theta;
			newCoords[ii*3+2] = 0.0;
		}
	}
};

/**
 * Projects a dataset
 */
class VTK_EXPORT vtkIzarCylindricalProjection : public vtkPointSetAlgorithm
{
public:
    static vtkIzarCylindricalProjection* New();
    vtkTypeMacro(vtkIzarCylindricalProjection, vtkPointSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    
    vtkGetMacro(Radius, double)
    vtkSetMacro(Radius, double)
    vtkGetMacro(ThetaRef, double)
    vtkSetMacro(ThetaRef, double)
    vtkGetMacro(FlattenNormals, int)
    vtkSetMacro(FlattenNormals, int)
    
protected:
	vtkIzarCylindricalProjection();
	virtual ~vtkIzarCylindricalProjection();
	
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
	/**
	 * The radius of the cylinder over which we project the data
	 */
	double Radius;
	
	/**
	 * Theta is computed modulo 2 pi. This reference value is used so that
	 * the values are between ThetaRef - pi and ThetaRef + pi
	 */
	double ThetaRef;
	
	/**
	 * Add normals to make the dataset look better
	 */
	int FlattenNormals;
	
	template<typename CoordsType>
	void CallSMPCylindricalProjection(vtkDataArray* oldCoords, vtkDataArray* newCoords, int nPts)
	{
		vtkSMPIzarCylindricalProjectionOp<CoordsType> worker;
		worker.oldCoords = reinterpret_cast<CoordsType*>(oldCoords->GetVoidPointer(0));
		worker.newCoords = reinterpret_cast<CoordsType*>(newCoords->GetVoidPointer(0));
		worker.ThetaRef = this->ThetaRef;
		worker.Radius = this->Radius;
		vtkSMPTools::For(0, nPts, worker);
	}
	
private:
    vtkIzarCylindricalProjection operator=(const vtkIzarCylindricalProjection&);
    vtkIzarCylindricalProjection(const vtkIzarCylindricalProjection&);
};

#endif
