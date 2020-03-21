#ifndef VTKIZARRTHETACALCULATOR_H
#define VTKIZARRTHETACALCULATOR_H

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

template<typename ValueType>
class SMPComputeRThetaFunctor
{
public:
	ValueType* coords;
	ValueType* r;
	ValueType* theta;
	
	void operator() (vtkIdType begin, vtkIdType end)
	{
		for(vtkIdType id = begin ; id < end ; id++)
		{
			r[id] = sqrt(coords[3*id+1]*coords[3*id+1] + coords[3*id+2]*coords[3*id+2]);
			theta[id] = atan2(coords[3*id+2], coords[3*id+1]);
		}
	}
};
	

/**
 * Cylindrical coordinates calculator
 */
class VTK_EXPORT vtkIzarRThetaCalculator : public vtkPointSetAlgorithm
{
public:
    static vtkIzarRThetaCalculator* New();
    vtkTypeMacro(vtkIzarRThetaCalculator, vtkPointSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    
protected:
	vtkIzarRThetaCalculator();
	virtual ~vtkIzarRThetaCalculator();
	
	virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);
	
	template<typename ValueType>
	void CallSMPComputeRTheta(vtkDataArray* coords, vtkDataArray* r,
		vtkDataArray* theta, int nPts)
	{
		SMPComputeRThetaFunctor<ValueType> functor;
		functor.coords = reinterpret_cast<ValueType*>(coords->GetVoidPointer(0));
		functor.r = reinterpret_cast<ValueType*>(r->GetVoidPointer(0));
		functor.theta = reinterpret_cast<ValueType*>(theta->GetVoidPointer(0));
		vtkSMPTools::For(0, nPts, functor);
	}
	
private:
    vtkIzarRThetaCalculator operator=(const vtkIzarRThetaCalculator&);
    vtkIzarRThetaCalculator(const vtkIzarRThetaCalculator&);
};

#endif
