#ifndef VTKIZAREXTRACTCYLINDRICALCOMPONENTS_H
#define VTKIZAREXTRACTCYLINDRICALCOMPONENTS_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkPointSetAlgorithm.h"
#include "vtkSMPTools.h"
#include "vtkMath.h"
#include "vtkDataArray.h"

#include <string>

class vtkInformation;
class vtkInformationVector;
class vtkDataSetAttributes;
class vtkPointSet;
class vtkDataSet;


/**
 * SMP functor
 */
template <typename CoordsType, typename ValueType>
class vtkSMPIzarExtractCylindricalComponentsOp
{
public:
	/* Inputs */
	CoordsType* coords;
	ValueType* inputVector;
	
	/* Outputs */
	ValueType* outputR;
	ValueType* outputTheta;
	
	void operator() (vtkIdType begin, vtkIdType end)
	{
		for(vtkIdType ii = begin; ii != end; ii++)
		{
			ValueType theta = atan2(coords[3*ii+2], coords[3*ii+1]);
			ValueType cost = cos(theta);
			ValueType sint = sin(theta);
			
			outputR[ii] = inputVector[3*ii+1]*cost + inputVector[3*ii+2]*sint;
			outputTheta[ii] = -inputVector[3*ii+1]*sint + inputVector[3*ii+2]*cost;
		}
	}
};


/**
 * Computes the radial and tangential components of a vector
 */
class VTK_EXPORT vtkIzarExtractCylindricalComponents : public vtkPointSetAlgorithm
{
public:
    static vtkIzarExtractCylindricalComponents* New();
    vtkTypeMacro(vtkIzarExtractCylindricalComponents, vtkPointSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
	
	void SetVectorName(int pipo1, int pipo2, int pipo3, int pipo4, const char* name) { this->VectorName = name; }
    std::string GetVectorName() { return this->VectorName; }
	
protected:
    vtkIzarExtractCylindricalComponents();
    ~vtkIzarExtractCylindricalComponents();

    virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);	
	
	template <typename CoordsType, typename ValueType>
	void CallSMPOp(vtkDataArray* coords, vtkDataArray* inputVector,
		vtkDataArray* outputR, vtkDataArray* outputTheta, int nPts)
	{
		vtkSMPIzarExtractCylindricalComponentsOp<CoordsType, ValueType> worker;
		
		worker.coords = reinterpret_cast<CoordsType*>(coords->GetVoidPointer(0));
		worker.inputVector = reinterpret_cast<ValueType*>(inputVector->GetVoidPointer(0));
		worker.outputR = reinterpret_cast<ValueType*>(outputR->GetVoidPointer(0));
		worker.outputTheta = reinterpret_cast<ValueType*>(outputTheta->GetVoidPointer(0));
		vtkSMPTools::For(0, nPts, worker);
	}
	
	std::string VectorName;
	
private:
    vtkIzarExtractCylindricalComponents operator=(const vtkIzarExtractCylindricalComponents&);
    vtkIzarExtractCylindricalComponents(const vtkIzarExtractCylindricalComponents&);
};

#endif
