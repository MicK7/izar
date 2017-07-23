
#include "vtkIzarRThetaCalculator.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"

vtkStandardNewMacro(vtkIzarRThetaCalculator)

vtkIzarRThetaCalculator::vtkIzarRThetaCalculator():
	vtkPointSetAlgorithm()
{
	IZAR_WARNING
}

vtkIzarRThetaCalculator::~vtkIzarRThetaCalculator()
{
}

void vtkIzarRThetaCalculator::PrintSelf(ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarRThetaCalculator\n";
}

int vtkIzarRThetaCalculator::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
{
	vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
	if(!inInfo)
	{
		return 0;
	}
	vtkPointSet* inData = vtkPointSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!inData)
	{
		return 0;
	}
	
	vtkInformation* outInfo = outVector->GetInformationObject(0);
	if(!outInfo)
	{
		return 0;
	}
	vtkPointSet* outData = vtkPointSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!outData)
	{
		return 0;
	}
	outData->ShallowCopy(inData);
	int nPts = outData->GetNumberOfPoints();
	vtkPoints* pts = outData->GetPoints();
	if(!pts) { return 1; }
	vtkDataArray* coords = pts->GetData();
	vtkDataArray* rAr = coords->NewInstance();
	rAr->SetNumberOfComponents(1);
	rAr->SetNumberOfTuples(nPts);
	rAr->SetName(RADIUS);
	vtkDataArray* thetaAr = coords->NewInstance();
	thetaAr->SetNumberOfComponents(1);
	thetaAr->SetNumberOfTuples(nPts);
	thetaAr->SetName(THETA);
	
	switch(coords->GetDataType())
	{
		vtkTemplateMacro(this->CallSMPComputeRTheta<VTK_TT>(coords, rAr, thetaAr, nPts));
		default:
			return 0;
			break;
	}
	
	if(outData->GetPointData()->HasArray(RADIUS))
	{
		outData->GetPointData()->RemoveArray(RADIUS);
	}
	if(outData->GetPointData()->HasArray(THETA))
	{
		outData->GetPointData()->RemoveArray(THETA);
	}
	
	outData->GetPointData()->AddArray(rAr);
	outData->GetPointData()->AddArray(thetaAr);
	
	return 1;
}
	
	
