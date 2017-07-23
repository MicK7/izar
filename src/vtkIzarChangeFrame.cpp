#include "vtkIzarChangeFrame.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointSet.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkSetGet.h"

#include <sstream>

vtkStandardNewMacro(vtkIzarChangeFrame)

vtkIzarChangeFrame::vtkIzarChangeFrame():
	FrameRotationSpeed(0.),
	FrameName("xx")
{
	IZAR_WARNING
}

vtkIzarChangeFrame::~vtkIzarChangeFrame()
{

}

void vtkIzarChangeFrame::PrintSelf(ostream& os, vtkIndent indent)
{
    os << indent << "vtkIzarChangeFrame";
}

int vtkIzarChangeFrame::RequestData(vtkInformation* request,
	vtkInformationVector** inVector, vtkInformationVector* outVector)
{
	vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
	if(!inInfo)
	{
		return 0;
	}
	vtkPointSet* inData = vtkPointSet::SafeDownCast(
		inInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!inData)
	{
		return 0;
	}
	
	vtkInformation* outInfo = outVector->GetInformationObject(0);
	vtkPointSet* outData = vtkPointSet::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!outData)
	{
		return 0;
	}
	
	outData->ShallowCopy(inData);
	
	if(outData->GetNumberOfPoints() == 0) { return 1; } // Do nothing
	
	vtkDataArray* coords = outData->GetPoints()->GetData();
	vtkDataArray* ts = outData->GetPointData()->GetArray(TS);
	vtkDataArray* ps = outData->GetPointData()->GetArray(PS);
	vtkDataArray* v = outData->GetPointData()->GetArray(VELOCITYABS);
	
	if(!coords) { vtkErrorMacro("Cannot find coordinates array"); return 0; }
	if(coords->GetNumberOfComponents() != 3) { vtkErrorMacro("Cannot find vector coordinates"); return 0; }
	if((!ts) || (ts->GetNumberOfComponents() != 1)) { vtkErrorMacro("Cannot find scalar Temperature field '" << TS << "'"); return 0; }
	if((!ps) || (ps->GetNumberOfComponents() != 1)) { vtkErrorMacro("Cannot find scalar Pressure field '" << PS "'"); return 0; }
	if((!v) || (v->GetNumberOfComponents() != 3)) { vtkErrorMacro("Cannot find Velocity vector field '" << VELOCITYABS << "'"); return 0; }
	
	int N = coords->GetNumberOfTuples();
	if((ts->GetNumberOfTuples() != N) || (ps->GetNumberOfTuples() != N) || (v->GetNumberOfTuples() != N))
	{
		vtkErrorMacro("Number of tuples mismatch");
		return 0;
	}
	
	// This filter could be extended to get rid of this assumption
	if((ps->GetDataType() != ts->GetDataType()) ||
		(ps->GetDataType() != v->GetDataType()))
	{
		vtkErrorMacro("The coordinates, Temperature, Pressure, and Velocity must have the same data type (eg. double)");
		return 0;
	}
	
	double gamma;
	double rgas;
	if(!Izar::GetDoubleFieldData(outData, GAMMA, gamma))
	{
		vtkErrorMacro("Cannot find gamma in the field data: assuming gamma = 1.4");
		gamma = 1.4;
	}
	if(!Izar::GetDoubleFieldData(outData, RGAS, rgas))
	{
		vtkErrorMacro("Cannot find rgas in the field data: assuming rgas = 287");
		rgas = 287.0;
	}
	
	vtkDataArray* ttrel = ts->NewInstance();
	ttrel->SetNumberOfComponents(1);
	ttrel->SetNumberOfTuples(N);
	vtkDataArray* ptrel = ts->NewInstance();
	ptrel->SetNumberOfComponents(1);
	ptrel->SetNumberOfTuples(N);
	vtkDataArray* w = ts->NewInstance();
	w->SetNumberOfComponents(3);
	w->SetNumberOfTuples(N);
	vtkDataArray* Mrel = ts->NewInstance();
	Mrel->SetNumberOfComponents(1);
	Mrel->SetNumberOfTuples(N);
	vtkDataArray* beta = ts->NewInstance();
	beta->SetNumberOfComponents(1);
	beta->SetNumberOfTuples(N);
	
	switch(vtkTemplate2PackMacro(coords->GetDataType(), ts->GetDataType()))
	{
		vtkTemplate2Macro((this->CallSMPChangeFrameOp<VTK_T1, VTK_T2>(
			coords, ts, ps, v, ttrel, ptrel, w, Mrel, beta, gamma, rgas, N)));
		default:
			vtkErrorMacro("Should not happen");
			return 0;
			break;
	}
	
	this->NameAndAddArray(outData, ttrel, GENERIC_TT);
	this->NameAndAddArray(outData, ptrel, GENERIC_PT);
	this->NameAndAddArray(outData, w, GENERIC_VELOCITY);
	this->NameAndAddArray(outData, Mrel, GENERIC_MACH);
	this->NameAndAddArray(outData, beta, GENERIC_BETA);
	
	ttrel->Delete();
	ptrel->Delete();
	w->Delete();
	Mrel->Delete();
	beta->Delete();
	
	return 1;
}

void vtkIzarChangeFrame::NameAndAddArray(vtkPointSet* data, vtkDataArray* ar, const std::string& name)
{
	std::stringstream ss;
	ss << name << "_" << this->FrameName;
	std::string newName = ss.str();
	if(data->GetPointData()->HasArray(newName.c_str()))
	{
		data->GetPointData()->RemoveArray(newName.c_str());
	}
	ar->SetName(newName.c_str());
	data->GetPointData()->AddArray(ar);
}
