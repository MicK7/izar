#include "vtkIzarReconstructFields.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointSet.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkSetGet.h"

vtkStandardNewMacro(vtkIzarReconstructFields)

vtkIzarReconstructFields::vtkIzarReconstructFields():
	ReplaceVars(1), AdvancedMode(0), SuthConst(110.4),
	SuthMuRef(1.717e-5), SuthTRef(273.0),
	ComputePs(1), ComputeTs(1), ComputePtabs(1), ComputeTtabs(1),
	ComputePtrel(1), ComputeTtrel(1), ComputeS(1), ComputeV(1),
	ComputeW(1), ComputeM(1), ComputeMrel(1), ComputePhi(1),
	ComputeAlpha(1), ComputeBeta(1), ComputeMu(1)
{
	IZAR_WARNING
}

vtkIzarReconstructFields::~vtkIzarReconstructFields()
{

}

void vtkIzarReconstructFields::PrintSelf(ostream& os, vtkIndent indent)
{
    os << indent << "vtkIzarReconstructFields";
}

int vtkIzarReconstructFields::RequestData(vtkInformation* request,
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
	vtkDataArray* rho = outData->GetPointData()->GetArray(RHO);
	vtkDataArray* rhov = outData->GetPointData()->GetArray(RHOV);
	vtkDataArray* rhoe = outData->GetPointData()->GetArray(RHOE);
	
	if(!coords) { vtkErrorMacro("Cannot find coordinates array"); return 0; }
	if(coords->GetNumberOfComponents() != 3) { vtkErrorMacro("Cannot find vector coordinates"); return 0; }
	if((!rho) || (rho->GetNumberOfComponents() != 1)) { vtkErrorMacro("Cannot find scalar Density field '" << RHO << "'"); return 0; }
	if((!rhov) || (rhov->GetNumberOfComponents() != 3)) { vtkErrorMacro("Cannot find vector Momentum field '" << RHOV "'"); return 0; }
	if((!rhoe) || (rhoe->GetNumberOfComponents() != 1)) { vtkErrorMacro("Cannot find scalar EnergyStagnationDensity field '" << RHOE << "'"); return 0; }
	
	int N = coords->GetNumberOfTuples();
	if((rho->GetNumberOfTuples() != N) || (rhov->GetNumberOfTuples() != N) || (rhoe->GetNumberOfTuples() != N))
	{
		vtkErrorMacro("Number of tuples mismatch");
		return 0;
	}
	
	// This filter could be extended to get rid of this assumption
	if((rho->GetDataType() != rhov->GetDataType()) ||
		(rho->GetDataType() != rhoe->GetDataType()))
	{
		vtkErrorMacro("The Density, Momentum, and EnergyStagnationDensity must have the same data type (eg. double)");
		return 0;
	}
	
	// Get the field data
	double omega;
	double gamma;
	double rgas;
	
	if(!Izar::GetDoubleFieldData(outData, OMEGA, omega))
	{
		vtkErrorMacro("Cannot find omega in the field data: assuming omega = 0");
		omega = 0.;
	}
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
	
	// Generate the new arrays
	vtkDataArray* ps = NULL;
	vtkDataArray* ts = NULL;
	vtkDataArray* ptabs = NULL;
	vtkDataArray* ttabs = NULL;
	vtkDataArray* ptrel = NULL;
	vtkDataArray* ttrel = NULL;
	vtkDataArray* s = NULL;
	vtkDataArray* v = NULL;
	vtkDataArray* w = NULL;
	vtkDataArray* M = NULL;
	vtkDataArray* Mrel = NULL;
	vtkDataArray* phi = NULL;
	vtkDataArray* alpha = NULL;
	vtkDataArray* beta = NULL;
	vtkDataArray* mu = NULL;
	
	if(ComputePs) { ps = rho->NewInstance(); ps->SetNumberOfComponents(1); ps->SetNumberOfTuples(N); }
	if(ComputeTs) { ts = rho->NewInstance(); ts->SetNumberOfComponents(1); ts->SetNumberOfTuples(N); }
	if(ComputePtabs) { ptabs = rho->NewInstance(); ptabs->SetNumberOfComponents(1); ptabs->SetNumberOfTuples(N); }
	if(ComputeTtabs) { ttabs = rho->NewInstance(); ttabs->SetNumberOfComponents(1); ttabs->SetNumberOfTuples(N); }
	if(ComputePtrel) { ptrel = rho->NewInstance(); ptrel->SetNumberOfComponents(1); ptrel->SetNumberOfTuples(N); }
	if(ComputeTtrel) { ttrel = rho->NewInstance(); ttrel->SetNumberOfComponents(1); ttrel->SetNumberOfTuples(N); }
	if(ComputeS) { s = rho->NewInstance(); s->SetNumberOfComponents(1); s->SetNumberOfTuples(N); }
	if(ComputeV) { v = rho->NewInstance(); v->SetNumberOfComponents(3); v->SetNumberOfTuples(N); }
	if(ComputeW) { w = rho->NewInstance(); w->SetNumberOfComponents(3); w->SetNumberOfTuples(N); }
	if(ComputeM) { M = rho->NewInstance(); M->SetNumberOfComponents(1); M->SetNumberOfTuples(N); }
	if(ComputeMrel) { Mrel = rho->NewInstance(); Mrel->SetNumberOfComponents(1); Mrel->SetNumberOfTuples(N); }
	if(ComputePhi) { phi = rho->NewInstance(); phi->SetNumberOfComponents(1); phi->SetNumberOfTuples(N); }
	if(ComputeAlpha) { alpha = rho->NewInstance(); alpha->SetNumberOfComponents(1); alpha->SetNumberOfTuples(N); }
	if(ComputeBeta) { beta = rho->NewInstance(); beta->SetNumberOfComponents(1); beta->SetNumberOfTuples(N); }
	if(ComputeMu) { mu = rho->NewInstance(); mu->SetNumberOfComponents(1); mu->SetNumberOfTuples(N); }
	
	switch(vtkTemplate2PackMacro(coords->GetDataType(), rho->GetDataType()))
	{
		vtkTemplate2Macro((this->CallSMPReconstructOp<VTK_T1, VTK_T2>(
			coords, rho, rhov, rhoe,
			ps, ts, ptabs, ttabs, ptrel, ttrel, s, v, w, M, Mrel, phi, alpha, beta, mu,
			ComputePs, ComputeTs, ComputePtabs, ComputeTtabs, ComputePtrel,
			ComputeTtrel, ComputeS, ComputeV, ComputeW, ComputeM, ComputeMrel, ComputePhi,
			ComputeAlpha, ComputeBeta, ComputeMu,
			gamma, rgas, omega, N)));
		default:
			vtkErrorMacro("Should not happen");
			return 0;
			break;
	}
	if(ps) { this->ReplaceVarByArrayIfNeeded(outData, ps, PS); ps->Delete(); }
	if(ts) { this->ReplaceVarByArrayIfNeeded(outData, ts, TS); ts->Delete(); }
	if(ptabs) { this->ReplaceVarByArrayIfNeeded(outData, ptabs, PTABS); ptabs->Delete(); }
	if(ttabs) { this->ReplaceVarByArrayIfNeeded(outData, ttabs, TTABS); ttabs->Delete(); }
	if(ptrel) { this->ReplaceVarByArrayIfNeeded(outData, ptrel, PTREL); ptrel->Delete(); }
	if(ttrel) { this->ReplaceVarByArrayIfNeeded(outData, ttrel, TTREL); ttrel->Delete(); }
	if(s) { this->ReplaceVarByArrayIfNeeded(outData, s, ENTROPY); s->Delete(); }
	if(v) { this->ReplaceVarByArrayIfNeeded(outData, v, VELOCITYABS); v->Delete(); }
	if(w) { this->ReplaceVarByArrayIfNeeded(outData, w, VELOCITYREL); w->Delete(); }
	if(M) { this->ReplaceVarByArrayIfNeeded(outData, M, MACHABS); M->Delete(); }
	if(Mrel) { this->ReplaceVarByArrayIfNeeded(outData, Mrel, MACHREL); Mrel->Delete(); }
	if(phi) { this->ReplaceVarByArrayIfNeeded(outData, phi, PHI); phi->Delete(); }
	if(alpha) { this->ReplaceVarByArrayIfNeeded(outData, alpha, ALPHA); alpha->Delete(); }
	if(beta) { this->ReplaceVarByArrayIfNeeded(outData, beta, BETA); beta->Delete(); }
	if(mu) { this->ReplaceVarByArrayIfNeeded(outData, mu, MU); mu->Delete(); }
	
	return 1;
}

void vtkIzarReconstructFields::ReplaceVarByArrayIfNeeded(vtkDataSet* data, vtkDataArray* ar, const char* name)
{
	if(data->GetPointData()->HasArray(name))
	{
		if(!this->ReplaceVars) { return; }
		data->GetPointData()->RemoveArray(name);
	}
	ar->SetName(name);
	data->GetPointData()->AddArray(ar);
}
