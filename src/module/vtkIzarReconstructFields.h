#ifndef VTKIZARRECONSTRUCTFIELDS_H
#define VTKIZARRECONSTRUCTFIELDS_H

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
 * Functor to compute the data with the vtkSMPTools framework
 */
template <typename CoordsValueType, typename ValueType>
class vtkSMPReconstructFieldsOp
{
public:
	// Inputs
	CoordsValueType* coords;
	ValueType* rho;
	ValueType* rhov;
	ValueType* rhoe;
	
	// Outputs
	ValueType* ps;
	ValueType* ts;
	ValueType* ptabs;
	ValueType* ttabs;
	ValueType* ptrel;
	ValueType* ttrel;
	ValueType* s;
	ValueType* v;
	ValueType* w;
	ValueType* M;
	ValueType* Mrel;
	ValueType* phi;
	ValueType* alpha;
	ValueType* beta;
	ValueType* mu;
	
	// Defines whether or not we have to compute the individual variables
	int computePs;
	int computeTs;
	int computePtabs;
	int computeTtabs;
	int computePtrel;
	int computeTtrel;
	int computeS;
	int computeV;
	int computeW;
	int computeM;
	int computeMrel;
	int computePhi;
	int computeAlpha;
	int computeBeta;
	int computeMu;
	
	// Constants
	double suthConst;
	double suthMuRef;
	double suthTRef;
	
	double cv;
	double cp;
	double rgas;
	double gamma;
	
	double omega;
	
	int nPts;
	
	/**
	 * Performs the actual computaiton of the values
	 */
	void operator()(vtkIdType begin, vtkIdType end)
	{
		for(vtkIdType ii = begin; ii != end; ii++)
		{
			// Velocities
			ValueType tmpW[3];
			ValueType tmpV[3];
			for(int dim = 0; dim < 3; dim++)
			{
				tmpW[dim] = rhov[ii*3+dim]/rho[ii];
			}
			ValueType tmpR = sqrt(coords[ii*3+1]*coords[ii*3+1] + coords[ii*3+2]*coords[ii*3+2]);
			ValueType tmpTheta = atan2(coords[ii*3+2], coords[ii*3+1]);
			ValueType tmpCos = cos(tmpTheta);
			ValueType tmpSin = sin(tmpTheta);
			ValueType tmpMagW2 = tmpW[0]*tmpW[0] + tmpW[1]*tmpW[1] + tmpW[2]*tmpW[2];
			ValueType tmpWr = tmpW[1]*tmpCos + tmpW[2]*tmpSin;
			ValueType tmpWt = -tmpW[1]*tmpSin + tmpW[2]*tmpCos;
			ValueType tmpVt = tmpWt + tmpR*omega;
			tmpV[0] = tmpW[0];
			tmpV[1] = tmpCos*tmpWr - tmpSin*tmpVt;
			tmpV[2] = tmpSin*tmpWr + tmpCos*tmpVt;
			ValueType tmpMagV2 = tmpV[0]*tmpV[0] + tmpV[1]*tmpV[1] + tmpV[2]*tmpV[2];
			
			// Angles
			ValueType tmpWm = sqrt(tmpW[0]*tmpW[0] + tmpWr*tmpWr);
			ValueType tmpPhi = atan2(tmpWr, tmpW[0]);
			ValueType tmpAlpha = vtkMath::DegreesFromRadians(atan2(tmpVt, tmpWm));
			ValueType tmpBeta = vtkMath::DegreesFromRadians(atan2(tmpWt, tmpWm));
			
			// Static pressure and temperature
			ValueType tmpTs = (rhoe[ii]/rho[ii] - 0.5*tmpMagW2)/cv;
			ValueType tmpPs = rho[ii]*rgas*tmpTs;
			ValueType tmpS = cp*log(tmpTs/TREF) - rgas*log(tmpPs/PREF);
			
			ValueType tmpMu = suthMuRef*sqrt(tmpTs/suthTRef)*(1.+suthConst/suthTRef)/(1.+suthConst/tmpTs);
			
			// Relative then absolute total values
			ValueType tmpA = sqrt(gamma*rgas*tmpTs);
			ValueType tmpMrel = sqrt(tmpMagW2)/tmpA;
			ValueType tempRatio_rel = (1.+0.5*(gamma-1.)*tmpMrel*tmpMrel);
			ValueType tmpTtrel = tmpTs*tempRatio_rel;
			ValueType tmpPtrel = tmpPs*pow(tempRatio_rel, gamma/(gamma-1.));
			
			ValueType tmpMabs = sqrt(tmpMagV2)/tmpA;
			ValueType tempRatio_abs = (1.+0.5*(gamma-1.)*tmpMabs*tmpMabs);
			ValueType tmpTtabs = tmpTs*tempRatio_abs;
			ValueType tmpPtabs = tmpPs*pow(tempRatio_abs, gamma/(gamma-1.));
			
			if(computePs) { ps[ii] = tmpPs; }
			if(computeTs) { ts[ii] = tmpTs; }
			if(computePtabs) { ptabs[ii] = tmpPtabs; }
			if(computeTtabs) { ttabs[ii] = tmpTtabs; }
			if(computePtrel) { ptrel[ii] = tmpPtrel; }
			if(computeTtrel) { ttrel[ii] = tmpTtrel; }
			if(computeS) { s[ii] = tmpS; }
			if(computeV) { for(int dim = 0; dim < 3; dim++) { v[ii*3+dim] = tmpV[dim]; } }
			if(computeW) { for(int dim = 0; dim < 3; dim++) { w[ii*3+dim] = tmpW[dim]; } }
			if(computeM) { M[ii] = tmpMabs; }
			if(computeMrel) { Mrel[ii] = tmpMrel; }
			if(computePhi) { phi[ii] = tmpPhi; }
			if(computeAlpha) { alpha[ii] = tmpAlpha; }
			if(computeBeta) { beta[ii] = tmpBeta; }
			if(computeMu) { mu[ii] = tmpMu; }
		}
	}
};

class VTK_EXPORT vtkIzarReconstructFields : public vtkPointSetAlgorithm
{
public:
    static vtkIzarReconstructFields* New();
    vtkTypeMacro(vtkIzarReconstructFields, vtkPointSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);


	vtkGetMacro(ReplaceVars, int)
	vtkSetMacro(ReplaceVars, int)
	vtkGetMacro(AdvancedMode, int)
	vtkSetMacro(AdvancedMode, int)

	vtkGetMacro(SuthConst, double)
	vtkSetMacro(SuthConst, double)
	vtkGetMacro(SuthMuRef, double)
	vtkSetMacro(SuthMuRef, double)
	vtkGetMacro(SuthTRef, double)
	vtkSetMacro(SuthTRef, double)
	
	vtkSetMacro(ComputePs, int)
	vtkSetMacro(ComputeTs, int)
	vtkSetMacro(ComputePtabs, int)
	vtkSetMacro(ComputeTtabs, int)
	vtkSetMacro(ComputePtrel, int)
	vtkSetMacro(ComputeTtrel, int)
	vtkSetMacro(ComputeS, int)
	vtkSetMacro(ComputeV, int)
	vtkSetMacro(ComputeW, int)
	vtkSetMacro(ComputeM, int)
	vtkSetMacro(ComputeMrel, int)
	vtkSetMacro(ComputePhi, int)
	vtkSetMacro(ComputeAlpha, int)
	vtkSetMacro(ComputeBeta, int)
	vtkSetMacro(ComputeMu, int)
	
protected:
    vtkIzarReconstructFields();
    ~vtkIzarReconstructFields();

    virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);	
	
	void ReplaceVarByArrayIfNeeded(vtkDataSet* data, vtkDataArray* ar, const char* name);
	
	/**
	 * If true, if existing variables already exist, they are replaced
	 * by the new ones. True by default
	 */
	int ReplaceVars;
	
	/**
	 * Does nothing, just used to hide the "advanced" options
	 */
	int AdvancedMode;
	
	/**
	 * Variables available in the "advanced" mode
	 */
	double SuthConst;
	double SuthMuRef;
	double SuthTRef;
	
	int ComputePs;
	int ComputeTs;
	int ComputePtabs;
	int ComputeTtabs;
	int ComputePtrel;
	int ComputeTtrel;
	int ComputeS;
	int ComputeV;
	int ComputeW;
	int ComputeM;
	int ComputeMrel;
	int ComputePhi;
	int ComputeAlpha;
	int ComputeBeta;
	int ComputeMu;
	
	
	/** Worker template: generates the vtkSMPTools functor and calls it
	 */
	template<typename CoordsValueType, typename ValueType>
	void CallSMPReconstructOp(
		vtkDataArray* coords, vtkDataArray* rho, vtkDataArray* rhov, vtkDataArray* rhoe,
		
		vtkDataArray* ps, vtkDataArray* ts, vtkDataArray* ptabs, vtkDataArray* ttabs,
		vtkDataArray* ptrel, vtkDataArray* ttrel, vtkDataArray* s, vtkDataArray* v, vtkDataArray* w,
		vtkDataArray* M, vtkDataArray* Mrel, vtkDataArray* phi, vtkDataArray* alpha, vtkDataArray* beta, vtkDataArray* mu,
		
		int computePs, int computeTs, int computePtabs, int computeTtabs,
		int computePtrel, int computeTtrel, int computeS, int computeV, int computeW,
		int computeM, int computeMrel, int computePhi, int computeAlpha, int computeBeta, int computeMu,
		
		double gamma, double rgas, double omega, int nPts)
	{
		vtkSMPReconstructFieldsOp<CoordsValueType, ValueType> worker;
		worker.coords = reinterpret_cast<CoordsValueType*>(coords->GetVoidPointer(0));
		worker.rho = reinterpret_cast<ValueType*>(rho->GetVoidPointer(0));
		worker.rhov = reinterpret_cast<ValueType*>(rhov->GetVoidPointer(0));
		worker.rhoe = reinterpret_cast<ValueType*>(rhoe->GetVoidPointer(0));
		
		if(computePs)
		{
			worker.ps = reinterpret_cast<ValueType*>(ps->GetVoidPointer(0));
		}
		else
		{
			worker.ps = NULL;
		}
		if(computeTs)
		{
			worker.ts = reinterpret_cast<ValueType*>(ts->GetVoidPointer(0));
		}
		else
		{
			worker.ts = NULL;
		}
		if(computePtabs)
		{
			worker.ptabs = reinterpret_cast<ValueType*>(ptabs->GetVoidPointer(0));
		}
		else
		{
			worker.ptabs = NULL;
		}
		if(computeTtabs)
		{
			worker.ttabs = reinterpret_cast<ValueType*>(ttabs->GetVoidPointer(0));
		}
		else
		{
			worker.ttabs = NULL;
		}
		if(computePtrel)
		{
			worker.ptrel = reinterpret_cast<ValueType*>(ptrel->GetVoidPointer(0));
		}
		else
		{
			worker.ptrel = NULL;
		}
		if(computeTtrel)
		{
			worker.ttrel = reinterpret_cast<ValueType*>(ttrel->GetVoidPointer(0));
		}
		else
		{
			worker.ttrel = NULL;
		}
		if(computeS)
		{
			worker.s = reinterpret_cast<ValueType*>(s->GetVoidPointer(0));
		}
		else
		{
			worker.s = NULL;
		}
		if(computeV)
		{
			worker.v = reinterpret_cast<ValueType*>(v->GetVoidPointer(0));
		}
		else
		{
			worker.v = NULL;
		}
		if(computeW)
		{
			worker.w = reinterpret_cast<ValueType*>(w->GetVoidPointer(0));
		}
		else
		{
			worker.w = NULL;
		}
		if(computeM)
		{
			worker.M = reinterpret_cast<ValueType*>(M->GetVoidPointer(0));
		}
		else
		{
			worker.M = NULL;
		}
		if(computeMrel)
		{
			worker.Mrel = reinterpret_cast<ValueType*>(Mrel->GetVoidPointer(0));
		}
		else
		{
			worker.Mrel = NULL;
		}
		if(computePhi)
		{
			worker.phi = reinterpret_cast<ValueType*>(phi->GetVoidPointer(0));
		}
		else
		{
			worker.phi = NULL;
		}
		if(computeAlpha)
		{
			worker.alpha = reinterpret_cast<ValueType*>(alpha->GetVoidPointer(0));
		}
		else
		{
			worker.alpha = NULL;
		}
		if(computeBeta)
		{
			worker.beta = reinterpret_cast<ValueType*>(beta->GetVoidPointer(0));
		}
		else
		{
			worker.beta = NULL;
		}
		if(computeMu)
		{
			worker.mu = reinterpret_cast<ValueType*>(mu->GetVoidPointer(0));
		}
		else
		{
			worker.mu = NULL;
		}
		worker.suthConst = this->SuthConst;
		worker.suthMuRef = this->SuthMuRef;
		worker.suthTRef = this->SuthTRef;
		
		worker.computePs = computePs;
		worker.computeTs = computeTs;
		worker.computePtabs = computePtabs;
		worker.computeTtabs = computeTtabs;
		worker.computePtrel = computePtrel;
		worker.computeTtrel = computeTtrel;
		worker.computeS = computeS;
		worker.computeV = computeV;
		worker.computeW = computeW;
		worker.computeM = computeM;
		worker.computeMrel = computeMrel;
		worker.computePhi = computePhi;
		worker.computeAlpha = computeAlpha;
		worker.computeBeta = computeBeta;
		worker.computeMu = computeMu;
		
		worker.cv = rgas/(gamma-1.);
		worker.cp = rgas*gamma/(gamma-1.);
		worker.rgas = rgas;
		worker.gamma = gamma;
		worker.omega = omega;
		worker.nPts = nPts;
		
		vtkSMPTools::For(0, nPts, worker);
	}
	
private:
    vtkIzarReconstructFields operator=(const vtkIzarReconstructFields&);
    vtkIzarReconstructFields(const vtkIzarReconstructFields&);

};



#endif
