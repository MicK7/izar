#ifndef VTKIZARCHANGEFRAME_H
#define VTKIZARCHANGEFRAME_H

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
class vtkSMPIzarChangeFrameOp
{
public:
	/* Inputs */
	CoordsType* coords;
	ValueType* ts;
	ValueType* ps;
	ValueType* v;
	
	/* Outputs */
	ValueType* ttrel;
	ValueType* ptrel;
	ValueType* w;
	ValueType* Mrel;
	ValueType* beta;
	
	/* Constants */
	double omega;
	
	double gamma;
	double rgas;
	
	void operator() (vtkIdType begin, vtkIdType end)
	{
		for(vtkIdType ii = begin; ii != end; ii++)
		{
			ValueType r = sqrt(coords[3*ii+1]*coords[3*ii+1] + coords[3*ii+2]*coords[3*ii+2]);
			ValueType theta = atan2(coords[3*ii+2], coords[3*ii+1]);
			ValueType cost = cos(theta);
			ValueType sint = sin(theta);
			ValueType vr = v[3*ii+1]*cost + v[3*ii+2]*sint;
			ValueType vt = -v[3*ii+1]*sint + v[3*ii+2]*cost;
			ValueType wt = vt - omega*r;
			ValueType wm = sqrt(v[3*ii]*v[3*ii] + vr*vr);
			beta[ii] = atan2(wt, wm)*180./vtkMath::Pi();
			w[3*ii] = v[3*ii];
			w[3*ii+1] = vr*cost - wt*sint;
			w[3*ii+2] = vr*sint + wt*cost;
			ValueType magW2 = w[3*ii]*w[3*ii] + w[3*ii+1]*w[3*ii+1] + w[3*ii+2]*w[3*ii+2];
			ValueType Mrel2 = magW2/(gamma*rgas*ts[ii]);
			Mrel[ii] = sqrt(Mrel2);
			ValueType tt_ts_rel_ratio = 1. + 0.5*(gamma-1.) * Mrel2;
			ttrel[ii] = ts[ii]*tt_ts_rel_ratio;
			ptrel[ii] = ps[ii]*pow(tt_ts_rel_ratio, gamma/(gamma-1.));
		}
	}
};


/**
 * Computes the variables in another frame of reference.
 * Needs the variables computed by ReconstructVariables
 */
class VTK_EXPORT vtkIzarChangeFrame : public vtkPointSetAlgorithm
{
public:
    static vtkIzarChangeFrame* New();
    vtkTypeMacro(vtkIzarChangeFrame, vtkPointSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

	vtkGetMacro(RPMInsteadOfRadPerSecond, int)
	vtkSetMacro(RPMInsteadOfRadPerSecond, int)

	vtkGetMacro(FrameRotationSpeed, double)
	vtkSetMacro(FrameRotationSpeed, double)
	
	/**
     * Standard Setter
     */
	void SetFrameName(const char* name) { this->FrameName = name; this->Modified(); }
    /**
     * Standard Getter
     */
    const char* GetFrameName() { return this->FrameName.c_str(); }
	
protected:
    vtkIzarChangeFrame();
    ~vtkIzarChangeFrame();

    virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);	
	void NameAndAddArray(vtkPointSet* data, vtkDataArray* ar, const std::string& name);
	
	int RPMInsteadOfRadPerSecond;
	double FrameRotationSpeed;
	std::string FrameName;
	
	template <typename CoordsType, typename ValueType>
	void CallSMPChangeFrameOp(vtkDataArray* coords, vtkDataArray* ts, vtkDataArray* ps, vtkDataArray* v,
		vtkDataArray* ttrel, vtkDataArray* ptrel, vtkDataArray* w, vtkDataArray* Mrel, vtkDataArray* beta,
		double gamma, double rgas,
		int nPts)
	{
		vtkSMPIzarChangeFrameOp<CoordsType, ValueType> worker;
		
		worker.coords = reinterpret_cast<CoordsType*>(coords->GetVoidPointer(0));
		worker.ts = reinterpret_cast<ValueType*>(ts->GetVoidPointer(0));
		worker.ps = reinterpret_cast<ValueType*>(ps->GetVoidPointer(0));
		worker.v = reinterpret_cast<ValueType*>(v->GetVoidPointer(0));
		
		worker.ttrel = reinterpret_cast<ValueType*>(ttrel->GetVoidPointer(0));
		worker.ptrel = reinterpret_cast<ValueType*>(ptrel->GetVoidPointer(0));
		worker.w = reinterpret_cast<ValueType*>(w->GetVoidPointer(0));
		worker.Mrel = reinterpret_cast<ValueType*>(Mrel->GetVoidPointer(0));
		worker.beta = reinterpret_cast<ValueType*>(beta->GetVoidPointer(0));
		
		worker.gamma = gamma;
		worker.rgas = rgas;
		if(this->RPMInsteadOfRadPerSecond)
		{
			worker.omega = this->FrameRotationSpeed * 2.*vtkMath::Pi()/60.;
		}
		else
		{
			worker.omega = this->FrameRotationSpeed;
		}
		
		vtkSMPTools::For(0, nPts, worker);
	}
	
private:
    vtkIzarChangeFrame operator=(const vtkIzarChangeFrame&);
    vtkIzarChangeFrame(const vtkIzarChangeFrame&);
};

#endif
