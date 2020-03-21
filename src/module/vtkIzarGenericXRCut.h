#ifndef VTKIZARGENERICXRCUT_H
#define VTKIZARGENERICXRCUT_H

#include "IzarDefines.h"

#include "vtkSetGet.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkSMPTools.h"
#include "vtkMath.h"
#include "vtkDataArray.h"
#include "vtkImplicitFunction.h"



class vtkInformation;
class vtkInformationVector;
class vtkDataSetAttributes;
class vtkPointSet;
class vtkDataSet;

/* Implicit function to perform the cutting */
class VTK_EXPORT vtkIzarConeXrFunction : public vtkImplicitFunction
{
public:
    static vtkIzarConeXrFunction* New();
    vtkTypeMacro(vtkIzarConeXrFunction, vtkImplicitFunction);
    void PrintSelf(ostream& os, vtkIndent indent);
    
    double EvaluateFunction(double x[3]);
    double EvaluateFunction(double x, double y, double z)
        {return this->vtkImplicitFunction::EvaluateFunction(x, y, z); };
    void EvaluateGradient(double x[3], double g[3]);
    
    vtkGetVector2Macro(Xr1, double);
    vtkSetVector2Macro(Xr1, double);
    vtkGetVector2Macro(Xr2, double);
    vtkSetVector2Macro(Xr2, double);
    
protected:
    vtkIzarConeXrFunction();
    ~vtkIzarConeXrFunction() {};
    
    double Xr1[2];
    double Xr2[2];
    
private:
    vtkIzarConeXrFunction(const vtkIzarConeXrFunction&);
    void operator=(const vtkIzarConeXrFunction&); 
    
};

/**
 * Generates an oriented cut over a revolution surface defined by
 * an XR line. Not to be used directly as a ParaView plugin, but as a
 * base class to several interface plugins.
 * 
 * Performs individual XR cuts, then clips them if needed and assembles them.
 * Uses a SMP strategy to perform the cutting and clipping in parallel.
 */
class VTK_EXPORT vtkIzarGenericXRCut : public vtkPolyDataAlgorithm
{
public:
    static vtkIzarGenericXRCut* New();
    vtkTypeMacro(vtkIzarGenericXRCut, vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

	vtkGetMacro(ClipFirstPoint, int)
	vtkSetMacro(ClipFirstPoint, int)
	vtkGetMacro(ClipLastPoint, int)
	vtkSetMacro(ClipLastPoint, int)
	vtkGetMacro(ComputeNormals, int)
	vtkSetMacro(ComputeNormals, int)
	vtkGetMacro(Parametrize, int)
	vtkSetMacro(Parametrize, int)
	
	
protected:
    vtkIzarGenericXRCut();
    ~vtkIzarGenericXRCut();

    virtual int RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector);	
    
    std::vector<double> XList;
    std::vector<double> RList;
    
    int ClipFirstPoint;
    int ClipLastPoint;
    
    int ComputeNormals;
    int Parametrize;
    
private:
    vtkIzarGenericXRCut operator=(const vtkIzarGenericXRCut&);
    vtkIzarGenericXRCut(const vtkIzarGenericXRCut&);
};

/**
 * SMP worker to perform the cuts and the clips
 */
class vtkIzarSMPXRCutter
{
public:
	vtkPointSet* input;
	std::vector<vtkPolyData*> output;
	
	const std::vector<double>* XList;
	const std::vector<double>* RList;
	
	int ClipFirstPoint;
	int ClipLastPoint;
	
	int ComputeNormals;
	int Parametrize;
	
	/**
	 * Worker
	 * 
	 * Assumes that XList.size() == RList.size() == output.size(),
	 * that output is made of NULLs, and that XList.size() >= 2
	 */
	void operator() (vtkIdType begin, vtkIdType end);
	
	/**
	 * Computes the intermediate lengths so that each thread can
	 * parametrize its surfaces
	 */
	void InitParametrization();

protected:
	std::vector<double> IntermediateLengths;
	double TotalLength;
};

#endif
