#include "vtkIzarGenericXRCut.h"
#include "IzarHelpers.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointSet.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkMath.h"
#include "vtkSMPTools.h"
#include "vtkAppendPolyData.h"
#include "vtkPlane.h"
#include "vtkCylinder.h"
#include "vtkClipPolyData.h"
#include "vtkCutter.h"

#include <sstream>

vtkStandardNewMacro(vtkIzarGenericXRCut)
vtkStandardNewMacro(vtkIzarConeXrFunction)

vtkIzarConeXrFunction::vtkIzarConeXrFunction():
	Xr1{0., 0.}, Xr2{0., 1.}
{
}

double vtkIzarConeXrFunction::EvaluateFunction(double pt[3])
{
    double r = sqrt(pt[1]*pt[1] + pt[2]*pt[2]);
    return (Xr1[1]-Xr2[1])*(pt[0]-Xr1[0]) + (Xr2[0]-Xr1[0])*(r - Xr1[1]);
}

void vtkIzarConeXrFunction::EvaluateGradient(double pt[3], double g[3])
{
    double r = sqrt(pt[1]*pt[1] + pt[2]*pt[2]);
    g[0] = Xr1[1]-Xr2[1];
    g[1] = (Xr2[0]-Xr1[0])*pt[1]/r;
    g[2] = (Xr2[0]-Xr1[0])*pt[2]/r;
}

void vtkIzarConeXrFunction::PrintSelf( ostream& os, vtkIndent indent )
{
    this->Superclass::PrintSelf( os, indent );
}

vtkIzarGenericXRCut::vtkIzarGenericXRCut():
	XList(), RList(), ClipFirstPoint(1), ClipLastPoint(1), ComputeNormals(1), Parametrize(1)
{
	IZAR_WARNING
}

vtkIzarGenericXRCut::~vtkIzarGenericXRCut()
{
}

void vtkIzarGenericXRCut::PrintSelf(ostream& os, vtkIndent indent)
{
    os << indent << "vtkIzarGenericXRCut";
}

int vtkIzarGenericXRCut::RequestData(vtkInformation* request,
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
		vtkErrorMacro("Error while retrieivng the input object. Make sure that this is a PointSet.");
		return 0;
	}
	
	vtkInformation* outInfo = outVector->GetInformationObject(0);
	vtkPolyData* outData = vtkPolyData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!outData)
	{
		return 0;
	}
	
	outData->Initialize();
	
	/* Perform some checks */
	if(this->XList.size() != this->RList.size())
	{
		vtkErrorMacro("X and R size mismatch");
		return 0;
	}
	if(this->XList.size() < 2)
	{
		vtkErrorMacro("Not enough points in the XR lists");
		return 0;
	}
	
	/* Perform the computation through a worker object */
	vtkIzarSMPXRCutter worker;
	worker.input = inData;
	worker.output = std::vector<vtkPolyData*>(this->XList.size()-1, NULL);
	worker.XList = &(this->XList);
	worker.RList = &(this->RList);
	worker.ClipFirstPoint = this->ClipFirstPoint;
	worker.ClipLastPoint = this->ClipLastPoint;
	worker.ComputeNormals = this->ComputeNormals;
	worker.Parametrize = this->Parametrize;
	if(this->Parametrize)
	{
		worker.InitParametrization();
	}
	vtkSMPTools::For(0, this->XList.size()-1, worker);
	
	/* Merge the resulting polydatas */
	vtkAppendPolyData* merger = vtkAppendPolyData::New();
	bool allEmpty = true;
	for(int ii = 0; ii < this->XList.size()-1; ii++)
	{
		// Check that we are not appending empty data
		if(worker.output[ii] && (worker.output[ii]->GetNumberOfPoints() != 0))
		{
			merger->AddInputData(worker.output[ii]);
			allEmpty = false;
		}
	}
	if(!allEmpty)
	{
		merger->Update();
		outData->ShallowCopy(merger->GetOutput());
	}
	
	merger->Delete();
	for(int ii = 0; ii < this->XList.size()-1; ii++)
	{
		if(worker.output[ii])
		{
			worker.output[ii]->Delete();
		}
	}
	return 1;
}

void vtkIzarSMPXRCutter::operator() (vtkIdType begin, vtkIdType end)
{
	int nPts = XList->size();
	int nCones = nPts-1;
	
	/* Implicit function to perform the slices */
	vtkIzarConeXrFunction* coneFunction = vtkIzarConeXrFunction::New();
	
	/* Implicit functions to perform the clippings */
	vtkPlane* planeMin = vtkPlane::New();
	planeMin->SetNormal(1., 0., 0.);
	vtkPlane* planeMax = vtkPlane::New();
	planeMax->SetNormal(1., 0., 0.);
	vtkCylinder* cylinderMin = vtkCylinder::New();
	cylinderMin->SetCenter(0., 0., 0.);
	cylinderMin->SetAxis(1., 0., 0.);
	vtkCylinder* cylinderMax = vtkCylinder::New();
	cylinderMax->SetCenter(0., 0., 0.);
	cylinderMax->SetAxis(1., 0., 0.);
	
	for(vtkIdType ii = begin; ii < end; ii++)
	{
		double x1 = (*XList)[ii];
		double x2 = (*XList)[ii+1];
		double r1 = (*RList)[ii];
		double r2 = (*RList)[ii+1];
		
		coneFunction->SetXr1(x1, r1);
		coneFunction->SetXr2(x2, r2);
		
		// Generates a mini pipeline to slice then optionally clip
		// Cutter
		vtkCutter* cutter = vtkCutter::New();
		vtkClipPolyData* clipperMin = NULL;
		vtkClipPolyData* clipperMax = NULL;
		cutter->SetCutFunction(coneFunction);
		cutter->SetNumberOfContours(1);
		cutter->SetValue(0, 0.0);
		cutter->SetInputData(input);
		vtkPolyDataAlgorithm* lastFilter = cutter;
		
		// First clip
		if((ii > 0) || ((ii == 0) && (ClipFirstPoint)))
		{
			clipperMin = vtkClipPolyData::New();
			clipperMin->SetInputConnection(lastFilter->GetOutputPort());
			// Detects whether it is better to clip along a plane or a cylinder
			if(fabs(x2-x1) > fabs(r2-r1))
			{
				// Clip along a plane
				planeMin->SetOrigin(x1, 0., 0.);
				clipperMin->SetClipFunction(planeMin);
				if(x2 > x1)
				{
					clipperMin->InsideOutOff(); // Keeps x >= xref
				}
				else
				{
					clipperMin->InsideOutOn(); // Keeps x <= xref
				}
			}
			else
			{
				// Clip along a cylinder
				cylinderMin->SetRadius(r1);
				clipperMin->SetClipFunction(cylinderMin);
				if(r2 > r1)
				{
					clipperMin->InsideOutOff();
				}
				else
				{
					clipperMin->InsideOutOn();
				}
			}
			lastFilter = clipperMin;
		}
		
		// Second clip
		if((ii < nCones-1) || ((ii == nCones-1) && (ClipLastPoint)))
		{
			clipperMax = vtkClipPolyData::New();
			clipperMax->SetInputConnection(lastFilter->GetOutputPort());
			
			if(fabs(x2-x1) > fabs(r2-r1))
			{
				// Plane clip
				planeMax->SetOrigin(x2, 0., 0.);
				clipperMax->SetClipFunction(planeMax);
				if(x2 > x1)
				{
					clipperMax->InsideOutOn();
				}
				else
				{
					clipperMax->InsideOutOff();
				}
			}
			else
			{
				// Cylinder clip
				cylinderMax->SetRadius(r2);
				clipperMax->SetClipFunction(cylinderMax);
				if(r2 > r1)
				{
					clipperMax->InsideOutOn();
				}
				else
				{
					clipperMax->InsideOutOff();
				}
			}
			lastFilter = clipperMax;
		}
		
		lastFilter->Update();
		output[ii] = lastFilter->GetOutput();
		// Take ownership of the newly generated object, as we are going to destroy the pipeline
		output[ii]->Register(NULL);
		
		// Clean the pipeline
		if(clipperMax)
		{
			clipperMax->Delete();
		}
		if(clipperMin)
		{
			clipperMin->Delete();
		}
		cutter->Delete();
		
		if(output[ii] && ComputeNormals)
		{
			// Computes the normals at the nodes and at the cell centers
			vtkIdType nDataPts = output[ii]->GetNumberOfPoints();
			vtkIdType nCells = output[ii]->GetNumberOfCells();
			vtkDoubleArray* normalsPts = vtkDoubleArray::New();
			normalsPts->SetName("Normals");
			normalsPts->SetNumberOfComponents(3);
			normalsPts->SetNumberOfTuples(nDataPts);
			double* normalsPtsPtr = reinterpret_cast<double*>(normalsPts->GetVoidPointer(0));
			vtkDoubleArray* normalsCells = vtkDoubleArray::New();
			normalsCells->SetName("Normals");
			normalsCells->SetNumberOfComponents(3);
			normalsCells->SetNumberOfTuples(nCells);
			double* normalsCellsPtr = reinterpret_cast<double*>(normalsCells->GetVoidPointer(0));
			double curPt[3];
			
			// Normals at the points
			for(vtkIdType kk = 0; kk < nDataPts; kk++)
			{
				output[ii]->GetPoint(kk, curPt);
				coneFunction->EvaluateGradient(curPt, &(normalsPtsPtr[3*kk]));
				vtkMath::Normalize(&(normalsPtsPtr[3*kk]));
			}
			
			// Normals at the cells
			vtkIdType nCellPts(0);
			const vtkIdType* cellPts(nullptr);
			vtkIdType idCell = 0;
			vtkCellArray* ca = output[ii]->GetPolys();
			ca->InitTraversal();
			while(ca->GetNextCell(nCellPts, cellPts))
			{
				// Get the cell barycenter
				double center[3] = { 0., 0., 0. };
				for(vtkIdType idPt = 0; idPt < nCellPts; idPt++)
				{
					output[ii]->GetPoint(cellPts[idPt], curPt);
					for(int dim = 0; dim < 3; dim++)
					{
						center[dim] += curPt[dim]/((double)(nCellPts));
					}
				}
				 
				coneFunction->EvaluateGradient(center, &(normalsCellsPtr[3*idCell]));
				vtkMath::Normalize(&(normalsCellsPtr[3*idCell]));
				idCell++;
			}
			output[ii]->GetPointData()->SetNormals(normalsPts);
			output[ii]->GetCellData()->SetNormals(normalsCells);
			
			normalsPts->Delete();
			normalsCells->Delete();
		}
		
		if(output[ii] && Parametrize)
		{
			vtkIdType nDataPts = output[ii]->GetNumberOfPoints();
			vtkDoubleArray* parameter = vtkDoubleArray::New();
			parameter->SetNumberOfComponents(1);
			parameter->SetNumberOfTuples(nDataPts);
			double* parameterPtr = reinterpret_cast<double*>(parameter->GetVoidPointer(0));
			double curPt[3];
			for(vtkIdType kk = 0;  kk < nDataPts; kk++)
			{
				output[ii]->GetPoint(kk, curPt);
				double curR = sqrt(curPt[1]*curPt[1] + curPt[2]*curPt[2]);
				parameterPtr[kk] = (IntermediateLengths[ii] + sqrt((curPt[0]-x1)*(curPt[0]-x1) + (curR-r1)*(curR-r1)))/TotalLength;
			}
			parameter->SetName("XRParameter");
			if(output[ii]->GetPointData()->HasArray("XRParameter"))
			{
				output[ii]->GetPointData()->RemoveArray("XRParameter");
			}
			output[ii]->GetPointData()->AddArray(parameter);
			parameter->Delete();
		}
	}
	
	coneFunction->Delete();
	planeMin->Delete();
	planeMax->Delete();
	cylinderMin->Delete();
	cylinderMax->Delete();
}

void vtkIzarSMPXRCutter::InitParametrization()
{
	int n = this->XList->size();
	this->IntermediateLengths.resize(n-1);
	this->TotalLength = 0.0;
	for(int ii = 0; ii < n-1; ii++)
	{
		double x1 = (*(this->XList))[ii];
		double x2 = (*(this->XList))[ii+1];
		double r1 = (*(this->RList))[ii];
		double r2 = (*(this->RList))[ii+1];
		this->IntermediateLengths[ii] = this->TotalLength;
		this->TotalLength += sqrt((x2-x1)*(x2-x1) + (r2-r1)*(r2-r1));
	}
}
