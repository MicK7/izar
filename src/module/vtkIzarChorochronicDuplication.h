#ifndef VTKIZARCHOROCHRONICDUPLICATION_H
#define VTKIZARCHOROCHRONICDUPLICATION_H

#include "IzarDefines.h"
#include "IzarHelpers.h"

#include "vtkSetGet.h"
#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkTransform.h"

#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>

class vtkInformation;
class vtkInformationVector;

/**
 * Duplicates a dataset by applying a temporal shift
 * 
 * Some information on the implementation : the request of several time steps
 * in a vtk pipeline is done by performing successive calls to RequestUpdateExtent
 * and RequestData : the RequestUpdateExtent is used to give the requested timestep, and
 * RequestData to retrieve this timestep. If we are not done yet and want to request
 * another time step, we have to add the CONTINUE_EXECUTING key to the request.
 * 
 * Optionnally, this filter can cache the data to avoid too much disk IO.
 * The user has to manually tune it so that it does not run out of RAM...
 */
class VTK_EXPORT vtkIzarChorochronicDuplication : public vtkMultiBlockDataSetAlgorithm
{
public:
    static vtkIzarChorochronicDuplication* New();
    vtkTypeMacro(vtkIzarChorochronicDuplication, vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    
    vtkGetMacro(ModeRotationSpeed, double)
    vtkSetMacro(ModeRotationSpeed, double)
    vtkGetMacro(RotationSpeedUnit, int)
    vtkSetMacro(RotationSpeedUnit, int)
    vtkGetMacro(NumberOfDiameters, int)
    vtkSetMacro(NumberOfDiameters, int)
    vtkGetMacro(ZSector, int)
    vtkSetMacro(ZSector, int)
    vtkGetMacro(IndexMin, int)
    vtkSetMacro(IndexMin, int)
    vtkGetMacro(IndexMax, int)
    vtkSetMacro(IndexMax, int)
    vtkGetMacro(TimeUnitInSeconds, double)
    vtkSetMacro(TimeUnitInSeconds, double)
    
protected:
	vtkIzarChorochronicDuplication();
	virtual ~vtkIzarChorochronicDuplication();
	virtual int RequestDataObject(vtkInformation* request,
		vtkInformationVector** inputVector, vtkInformationVector* outputVector);
	virtual int RequestInformation(vtkInformation* request,
		vtkInformationVector** inputVector, vtkInformationVector* outputVector);
	virtual int RequestUpdateExtent(vtkInformation* request,
		vtkInformationVector** inputVector, vtkInformationVector* outputVector);
	virtual int RequestData(vtkInformation* request,
		vtkInformationVector** inputVector, vtkInformationVector* outputVector);
	virtual int FillInputPortInformation(int port, vtkInformation* info);
	
	/**
	 * Check that the input data and the parameters provide a correct
	 * setup for a chorochronic reconstruction, and computes some
	 * important variables
	 */
	void SetAndCheckTemporalData();
	/**
	 * Called at each new time step : determines which time steps are
	 * used for the reconstruction and populates the DuplicationsTimesteps
	 * and TimestepsRequests vectors
	 */
	void InitializeRequest(double timestep);
	
	inline bool TimestepsAreEqual(const double& a, const double& b) { return fabs(a-b) < this->TimeTolerance; }
	
	// Parameters
	// Characteristics of the unsteady fluctuations
	/**
	 * Rotation speed in the current dataset frame of reference of the
	 * rotating wave driving the unsteadiness
	 * 
	 * In the case of a stage computation, this is the opposite row
	 * relative rotation speed.
	 */
	double ModeRotationSpeed;
	/**
	 * Boolean : true if the rotation speed is given in rpm, false
	 * if it is given in rad/s
	 */
	int RotationSpeedUnit;
	/**
	 * Number of diameters of the rotating wave driving the unsteadiness
	 * 
	 * In the case of a stage computation, this is the opposite number
	 * of blades
	 */
	int NumberOfDiameters;
	/**
	 * Periodicity of the simulated sector (not retrieved from
	 * the ZSector FieldData from technical reasons : we do not have the
	 * input dataset at the REQUEST_INFORMATION step)
	 */
	int ZSector;
	
	// Duplication parameters
	/**
	 * We are duplicating from IndexMin to IndexMax
	 */
	int IndexMin;
	/**
	 * We are duplicating from IndexMin to IndexMax
	 */
	int IndexMax;
	
	// Unit conversion parameters
	/**
	 * Duration of a VTK time unit in seconds
	 */
	double TimeUnitInSeconds;
	
	// Internal data. Most of them are set up during the RequestInformation pass.
	/**
	 * Rotation speed in rad/time unit
	 */
	double RotationSpeedDim;
	/**
	 * Stores the available timesteps, ordered
	 */
	std::vector<double> InputTimesteps;
	/**
	 * Stores the output timesteps, ordered. This is a subset of the InputTimesteps.
	 */
	std::vector<double> OutputTimesteps;
	/**
	 * Time period between two consecutive time steps
	 */
	double TimestepsPeriod;
	/**
	 * f(theta, t) = f(theta + DeltaTheta, t + DeltaTChoro)
	 */
	double DeltaTChoro;
	/**
	 * f(theta, t) = f(theta + DeltaTheta, t + DeltaTChoro)
	 */
	double DeltaTheta;
	/**
	 * Temporal period used in the modulo function
	 */
	double TemporalPeriod;
	/**
	 * Number of timesteps per period
	 */
	int NTimestepsPerPeriod;
	/**
	 * Signed number of timesteps representing a chorochronic period
	 */
	int NTimestepsPerChoroPeriod;
	/**
	 * Tolerance used to compare two timesteps. Dynamically computed
	 * from the time steps provided as an input. It is 1.e-3 times
	 * the period between two timesteps.
	 */
	double TimeTolerance;
	
	// To handle the retrieval of several time steps
	/**
	 * Timesteps requested for each duplication. May contain the same timestep several times
	 */
	std::map<int, double> DuplicationsTimesteps;
	/**
	 * Vector of unique timesteps to request, corresponding to all the timesteps represented in DuplicationsTimesteps
	 */
	std::vector<double> TimestepsRequests;
	/**
	 * Index corresponding to the currently requested timestep
	 * -1 if we are done
	 */
	int CurrentRequestIndex;
	
	// To handle a RequestInformation() fail
	bool TimestepsAreOK;
	std::string InfoErrorMsg;
	
private:
    vtkIzarChorochronicDuplication operator=(const vtkIzarChorochronicDuplication&);
    vtkIzarChorochronicDuplication(const vtkIzarChorochronicDuplication&);
};

#endif
