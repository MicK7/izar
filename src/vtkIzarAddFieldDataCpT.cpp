
#include "vtkIzarAddFieldDataCpT.h"

#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkIzarAddFieldDataCpT)

vtkIzarAddFieldDataCpT::vtkIzarAddFieldDataCpT() : vtkIzarAddFieldData()
{
	this->intDataNames.resize(1);
	this->intDataValues.resize(1);
	this->doubleDataNames.resize(6);
	this->doubleDataValues.resize(6);
	
	this->intDataNames[0] = ZSECTOR;
	this->doubleDataNames[0] = RGAS;
	this->doubleDataNames[1] = RFUEL;
	this->doubleDataNames[2] = RWATER;
	this->doubleDataNames[3] = FAR;
	this->doubleDataNames[4] = WAR;
	this->doubleDataNames[5] = OMEGA;
	
	IZAR_WARNING
}

vtkIzarAddFieldDataCpT::~vtkIzarAddFieldDataCpT()
{
}

void vtkIzarAddFieldDataCpT::PrintSelf(ostream& os, vtkIndent indent)
{
	os << "vtkIzarAddFieldDataCpT\n";
	this->Superclass::PrintSelf(os, indent.GetNextIndent());
}
