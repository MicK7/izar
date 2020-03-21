
#include "vtkIzarAddFieldDataPG.h"

#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkIzarAddFieldDataPG)

vtkIzarAddFieldDataPG::vtkIzarAddFieldDataPG() : vtkIzarAddFieldData()
{
	this->intDataNames.resize(1);
	this->intDataValues.resize(1);
	this->doubleDataNames.resize(3);
	this->doubleDataValues.resize(3);
	
	this->intDataNames[0] = ZSECTOR;
	this->doubleDataNames[0] = GAMMA;
	this->doubleDataNames[1] = RGAS;
	this->doubleDataNames[2] = OMEGA;
	
	IZAR_WARNING
}

vtkIzarAddFieldDataPG::~vtkIzarAddFieldDataPG()
{
}

void vtkIzarAddFieldDataPG::PrintSelf(ostream& os, vtkIndent indent)
{
	os << "vtkIzarAddFieldDataPG\n";
	this->Superclass::PrintSelf(os, indent.GetNextIndent());
}
