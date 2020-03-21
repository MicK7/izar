#include "vtkIzarConeCut.h"

#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkIzarConeCut)

vtkIzarConeCut::vtkIzarConeCut(): vtkIzarGenericXRCut()
{
	this->XList.resize(2);
	this->XList[0] = 0.0;
	this->XList[1] = 0.0;
	this->RList.resize(2);
	this->RList[0] = 0.0;
	this->RList[1] = 1.0;
	
	this->ClipFirstPoint = 0;
	this->ClipLastPoint = 0;
}

vtkIzarConeCut::~vtkIzarConeCut()
{
}

void vtkIzarConeCut::PrintSelf(ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarConeCut\n";
}
