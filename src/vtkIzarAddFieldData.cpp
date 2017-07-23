
#include "vtkIzarAddFieldData.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"


vtkStandardNewMacro(vtkIzarAddFieldData)

vtkIzarAddFieldData::vtkIzarAddFieldData()
{
}

vtkIzarAddFieldData::~vtkIzarAddFieldData()
{
}

void vtkIzarAddFieldData::PrintSelf(ostream& os, vtkIndent indent)
{
	os << indent << "vtkIzarAddFieldData\n";
	os << indent << "List of FieldData to add\n";
	vtkIndent nextIndent = indent.GetNextIndent();
	os << indent << this->intDataNames.size() << " int values:\n";
	for(int ii = 0; ii < this->intDataNames.size(); ii++)
	{
		os << nextIndent << this->intDataNames[ii] << " -> " << this->intDataValues[ii] << "\n";
	}
	os << indent << this->doubleDataNames.size() << " double values:\n";
	for(int ii = 0; ii < this->doubleDataNames[ii].size(); ii++)
	{
		os << nextIndent << this->doubleDataNames[ii] << " -> " << this->doubleDataValues[ii] << "\n";
	}
}

void vtkIzarAddFieldData::AddIntData(const std::string& name, int data)
{
	this->intDataNames.push_back(name);
	this->intDataValues.push_back(data);
}

void vtkIzarAddFieldData::AddDoubleData(const std::string& name, double data)
{
	this->doubleDataNames.push_back(name);
	this->doubleDataValues.push_back(data);
}

void vtkIzarAddFieldData::clearData()
{
	this->intDataNames.clear();
	this->intDataValues.clear();
	this->doubleDataNames.clear();
	this->doubleDataValues.clear();
}

int vtkIzarAddFieldData::RequestData(vtkInformation* request, vtkInformationVector** inVector, vtkInformationVector* outVector)
{
	vtkInformation* inInfo = inVector[0]->GetInformationObject(0);
	if(!inInfo)
	{
		return 0;
	}
	vtkDataSet* inData = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!inData)
	{
		return 0;
	}
	
	vtkInformation* outInfo = outVector->GetInformationObject(0);
	if(!outInfo)
	{
		return 0;
	}
	vtkDataSet* outData = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	if(!outData)
	{
		return 0;
	}
	outData->ShallowCopy(inData);
	for(int ii = 0; ii < this->intDataNames.size(); ii++)
	{
		vtkIntArray* ar = vtkIntArray::New();
		ar->SetNumberOfComponents(1);
		ar->SetNumberOfTuples(1);
		int* data = reinterpret_cast<int*>(ar->GetVoidPointer(0));
		data[0] = this->intDataValues[ii];
		ar->SetName(this->intDataNames[ii].c_str());
		outData->GetFieldData()->AddArray(ar);
		ar->Delete();
	}
	
	for(int ii = 0; ii < this->doubleDataNames.size(); ii++)
	{
		vtkDoubleArray* ar = vtkDoubleArray::New();
		ar->SetNumberOfComponents(1);
		ar->SetNumberOfTuples(1);
		double* data = reinterpret_cast<double*>(ar->GetVoidPointer(0));
		data[0] = this->doubleDataValues[ii];
		ar->SetName(this->doubleDataNames[ii].c_str());
		outData->GetFieldData()->AddArray(ar);
		ar->Delete();
	}
	return 1;
}
